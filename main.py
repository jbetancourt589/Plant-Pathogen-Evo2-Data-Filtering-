#!/usr/bin/env python3
"""
Checks if each UC IPM pathogen appears in the Evo2 organism file.

RUN: python main.py --evo2-file all_bacterial_species_in_evo2.txt --output results.txt
"""

import argparse
import csv
import html
import re
import sys
import urllib.request
from dataclasses import dataclass
from html.parser import HTMLParser
from pathlib import Path


UC_IPM_URL = "https://ipm.ucanr.edu/PMG/diseases/diseaseslist.html"

# Skip non-organism names.
SKIP_NAMES = {"", "none", "unknown", "various"}

# Stop matching before strain details.
STRAIN_WORDS = {
    "subsp",
    "subspecies",
    "ssp",
    "str",
    "strain",
    "var",
    "pv",
    "serovar",
    "biovar",
    "genomovar",
    "isolate",
    "clone",
}


@dataclass(frozen=True)
class Evo2Organism:
    """One Evo2 row."""

    assembly_id: str
    species_name: str


class SimpleTableParser(HTMLParser):
    """Collect rows from an HTML table."""

    def __init__(self):
        super().__init__()
        self.rows = []
        self.current_row = []
        self.current_cell = []
        self.inside_cell = False

    def handle_starttag(self, tag, attrs):
        if tag.lower() in {"td", "th"}:
            self.inside_cell = True
            self.current_cell = []

    def handle_data(self, data):
        if self.inside_cell:
            self.current_cell.append(data)

    def handle_endtag(self, tag):
        tag = tag.lower()

        if tag in {"td", "th"} and self.inside_cell:
            cell_text = clean_spaces(" ".join(self.current_cell))
            self.current_row.append(cell_text)
            self.inside_cell = False

        if tag == "tr" and self.current_row:
            self.rows.append(self.current_row)
            self.current_row = []


def clean_spaces(text):
    """Collapse extra whitespace."""
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def normalize_name(name):
    """Make names easier to compare."""
    name = html.unescape(name).lower()
    name = name.replace("&", " and ")
    name = re.sub(r"\b(spp?|sp)\.", r"\1", name)
    name = re.sub(r"[^a-z0-9]+", " ", name)
    return clean_spaces(name)


def remove_parentheses(name):
    """Remove parenthetical notes."""
    return clean_spaces(re.sub(r"\([^)]*\)", " ", name))


def possible_name_versions(name):
    """Return synonym/name variants."""
    versions = {name, remove_parentheses(name)}

    synonym_match = re.search(r"^(\S+)\s+\((?:syn\.?|synonym)\s+([^)]+)\)\s+(.+)$", name, re.I)
    if synonym_match:
        original_genus, synonym_text, rest_of_name = synonym_match.groups()
        versions.add(f"{original_genus} {rest_of_name}")

        for synonym_genus in re.split(r"[,;/]|\bor\b|\band\b", synonym_text):
            synonym_genus = clean_spaces(synonym_genus)
            if synonym_genus:
                versions.add(f"{synonym_genus} {rest_of_name}")

    return {clean_spaces(version) for version in versions if clean_spaces(version)}


def main_part_of_name(name):
    """Keep the main species-level name."""
    words = normalize_name(remove_parentheses(name)).split()
    if not words:
        return ""

    # Candidatus names usually need three words.
    if words[0] == "candidatus" and len(words) >= 3:
        return " ".join(words[:3])

    kept_words = []
    for word in words:
        if word in STRAIN_WORDS:
            break
        if re.fullmatch(r"[a-z]*\d+[a-z]*", word):
            break

        kept_words.append(word)
        if len(kept_words) == 2:
            break

    return " ".join(kept_words)


def is_genus_level_name(name):
    """Check for broad names like 'Pseudomonas spp.'."""
    normalized = normalize_name(name)
    return bool(re.search(r"\b(spp|sp)\b", normalized)) and len(normalized.split()) >= 2


def first_word(name):
    """Return the first word."""
    words = normalize_name(name).split()
    return words[0] if words else ""


def read_evo2_file(path):
    """Read Evo2 organism rows."""
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        sample = file.read(4096)
        file.seek(0)

        delimiter = "\t" if "\t" in sample else ","
        reader = csv.DictReader(file, delimiter=delimiter)

        # Normal input format.
        if reader.fieldnames and "Species_Name" in reader.fieldnames:
            return [
                Evo2Organism(
                    assembly_id=clean_spaces(row.get("Assembly_ID", "")),
                    species_name=clean_spaces(row["Species_Name"]),
                )
                for row in reader
                if row.get("Species_Name")
            ]

        # Fallback: one organism name per line.
        file.seek(0)
        return [
            Evo2Organism(assembly_id="", species_name=clean_spaces(line))
            for line in file
            if clean_spaces(line) and not line.lower().startswith("assembly_id")
        ]


def download_uc_ipm_page(url):
    """Download the UC IPM page."""
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "GenomeModeling Evo2 pathogen checker"},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8", errors="replace")


def get_pathogens_from_uc_ipm(page_html):
    """Get pathogen names from the UC IPM table."""
    parser = SimpleTableParser()
    parser.feed(page_html)

    pathogens = []
    for row in parser.rows:
        if len(row) < 4:
            continue

        host, common_name, scientific_name, disease_type = row[:4]

        if row[:4] == ["Plant or crop host", "Common name", "Scientific name", "Type"]:
            continue

        if normalize_name(scientific_name) in SKIP_NAMES:
            continue

        pathogens.append(scientific_name)

    return pathogens


def make_training_lookup_maps(evo2_organisms):
    """Build name-to-assembly lookup maps."""
    normalized_names = {}
    main_names = {}
    genera = {}

    for organism in evo2_organisms:
        species_name = organism.species_name
        assembly_id = organism.assembly_id

        normalized_names.setdefault(normalize_name(species_name), set()).add(assembly_id)

        main_name = main_part_of_name(species_name)
        if main_name:
            main_names.setdefault(main_name, set()).add(assembly_id)

        genus = first_word(species_name)
        if genus:
            genera.setdefault(genus, set()).add(assembly_id)

    return normalized_names, main_names, genera


def matching_assembly_ids(pathogen_name, normalized_evo2_names, main_evo2_names, evo2_genera):
    """Find matching assembly IDs."""
    matches = set()

    for version in sorted(possible_name_versions(pathogen_name)):
        normalized_version = normalize_name(version)

        # Exact name match.
        if normalized_version in normalized_evo2_names:
            matches.update(normalized_evo2_names[normalized_version])

        # Genus-level match.
        if is_genus_level_name(version):
            genus = first_word(version)
            if genus in evo2_genera:
                matches.update(evo2_genera[genus])

        # Species-level match.
        main_name = main_part_of_name(version)
        if main_name in main_evo2_names:
            matches.update(main_evo2_names[main_name])

        # Prefix match for names with extra details.
        for evo2_name, assembly_ids in normalized_evo2_names.items():
            if evo2_name.startswith(f"{normalized_version} "):
                matches.update(assembly_ids)

    return matches


def write_output_file(matches, output_path):
    """Write the output file."""
    with output_path.open("w", encoding="utf-8", newline="\n") as file:
        file.write("# Y = found in the Evo2 organism file; N = not found.\n")
        file.write("# Multiple assembly IDs mean multiple Evo2 assemblies matched the same pathogen.\n")
        file.write("# Pathogen\tY/N\tAssembly_IDs\n\n\n")

        for pathogen in sorted(matches):
            assembly_ids = matches[pathogen]
            status = "Y" if assembly_ids else "N"
            assembly_id_text = ";".join(sorted(assembly_id for assembly_id in assembly_ids if assembly_id))
            file.write(f"{pathogen}\t{status}\t{assembly_id_text}\n")


def get_command_line_args():
    """Read command-line options."""
    parser = argparse.ArgumentParser(
        description="Check whether UC IPM plant pathogens were used in Evo2 training."
    )
    parser.add_argument("--evo2-file", required=True, type=Path, help="Path to the Evo2 organism text file.")
    parser.add_argument("--output", required=True, type=Path, help="Path for the output text file.")
    parser.add_argument("--url", default=UC_IPM_URL, help="UC IPM disease-list URL.")
    return parser.parse_args()


def main():
    args = get_command_line_args()

    if not args.evo2_file.exists():
        print(f"Error: Evo2 file does not exist: {args.evo2_file}", file=sys.stderr)
        return 2

    # Read Evo2 organisms.
    evo2_organisms = read_evo2_file(args.evo2_file)
    normalized_evo2_names, main_evo2_names, evo2_genera = make_training_lookup_maps(evo2_organisms)

    # Get UC IPM pathogens.
    page_html = download_uc_ipm_page(args.url)
    pathogens = get_pathogens_from_uc_ipm(page_html)
    if not pathogens:
        print("Error: no pathogen names were found on the UC IPM page.", file=sys.stderr)
        return 1

    # Match pathogens to Evo2 rows.
    matches = {}
    for pathogen in pathogens:
        if pathogen not in matches:
            matches[pathogen] = matching_assembly_ids(
                pathogen,
                normalized_evo2_names,
                main_evo2_names,
                evo2_genera,
            )

    # Write results.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_output_file(matches, args.output)

    yes_count = sum(bool(assembly_ids) for assembly_ids in matches.values())
    no_count = sum(not assembly_ids for assembly_ids in matches.values())
    print(f"Wrote {len(matches)} pathogen results to {args.output}")
    print(f"Y: {yes_count}")
    print(f"N: {no_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
