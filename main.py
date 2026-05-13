#!/usr/bin/env python3
"""
Check plant pathogen names against the Evo2 eukaryotic dataset and full Evo2 list.

RUN: python main.py
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
DEFAULT_COMBINED_PATHOGEN_FILE = Path("combined_plant_pathogen_list.txt")
DEFAULT_EUKARYOTIC_PATHOGEN_FILE = Path("evo2_eukaryotic_dataset.txt")
DEFAULT_EVO2_TRAINING_FILE = Path("evo2_full_training_dataset.txt")
DEFAULT_EUKARYOTIC_OUTPUT_FILE = Path("plant_pathogens_vs._eukaryotes_evo2")
DEFAULT_ENTIRE_EVO2_OUTPUT_FILE = Path("plant_pathogen_vs._entire_evo2")

# Names from the UC IPM table that are not actual organisms.
SKIP_NAMES = {"", "none", "unknown", "various"}

# Words after these are usually strain/pathovar details.
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
    """One row from an assembly/species file."""

    assembly_id: str
    species_name: str


class SimpleTableParser(HTMLParser):
    """Small table parser for the UC IPM page."""

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
    """Clean up repeated whitespace."""
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def normalize_name(name):
    """Lowercase and remove punctuation for matching."""
    name = html.unescape(name).lower()
    name = name.replace("&", " and ")
    name = re.sub(r"\b(spp?|sp)\.", r"\1", name)
    name = re.sub(r"[^a-z0-9]+", " ", name)
    return clean_spaces(name)


def remove_parentheses(name):
    """Drop notes in parentheses."""
    return clean_spaces(re.sub(r"\([^)]*\)", " ", name))


def possible_name_versions(name):
    """Handle a few simple synonym formats."""
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
    """Keep the genus/species part used for matching."""
    words = normalize_name(remove_parentheses(name)).split()
    if not words:
        return ""

    # Candidatus names usually need the third word too.
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
    """True for names like Pseudomonas spp."""
    normalized = normalize_name(name)
    return bool(re.search(r"\b(spp|sp)\b", normalized)) and len(normalized.split()) >= 2


def first_word(name):
    """Get the genus-ish first word."""
    words = normalize_name(name).split()
    return words[0] if words else ""


def read_evo2_file(path):
    """Read an assembly/species file."""
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        sample = file.read(4096)
        file.seek(0)

        delimiter = "\t" if "\t" in sample else ","
        reader = csv.DictReader(file, delimiter=delimiter)

        # Tab/csv files with headers.
        if reader.fieldnames and "Species_Name" in reader.fieldnames:
            return [
                Evo2Organism(
                    assembly_id=clean_spaces(row.get("Assembly_ID", "")),
                    species_name=clean_spaces(row["Species_Name"]),
                )
                for row in reader
                if row.get("Species_Name")
            ]

        # Plain one-name-per-line files.
        file.seek(0)
        return [
            Evo2Organism(assembly_id=clean_spaces(line), species_name=clean_spaces(line))
            for line in file
            if clean_spaces(line) and not line.lower().startswith("assembly_id")
        ]


def read_pathogen_name_file(path):
    """Read a plain list of pathogen names."""
    with path.open("r", encoding="utf-8-sig") as file:
        return [
            clean_spaces(line)
            for line in file
            if clean_spaces(line) and not clean_spaces(line).startswith("#")
        ]


def download_uc_ipm_page(url):
    """Fetch the UC IPM disease list."""
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "GenomeModeling Evo2 pathogen checker"},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8", errors="replace")


def get_pathogens_from_uc_ipm(page_html):
    """Pull scientific names out of the UC IPM table."""
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
    """Build lookup tables for matching."""
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
    """Find Evo2 assembly IDs for one name."""
    matches = set()

    for version in sorted(possible_name_versions(pathogen_name)):
        normalized_version = normalize_name(version)

        # Exact match.
        if normalized_version in normalized_evo2_names:
            matches.update(normalized_evo2_names[normalized_version])

        # Genus-level names, like Pseudomonas spp.
        if is_genus_level_name(version):
            genus = first_word(version)
            if genus in evo2_genera:
                matches.update(evo2_genera[genus])

        # Genus + species match.
        main_name = main_part_of_name(version)
        if main_name in main_evo2_names:
            matches.update(main_evo2_names[main_name])

        # Match names that have strain details after the species.
        for evo2_name, assembly_ids in normalized_evo2_names.items():
            if evo2_name.startswith(f"{normalized_version} "):
                matches.update(assembly_ids)

    return matches


def write_output_file(matches, output_path, source_name):
    """Write the Y/N result file."""
    with output_path.open("w", encoding="utf-8", newline="\n") as file:
        file.write(f"# Y = found in {source_name}; N = not found.\n")
        file.write("# Assembly_IDs contains matching Evo2 assembly IDs when Y; N rows leave it blank.\n")
        file.write("# Assembly_IDs\tPathogen\tY/N\n\n\n")

        for pathogen in sorted(matches):
            assembly_ids = matches[pathogen]
            status = "Y" if assembly_ids else "N"
            assembly_id_text = ";".join(sorted(assembly_id for assembly_id in assembly_ids if assembly_id))
            file.write(f"{assembly_id_text}\t{pathogen}\t{status}\n")


def get_command_line_args():
    """Command-line options."""
    parser = argparse.ArgumentParser(
        description="Check whether UC IPM plant pathogens were used in Evo2 training."
    )
    parser.add_argument(
        "--combined-pathogens",
        default=DEFAULT_COMBINED_PATHOGEN_FILE,
        type=Path,
        help=f"Path to the extra plant pathogen list. Default: {DEFAULT_COMBINED_PATHOGEN_FILE}",
    )
    parser.add_argument(
        "--eukaryotic-dataset",
        default=DEFAULT_EUKARYOTIC_PATHOGEN_FILE,
        type=Path,
        help=(
            "Path to the Evo2 eukaryotic dataset. "
            f"Default: {DEFAULT_EUKARYOTIC_PATHOGEN_FILE}"
        ),
    )
    parser.add_argument(
        "--evo2-training-file",
        default=DEFAULT_EVO2_TRAINING_FILE,
        type=Path,
        help=f"Path to the full Evo2 training dataset. Default: {DEFAULT_EVO2_TRAINING_FILE}",
    )
    parser.add_argument(
        "--eukaryotic-output",
        default=DEFAULT_EUKARYOTIC_OUTPUT_FILE,
        type=Path,
        help=f"Path for the eukaryotic comparison output. Default: {DEFAULT_EUKARYOTIC_OUTPUT_FILE}",
    )
    parser.add_argument(
        "--entire-evo2-output",
        default=DEFAULT_ENTIRE_EVO2_OUTPUT_FILE,
        type=Path,
        help=f"Path for the full Evo2 comparison output. Default: {DEFAULT_ENTIRE_EVO2_OUTPUT_FILE}",
    )
    parser.add_argument("--url", default=UC_IPM_URL, help="UC IPM disease-list URL.")
    return parser.parse_args()


def combined_pathogen_names(url, combined_pathogens_path):
    """Website names plus the local pathogen list."""
    page_html = download_uc_ipm_page(url)
    pathogens = get_pathogens_from_uc_ipm(page_html)
    if not pathogens:
        raise RuntimeError("no pathogen names were found on the UC IPM page")

    extra_pathogens = read_pathogen_name_file(combined_pathogens_path)
    return list(dict.fromkeys(pathogens + extra_pathogens))


def write_matches_for_names(pathogen_names, lookup_maps, output_path, source_name):
    """Make one result file."""
    assembly_names, assembly_main_names, assembly_genera = lookup_maps
    matches = {}
    for pathogen in pathogen_names:
        if pathogen not in matches:
            matches[pathogen] = matching_assembly_ids(
                pathogen,
                assembly_names,
                assembly_main_names,
                assembly_genera,
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_output_file(matches, output_path, source_name)

    yes_count = sum(bool(assembly_ids) for assembly_ids in matches.values())
    no_count = sum(not assembly_ids for assembly_ids in matches.values())
    print(f"Wrote {len(matches)} pathogen results to {output_path}")
    print(f"Y: {yes_count}")
    print(f"N: {no_count}")
    return matches


def main():
    args = get_command_line_args()

    if not args.combined_pathogens.exists():
        print(f"Error: combined pathogen file does not exist: {args.combined_pathogens}", file=sys.stderr)
        return 2
    if not args.eukaryotic_dataset.exists():
        print(f"Error: eukaryotic dataset file does not exist: {args.eukaryotic_dataset}", file=sys.stderr)
        return 2
    if not args.evo2_training_file.exists():
        print(f"Error: Evo2 training file does not exist: {args.evo2_training_file}", file=sys.stderr)
        return 2

    eukaryotic_rows = read_evo2_file(args.eukaryotic_dataset)
    eukaryotic_lookup_maps = make_training_lookup_maps(eukaryotic_rows)

    evo2_rows = read_evo2_file(args.evo2_training_file)
    evo2_lookup_maps = make_training_lookup_maps(evo2_rows)

    try:
        pathogen_names = combined_pathogen_names(
            args.url,
            args.combined_pathogens,
        )
    except RuntimeError as error:
        print(f"Error: {error}.", file=sys.stderr)
        return 1

    write_matches_for_names(
        pathogen_names,
        eukaryotic_lookup_maps,
        args.eukaryotic_output,
        args.eukaryotic_dataset.name,
    )

    write_matches_for_names(
        pathogen_names,
        evo2_lookup_maps,
        args.entire_evo2_output,
        args.evo2_training_file.name,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
