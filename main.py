#!/usr/bin/env python3
"""
Make a Y/N text file for the organisms that have/haven't been trained by evo2

RUN: python main.py --evo2-file all_bacterial_species_in_evo2.txt --output results.txt

"""

import argparse
import csv
import html
import re
import sys
import urllib.request
from html.parser import HTMLParser
from pathlib import Path


UC_IPM_URL = "https://ipm.ucanr.edu/PMG/diseases/diseaseslist.html"

# These are not real organism names, so the script skips them.
SKIP_NAMES = {"", "none", "unknown", "various"}

# Words that usually mean "extra strain details start here".
# For matching, we mostly care about genus + species, not the strain ID.
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


class SimpleTableParser(HTMLParser):
    """
    Reads HTML tables and stores each row as a list of cell values.

    The UC IPM page has one big table. Python's standard library does not have
    a table reader, so this tiny parser collects the table cells for us.
    """

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
    """Turn repeated spaces/newlines into one normal space."""
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def normalize_name(name):
    """
    Make names easier to compare.

    Example:
    "Pythium spp." and "pythium spp" become the same text.
    """
    name = html.unescape(name).lower()
    name = name.replace("&", " and ")
    name = re.sub(r"\b(spp?|sp)\.", r"\1", name)
    name = re.sub(r"[^a-z0-9]+", " ", name)
    return clean_spaces(name)


def remove_parentheses(name):
    """Remove notes like '(syn. Cladosporium)' from a name."""
    return clean_spaces(re.sub(r"\([^)]*\)", " ", name))


def possible_name_versions(name):
    """
    Return different versions of one pathogen name.

    This helps with names that include synonyms, like:
    "Fusicladium (syn. Cladosporium) carpophilum"
    """
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
    """
    Keep the main part of a scientific name for matching.

    Most exact species matches only need genus + species:
    "Pseudomonas syringae pv. tomato" becomes "pseudomonas syringae".
    """
    words = normalize_name(remove_parentheses(name)).split()
    if not words:
        return ""

    # "Candidatus Liberibacter asiaticus" needs three words, not two.
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
    """
    True for broad names like 'Pseudomonas spp.'.

    For these, the script checks whether any Evo2 organism starts with that
    genus instead of requiring one exact species.
    """
    normalized = normalize_name(name)
    return bool(re.search(r"\b(spp|sp)\b", normalized)) and len(normalized.split()) >= 2


def first_word(name):
    """Return the genus, which is usually the first word."""
    words = normalize_name(name).split()
    return words[0] if words else ""


def read_evo2_file(path):
    """Read the Evo2 training file and return only the organism names."""
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        sample = file.read(4096)
        file.seek(0)

        delimiter = "\t" if "\t" in sample else ","
        reader = csv.DictReader(file, delimiter=delimiter)

        # Your file has a Species_Name column, so this is the normal path.
        if reader.fieldnames and "Species_Name" in reader.fieldnames:
            return [clean_spaces(row["Species_Name"]) for row in reader if row.get("Species_Name")]

        # Backup path for a simpler file with one organism name per line.
        file.seek(0)
        return [
            clean_spaces(line)
            for line in file
            if clean_spaces(line) and not line.lower().startswith("assembly_id")
        ]


def download_uc_ipm_page(url):
    """Download the UC IPM disease-list page."""
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "GenomeModeling Evo2 pathogen checker"},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8", errors="replace")


def get_pathogens_from_uc_ipm(page_html):
    """Pull the scientific-name column out of the UC IPM table."""
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


def make_training_lookup_sets(evo2_names):
    """
    Make lookup sets so matching is fast.

    normalized_names: full normalized organism names
    main_names: genus + species style names
    genera: just the genus names
    """
    normalized_names = {normalize_name(name) for name in evo2_names}
    main_names = {main_part_of_name(name) for name in evo2_names if main_part_of_name(name)}
    genera = {first_word(name) for name in evo2_names if first_word(name)}
    return normalized_names, main_names, genera


def was_used_in_evo2(pathogen_name, normalized_evo2_names, main_evo2_names, evo2_genera):
    """Return True if the pathogen looks like it appears in the Evo2 file."""
    for version in possible_name_versions(pathogen_name):
        normalized_version = normalize_name(version)

        # Best case: exact normalized name match.
        if normalized_version in normalized_evo2_names:
            return True

        # Broad names like "Erwinia spp." match if Evo2 has any Erwinia.
        if is_genus_level_name(version):
            return first_word(version) in evo2_genera

        # Species-level names can match even when Evo2 has extra strain words.
        if main_part_of_name(version) in main_evo2_names:
            return True

        # Example: pathogen is "Xylella fastidiosa" and Evo2 has
        # "Xylella fastidiosa 9a5c".
        if any(evo2_name.startswith(f"{normalized_version} ") for evo2_name in normalized_evo2_names):
            return True

    return False


def write_output_file(pathogens, statuses, output_path):
    """Write the final shareable text file."""
    with output_path.open("w", encoding="utf-8", newline="\n") as file:
        for pathogen in sorted(statuses):
            file.write(f"{pathogen}\t{statuses[pathogen]}\n")


def get_command_line_args():
    """Set up the command-line options."""
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

    # Step 1: read the Evo2 training organisms from your text file.
    evo2_names = read_evo2_file(args.evo2_file)
    normalized_evo2_names, main_evo2_names, evo2_genera = make_training_lookup_sets(evo2_names)

    # Step 2: download the UC IPM page and get the plant pathogen names.
    page_html = download_uc_ipm_page(args.url)
    pathogens = get_pathogens_from_uc_ipm(page_html)
    if not pathogens:
        print("Error: no pathogen names were found on the UC IPM page.", file=sys.stderr)
        return 1

    # Step 3: compare each pathogen against the Evo2 names.
    statuses = {}
    for pathogen in pathogens:
        if pathogen not in statuses:
            statuses[pathogen] = (
                "Y"
                if was_used_in_evo2(pathogen, normalized_evo2_names, main_evo2_names, evo2_genera)
                else "N"
            )

    # Step 4: write the final txt file.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_output_file(pathogens, statuses, args.output)

    yes_count = sum(status == "Y" for status in statuses.values())
    no_count = sum(status == "N" for status in statuses.values())
    print(f"Wrote {len(statuses)} pathogen results to {args.output}")
    print(f"Y: {yes_count}")
    print(f"N: {no_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
