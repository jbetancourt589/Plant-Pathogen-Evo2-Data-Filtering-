#!/usr/bin/env python3
"""Count non-ACGT letters in the Evo2 Data Reproduction 3 OpenGenome2 FASTA."""

from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
FASTA_PATH = (
    PROJECT_ROOT
    / "Datasets"
    / "Evo2 Data Reproduction 3"
    / "open_genome2_filtered_GCA_001599115.1.fasta"
)
OUTPUT_PATH = (
    PROJECT_ROOT
    / "Results"
    / "Evo2 Data Reproduction #3"
    / "non_acgt_counts.txt"
)
ACGT = {"A", "C", "G", "T"}


def count_non_acgt(path: Path) -> Counter[str]:
    """Return counts of non-ACGT sequence characters in a FASTA file."""
    counts: Counter[str] = Counter()

    with path.open("r", encoding="utf-8") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line or line.startswith(">"):
                continue

            for base in line:
                base_upper = base.upper()
                if base_upper not in ACGT:
                    counts[base_upper] += 1

    return counts


def main() -> None:
    """Print non-ACGT character counts."""
    if not FASTA_PATH.is_file():
        raise FileNotFoundError(f"FASTA file was not found: {FASTA_PATH}")

    counts = count_non_acgt(FASTA_PATH)
    total = sum(counts.values())
    output_lines = [
        f"File: {FASTA_PATH}",
        f"Total non-ACGT letters: {total}",
    ]
    output_lines.extend(f"{base}: {count}" for base, count in sorted(counts.items()))

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text("\n".join(output_lines) + "\n", encoding="utf-8")

    print(f"Wrote counts: {OUTPUT_PATH}")
    print("\n".join(output_lines))


if __name__ == "__main__":
    main()
