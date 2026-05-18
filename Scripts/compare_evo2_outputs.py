#!/usr/bin/env python3
"""
Compare John-filtered and OpenGenome/Evo2-filtered FASTA files.

This script uses only the Python standard library. It:
- summarizes each FASTA file
- compares John output to OpenGenome output
- validates OpenGenome split/chunk records against the original FASTA
- recreates an OpenGenome-like A/C/G/T-only chunk filter
- compares the recreated FASTA to the actual OpenGenome FASTA
"""

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path


MIN_CHUNK_LENGTH = 10_000
ACGT = {"A", "C", "G", "T"}
SPLIT_ID_RE = re.compile(r"^(.+):(\d+)-(\d+)$")
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ORIGINAL_PATH = PROJECT_ROOT / "Datasets" / "orginial_from_ncbi_GCA_000359685.2.fasta"
DEFAULT_JOHN_PATH = PROJECT_ROOT / "Datasets" / "john_filtered_GCA_000359685.2.fasta"
DEFAULT_OPENGENOME_PATH = PROJECT_ROOT / "Datasets" / "open_genome2_filtered_GCA_000359685.2.fasta"
DEFAULT_OUTDIR = PROJECT_ROOT / "Results" / "Evo2 Data Reproduction"


@dataclass
class FastaRecord:
    """One FASTA record."""

    record_id: str
    header: str
    sequence: str


@dataclass
class FastaStats:
    """Basic summary statistics for a FASTA file."""

    record_count: int
    total_bp: int
    total_n_count: int
    total_non_acgt_count: int
    min_contig_length: int
    max_contig_length: int
    mean_contig_length: float


@dataclass
class SplitId:
    """Coordinates parsed from an OpenGenome-style split record ID."""

    original_contig: str
    start: int
    end: int


def parse_fasta(path: Path) -> dict[str, FastaRecord]:
    """Read a FASTA file into a dictionary keyed by the first header token."""
    records: dict[str, FastaRecord] = {}
    current_header = None
    current_id = None
    sequence_lines: list[str] = []

    with path.open("r", encoding="utf-8") as fasta_file:
        for line_number, raw_line in enumerate(fasta_file, start=1):
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None and current_id is not None:
                    add_fasta_record(records, current_id, current_header, sequence_lines)

                current_header = line[1:].strip()
                if not current_header:
                    raise ValueError(f"Empty FASTA header in {path} on line {line_number}")

                current_id = current_header.split()[0]
                sequence_lines = []
            else:
                if current_header is None:
                    raise ValueError(f"Sequence found before first header in {path} on line {line_number}")

                sequence_lines.append(line)

    if current_header is not None and current_id is not None:
        add_fasta_record(records, current_id, current_header, sequence_lines)

    return records


def add_fasta_record(
    records: dict[str, FastaRecord],
    record_id: str,
    header: str,
    sequence_lines: list[str],
) -> None:
    """Add one parsed FASTA record and reject duplicate IDs."""
    if record_id in records:
        raise ValueError(f"Duplicate FASTA record ID found: {record_id}")

    sequence = "".join(sequence_lines).upper()
    records[record_id] = FastaRecord(record_id=record_id, header=header, sequence=sequence)


def calculate_stats(records: dict[str, FastaRecord]) -> FastaStats:
    """Calculate basic FASTA statistics."""
    lengths = [len(record.sequence) for record in records.values()]
    total_bp = sum(lengths)
    total_n_count = sum(record.sequence.count("N") for record in records.values())
    total_non_acgt_count = sum(
        count_non_acgt(record.sequence) for record in records.values()
    )

    if lengths:
        min_length = min(lengths)
        max_length = max(lengths)
        mean_length = total_bp / len(lengths)
    else:
        min_length = 0
        max_length = 0
        mean_length = 0.0

    return FastaStats(
        record_count=len(records),
        total_bp=total_bp,
        total_n_count=total_n_count,
        total_non_acgt_count=total_non_acgt_count,
        min_contig_length=min_length,
        max_contig_length=max_length,
        mean_contig_length=mean_length,
    )


def count_non_acgt(sequence: str) -> int:
    """Count characters that are not A, C, G, or T. This includes N."""
    return sum(1 for base in sequence.upper() if base not in ACGT)


def parse_split_id(record_id: str) -> SplitId | None:
    """Parse IDs like AOUM02000103.1:2471-24374."""
    match = SPLIT_ID_RE.match(record_id)
    if not match:
        return None

    return SplitId(
        original_contig=match.group(1),
        start=int(match.group(2)),
        end=int(match.group(3)),
    )


def compare_record_sets(
    left_records: dict[str, FastaRecord],
    right_records: dict[str, FastaRecord],
) -> dict[str, object]:
    """Compare two FASTA dictionaries by record ID and sequence."""
    left_ids = set(left_records)
    right_ids = set(right_records)
    shared_ids = left_ids & right_ids

    matching_sequences = sorted(
        record_id
        for record_id in shared_ids
        if left_records[record_id].sequence == right_records[record_id].sequence
    )
    sequence_mismatches = sorted(
        record_id
        for record_id in shared_ids
        if left_records[record_id].sequence != right_records[record_id].sequence
    )

    return {
        "only_left": sorted(left_ids - right_ids),
        "only_right": sorted(right_ids - left_ids),
        "shared_ids": sorted(shared_ids),
        "matching_sequences": matching_sequences,
        "sequence_mismatches": sequence_mismatches,
        "left_total_bp": sum(len(record.sequence) for record in left_records.values()),
        "right_total_bp": sum(len(record.sequence) for record in right_records.values()),
    }


def validate_opengenome_chunks(
    original_records: dict[str, FastaRecord],
    opengenome_records: dict[str, FastaRecord],
) -> list[dict[str, object]]:
    """Validate OpenGenome-style chunk records against the original FASTA."""
    validation_rows: list[dict[str, object]] = []

    for record_id, record in sorted(opengenome_records.items()):
        split_id = parse_split_id(record_id)
        if split_id is None:
            continue

        original_record = original_records.get(split_id.original_contig)
        chunk_length = len(record.sequence)
        contains_n = "N" in record.sequence
        contains_non_acgt = count_non_acgt(record.sequence) > 0

        original_found = original_record is not None
        coordinate_match = check_coordinate_match(original_record, split_id, record.sequence)
        expected_length = coordinate_match["expected_length"]
        coordinates_valid = coordinate_match["coordinates_valid"]
        sequence_matches_original = coordinate_match["sequence_matches_original"]
        coordinate_convention = coordinate_match["coordinate_convention"]
        note = ""

        if not original_found:
            note = "original contig not found"
        elif not coordinates_valid:
            note = "coordinates outside original contig"
        elif expected_length != chunk_length:
            note = "chunk length does not match coordinates"
        elif not sequence_matches_original:
            note = "sequence differs from original subsequence"

        validation_rows.append(
            {
                "record_id": record_id,
                "original_contig": split_id.original_contig,
                "start": split_id.start,
                "end": split_id.end,
                "chunk_length": chunk_length,
                "expected_length_from_coordinates": expected_length,
                "original_contig_found": original_found,
                "coordinates_valid": coordinates_valid,
                "matching_coordinate_convention": coordinate_convention,
                "sequence_matches_original": sequence_matches_original,
                "contains_n": contains_n,
                "contains_non_acgt": contains_non_acgt,
                "note": note,
            }
        )

    return validation_rows


def check_coordinate_match(
    original_record: FastaRecord | None,
    split_id: SplitId,
    chunk_sequence: str,
) -> dict[str, object]:
    """Check whether parsed coordinates match the original under common conventions."""
    if original_record is None:
        return {
            "coordinates_valid": False,
            "expected_length": split_id.end - split_id.start + 1,
            "sequence_matches_original": False,
            "coordinate_convention": "",
        }

    one_based = coordinate_candidate(
        original_record.sequence,
        split_id.start,
        split_id.end,
        chunk_sequence,
        "1-based inclusive",
    )
    zero_based = coordinate_candidate(
        original_record.sequence,
        split_id.start,
        split_id.end,
        chunk_sequence,
        "0-based half-open",
    )

    for candidate in [one_based, zero_based]:
        if candidate["sequence_matches_original"]:
            return candidate

    if one_based["coordinates_valid"]:
        return one_based
    if zero_based["coordinates_valid"]:
        return zero_based
    return one_based


def coordinate_candidate(
    original_sequence: str,
    start: int,
    end: int,
    chunk_sequence: str,
    convention: str,
) -> dict[str, object]:
    """Build one coordinate-convention validation candidate."""
    if convention == "1-based inclusive":
        coordinates_valid = start >= 1 and end >= start and end <= len(original_sequence)
        expected_length = end - start + 1
        original_subsequence = original_sequence[start - 1 : end] if coordinates_valid else ""
    else:
        coordinates_valid = start >= 0 and end >= start and end <= len(original_sequence)
        expected_length = end - start
        original_subsequence = original_sequence[start:end] if coordinates_valid else ""

    return {
        "coordinates_valid": coordinates_valid,
        "expected_length": expected_length,
        "sequence_matches_original": coordinates_valid and chunk_sequence == original_subsequence,
        "coordinate_convention": convention if coordinates_valid else "",
    }


def recreate_opengenome_like_records(
    original_records: dict[str, FastaRecord],
    min_chunk_length: int,
) -> dict[str, FastaRecord]:
    """Split original records into continuous A/C/G/T-only chunks."""
    recreated_records: dict[str, FastaRecord] = {}

    for original_id, original_record in original_records.items():
        sequence = original_record.sequence
        chunk_start_index = None

        for index, base in enumerate(sequence):
            if base in ACGT:
                if chunk_start_index is None:
                    chunk_start_index = index
            elif chunk_start_index is not None:
                add_recreated_chunk(
                    recreated_records,
                    original_id,
                    sequence,
                    chunk_start_index,
                    index - 1,
                    min_chunk_length,
                )
                chunk_start_index = None

        if chunk_start_index is not None:
            add_recreated_chunk(
                recreated_records,
                original_id,
                sequence,
                chunk_start_index,
                len(sequence) - 1,
                min_chunk_length,
            )

    return recreated_records


def add_recreated_chunk(
    recreated_records: dict[str, FastaRecord],
    original_id: str,
    sequence: str,
    start_index: int,
    end_index: int,
    min_chunk_length: int,
) -> None:
    """Add one recreated chunk if it is long enough."""
    chunk_length = end_index - start_index + 1
    if chunk_length < min_chunk_length:
        return

    start_coordinate = start_index + 1
    end_coordinate = end_index + 1
    record_id = f"{original_id}:{start_coordinate}-{end_coordinate}"
    chunk_sequence = sequence[start_index : end_index + 1]
    recreated_records[record_id] = FastaRecord(
        record_id=record_id,
        header=record_id,
        sequence=chunk_sequence,
    )


def write_fasta(records: dict[str, FastaRecord], path: Path) -> None:
    """Write records to a FASTA file with 80 bases per line."""
    with path.open("w", encoding="utf-8", newline="\n") as output_file:
        for record in records.values():
            output_file.write(f">{record.header}\n")
            for start in range(0, len(record.sequence), 80):
                output_file.write(record.sequence[start : start + 80] + "\n")


def write_comparison_report(
    john_records: dict[str, FastaRecord],
    opengenome_records: dict[str, FastaRecord],
    recreated_records: dict[str, FastaRecord],
    path: Path,
) -> None:
    """Write record-level comparisons to CSV."""
    with path.open("w", encoding="utf-8", newline="") as csv_file:
        fieldnames = [
            "comparison",
            "record_id",
            "status",
            "left_name",
            "left_length",
            "right_name",
            "right_length",
            "sequences_match",
            "note",
        ]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        write_pairwise_rows(
            writer,
            "john_vs_opengenome",
            "john",
            john_records,
            "opengenome",
            opengenome_records,
        )
        write_pairwise_rows(
            writer,
            "recreated_vs_opengenome",
            "recreated",
            recreated_records,
            "opengenome",
            opengenome_records,
        )


def write_pairwise_rows(
    writer: csv.DictWriter,
    comparison_name: str,
    left_name: str,
    left_records: dict[str, FastaRecord],
    right_name: str,
    right_records: dict[str, FastaRecord],
) -> None:
    """Write rows for one pairwise FASTA comparison."""
    all_ids = sorted(set(left_records) | set(right_records))

    for record_id in all_ids:
        left_record = left_records.get(record_id)
        right_record = right_records.get(record_id)

        if left_record is None:
            status = f"only_in_{right_name}"
            sequences_match = ""
            note = f"missing from {left_name}"
        elif right_record is None:
            status = f"only_in_{left_name}"
            sequences_match = ""
            note = f"missing from {right_name}"
        else:
            sequences_match = left_record.sequence == right_record.sequence
            status = "shared_id_matching_sequence" if sequences_match else "shared_id_sequence_mismatch"
            note = ""

        writer.writerow(
            {
                "comparison": comparison_name,
                "record_id": record_id,
                "status": status,
                "left_name": left_name,
                "left_length": len(left_record.sequence) if left_record else "",
                "right_name": right_name,
                "right_length": len(right_record.sequence) if right_record else "",
                "sequences_match": sequences_match,
                "note": note,
            }
        )


def write_chunk_validation_csv(validation_rows: list[dict[str, object]], path: Path) -> None:
    """Write OpenGenome chunk validation rows to CSV."""
    fieldnames = [
        "record_id",
        "original_contig",
        "start",
        "end",
        "chunk_length",
        "expected_length_from_coordinates",
        "original_contig_found",
        "coordinates_valid",
        "matching_coordinate_convention",
        "sequence_matches_original",
        "contains_n",
        "contains_non_acgt",
        "note",
    ]

    with path.open("w", encoding="utf-8", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(validation_rows)


def write_summary(
    original_path: Path,
    john_path: Path,
    opengenome_path: Path,
    original_stats: FastaStats,
    john_stats: FastaStats,
    opengenome_stats: FastaStats,
    recreated_stats: FastaStats,
    john_vs_opengenome: dict[str, object],
    recreated_vs_opengenome: dict[str, object],
    validation_rows: list[dict[str, object]],
    path: Path,
) -> None:
    """Write a plain-text summary of the analysis."""
    split_records = len(validation_rows)
    split_records_matching_original = sum(
        1 for row in validation_rows if row["sequence_matches_original"]
    )
    split_records_with_non_acgt = sum(
        1 for row in validation_rows if row["contains_non_acgt"]
    )
    coordinate_conventions = summarize_coordinate_conventions(validation_rows)

    recreated_missing = len(recreated_vs_opengenome["only_right"])
    recreated_extra = len(recreated_vs_opengenome["only_left"])
    recreated_mismatches = len(recreated_vs_opengenome["sequence_mismatches"])
    recreated_matches_actual = (
        recreated_missing == 0 and recreated_extra == 0 and recreated_mismatches == 0
    )

    opengenome_appears_split = (
        split_records > 0
        and split_records_matching_original == split_records
        and split_records_with_non_acgt == 0
    )

    with path.open("w", encoding="utf-8", newline="\n") as summary_file:
        summary_file.write("Evo2/OpenGenome Filtering Reproducibility Summary\n")
        summary_file.write("================================================\n\n")

        summary_file.write("Input files\n")
        summary_file.write(f"- Original NCBI FASTA: {original_path}\n")
        summary_file.write(f"- John's filtered FASTA: {john_path}\n")
        summary_file.write(f"- OpenGenome/Evo2 filtered FASTA: {opengenome_path}\n\n")

        summary_file.write("FASTA summary statistics\n")
        write_stats_line(summary_file, "Original NCBI", original_stats)
        write_stats_line(summary_file, "John filtered", john_stats)
        write_stats_line(summary_file, "OpenGenome/Evo2 filtered", opengenome_stats)
        write_stats_line(summary_file, "Recreated OpenGenome-like", recreated_stats)
        summary_file.write("\n")

        summary_file.write("John vs OpenGenome/Evo2\n")
        summary_file.write(f"- Records only in John: {len(john_vs_opengenome['only_left'])}\n")
        summary_file.write(f"- Records only in OpenGenome/Evo2: {len(john_vs_opengenome['only_right'])}\n")
        summary_file.write(f"- Records present in both: {len(john_vs_opengenome['shared_ids'])}\n")
        summary_file.write(
            "- Total bp difference (John minus OpenGenome/Evo2): "
            f"{john_vs_opengenome['left_total_bp'] - john_vs_opengenome['right_total_bp']}\n\n"
        )

        summary_file.write("N and non-ACGT content\n")
        summary_file.write(boolean_line("John's output still contains N bases", john_stats.total_n_count > 0))
        summary_file.write(f"  N count: {john_stats.total_n_count}\n")
        summary_file.write(boolean_line("OpenGenome/Evo2 contains N bases", opengenome_stats.total_n_count > 0))
        summary_file.write(f"  N count: {opengenome_stats.total_n_count}\n")
        summary_file.write(
            f"- John non-ACGT count, including N: {john_stats.total_non_acgt_count}\n"
        )
        summary_file.write(
            "- OpenGenome/Evo2 non-ACGT count, including N: "
            f"{opengenome_stats.total_non_acgt_count}\n\n"
        )

        summary_file.write("OpenGenome/Evo2 chunk validation\n")
        summary_file.write(f"- OpenGenome-style split records detected: {split_records}\n")
        summary_file.write(
            "- Split records matching the original NCBI subsequence exactly: "
            f"{split_records_matching_original} of {split_records}\n"
        )
        summary_file.write(
            f"- Split records containing non-ACGT characters: {split_records_with_non_acgt}\n"
        )
        summary_file.write(f"- Coordinate convention observed: {coordinate_conventions}\n")
        summary_file.write(
            boolean_line(
                "OpenGenome/Evo2 appears to split contigs around N/non-ACGT regions",
                opengenome_appears_split,
            )
        )
        summary_file.write("\n")

        summary_file.write("Recreated OpenGenome-like comparison\n")
        summary_file.write(
            f"- Exact matching record IDs: {len(recreated_vs_opengenome['shared_ids'])}\n"
        )
        summary_file.write(
            "- Exact matching sequences among matching IDs: "
            f"{len(recreated_vs_opengenome['matching_sequences'])}\n"
        )
        summary_file.write(f"- Records missing from recreated output: {recreated_missing}\n")
        summary_file.write(f"- Extra records in recreated output: {recreated_extra}\n")
        summary_file.write(
            f"- Sequence mismatches where IDs match: {recreated_mismatches}\n"
        )
        summary_file.write(
            boolean_line(
                "Recreated OpenGenome-like filtering matches the actual OpenGenome/Evo2 file",
                recreated_matches_actual,
            )
        )
        summary_file.write(
            "- Note: recreated chunk names use the requested 1-based inclusive coordinates. "
            "If the actual OpenGenome/Evo2 file uses a different coordinate convention, "
            "record IDs can differ even when chunk sequences and total bp agree.\n"
        )


def write_stats_line(summary_file, label: str, stats: FastaStats) -> None:
    """Write one FASTA statistics line to the summary file."""
    summary_file.write(
        f"- {label}: records={stats.record_count}, total_bp={stats.total_bp}, "
        f"N_count={stats.total_n_count}, non_ACGT_count={stats.total_non_acgt_count}, "
        f"min_len={stats.min_contig_length}, max_len={stats.max_contig_length}, "
        f"mean_len={stats.mean_contig_length:.2f}\n"
    )


def boolean_line(label: str, value: bool) -> str:
    """Format a yes/no line for the text summary."""
    return f"- {label}: {'YES' if value else 'NO'}\n"


def summarize_coordinate_conventions(validation_rows: list[dict[str, object]]) -> str:
    """Summarize coordinate conventions seen in validated split records."""
    if not validation_rows:
        return "none"

    counts: dict[str, int] = {}
    for row in validation_rows:
        convention = str(row["matching_coordinate_convention"] or "unmatched")
        counts[convention] = counts.get(convention, 0) + 1

    return ", ".join(f"{name}={count}" for name, count in sorted(counts.items()))


def build_arg_parser() -> argparse.ArgumentParser:
    """Create the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Compare John and OpenGenome/Evo2 filtered FASTA outputs."
    )
    parser.add_argument(
        "--original",
        type=Path,
        default=DEFAULT_ORIGINAL_PATH,
        help=f"Original NCBI input FASTA. Default: {DEFAULT_ORIGINAL_PATH}",
    )
    parser.add_argument(
        "--john",
        type=Path,
        default=DEFAULT_JOHN_PATH,
        help=f"John's filtered output FASTA. Default: {DEFAULT_JOHN_PATH}",
    )
    parser.add_argument(
        "--opengenome",
        type=Path,
        default=DEFAULT_OPENGENOME_PATH,
        help=f"OpenGenome/Evo2 filtered output FASTA. Default: {DEFAULT_OPENGENOME_PATH}",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help=f"Directory where analysis outputs will be written. Default: {DEFAULT_OUTDIR}",
    )
    return parser


def require_existing_file(path: Path, label: str) -> None:
    """Give a clear error if an expected input FASTA file is missing."""
    if not path.exists():
        raise FileNotFoundError(f"{label} file was not found: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"{label} path is not a file: {path}")


def main() -> None:
    """Run the FASTA comparison workflow."""
    args = build_arg_parser().parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    require_existing_file(args.original, "Original NCBI FASTA")
    require_existing_file(args.john, "John filtered FASTA")
    require_existing_file(args.opengenome, "OpenGenome/Evo2 filtered FASTA")

    original_records = parse_fasta(args.original)
    john_records = parse_fasta(args.john)
    opengenome_records = parse_fasta(args.opengenome)

    original_stats = calculate_stats(original_records)
    john_stats = calculate_stats(john_records)
    opengenome_stats = calculate_stats(opengenome_records)

    john_vs_opengenome = compare_record_sets(john_records, opengenome_records)
    validation_rows = validate_opengenome_chunks(original_records, opengenome_records)

    recreated_records = recreate_opengenome_like_records(
        original_records,
        min_chunk_length=MIN_CHUNK_LENGTH,
    )
    recreated_stats = calculate_stats(recreated_records)
    recreated_vs_opengenome = compare_record_sets(recreated_records, opengenome_records)

    recreated_fasta_path = outdir / "recreated_opengenome_like.fasta"
    summary_path = outdir / "summary.txt"
    comparison_report_path = outdir / "comparison_report.csv"
    chunk_validation_path = outdir / "opengenome_chunk_validation.csv"

    write_fasta(recreated_records, recreated_fasta_path)
    write_comparison_report(
        john_records,
        opengenome_records,
        recreated_records,
        comparison_report_path,
    )
    write_chunk_validation_csv(validation_rows, chunk_validation_path)
    write_summary(
        args.original,
        args.john,
        args.opengenome,
        original_stats,
        john_stats,
        opengenome_stats,
        recreated_stats,
        john_vs_opengenome,
        recreated_vs_opengenome,
        validation_rows,
        summary_path,
    )

    print(f"Wrote summary: {summary_path}")
    print(f"Wrote comparison report: {comparison_report_path}")
    print(f"Wrote chunk validation: {chunk_validation_path}")
    print(f"Wrote recreated FASTA: {recreated_fasta_path}")


if __name__ == "__main__":
    main()
