#!/usr/bin/env python3
"""
Shared Evo2/OpenGenome2 FASTA reproduction workflow.

The recreated FASTA is intentionally strict:
- valid sequence characters are only A, C, G, and T
- N and all other ambiguity letters are break points
- kept chunks must be at least 10,000 bp
- chunk IDs use 0-based half-open coordinates

Optional GTF/GFF annotation filtering removes only conservative, true
centromere-region annotations before the strict sequence filtering runs.
"""

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path


MIN_CHUNK_LENGTH = 10_000
ACGT = {"A", "C", "G", "T"}
SPLIT_ID_RE = re.compile(r"^(.+):(\d+)-(\d+)$")

GENE_LIKE_FEATURES = {
    "gene",
    "transcript",
    "mrna",
    "exon",
    "cds",
    "start_codon",
    "stop_codon",
    "five_prime_utr",
    "three_prime_utr",
    "utr",
}


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


@dataclass
class AnnotationRecord:
    """One parsed GTF/GFF annotation row."""

    seqid: str
    source: str
    feature: str
    start: int
    end: int
    attributes: str
    selected_for_removal: bool
    reason: str
    bp_removed: int


@dataclass
class AnnotationResult:
    """Parsed annotation rows plus merged removal intervals."""

    records: list[AnnotationRecord]
    intervals_by_seqid: dict[str, list[tuple[int, int]]]
    total_records: int
    true_centromere_regions: int
    total_centromere_bp_removed: int


def optional_path(value: str) -> Path | None:
    """Convert a CLI path argument, treating an empty string as omitted."""
    if value is None or str(value).strip() == "":
        return None
    return Path(value)


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

    records[record_id] = FastaRecord(
        record_id=record_id,
        header=header,
        sequence="".join(sequence_lines),
    )


def calculate_stats(records: dict[str, FastaRecord]) -> FastaStats:
    """Calculate basic FASTA statistics."""
    lengths = [len(record.sequence) for record in records.values()]
    total_bp = sum(lengths)
    total_n_count = sum(record.sequence.upper().count("N") for record in records.values())
    total_non_acgt_count = sum(count_non_acgt(record.sequence) for record in records.values())

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
    return sum(1 for base in sequence if base.upper() not in ACGT)


def parse_annotation_file(path: Path, original_records: dict[str, FastaRecord]) -> AnnotationResult:
    """Parse a GTF/GFF file and collect conservative centromere intervals."""
    annotation_records: list[AnnotationRecord] = []
    selected_intervals: dict[str, list[tuple[int, int]]] = {}

    with path.open("r", encoding="utf-8") as annotation_file:
        for line_number, raw_line in enumerate(annotation_file, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            columns = line.split("\t")
            if len(columns) < 9:
                annotation_records.append(
                    AnnotationRecord(
                        seqid="",
                        source="",
                        feature="",
                        start=0,
                        end=0,
                        attributes=line,
                        selected_for_removal=False,
                        reason=f"skipped malformed line {line_number}",
                        bp_removed=0,
                    )
                )
                continue

            seqid, source, feature, start_text, end_text = columns[:5]
            attributes = columns[8]
            selected = False
            reason = "not a true centromere-region annotation"
            bp_removed = 0

            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                start = 0
                end = 0
                reason = f"skipped invalid coordinates on line {line_number}"
            else:
                selected, reason = is_true_centromere_region(feature, attributes)
                if selected:
                    start_index = start - 1
                    end_index = end
                    bp_removed = max(0, end_index - start_index)
                    selected_intervals.setdefault(seqid, []).append((start_index, end_index))

            annotation_records.append(
                AnnotationRecord(
                    seqid=seqid,
                    source=source,
                    feature=feature,
                    start=start,
                    end=end,
                    attributes=attributes,
                    selected_for_removal=selected,
                    reason=reason,
                    bp_removed=bp_removed,
                )
            )

    merged_intervals = merge_and_clip_intervals(selected_intervals, original_records)
    total_bp_removed = sum(
        end - start
        for intervals in merged_intervals.values()
        for start, end in intervals
    )

    return AnnotationResult(
        records=annotation_records,
        intervals_by_seqid=merged_intervals,
        total_records=len(annotation_records),
        true_centromere_regions=sum(1 for record in annotation_records if record.selected_for_removal),
        total_centromere_bp_removed=total_bp_removed,
    )


def is_true_centromere_region(feature: str, attributes: str) -> tuple[bool, str]:
    """Return whether an annotation row marks the region itself as centromere."""
    feature_lower = feature.strip().lower()
    if feature_lower == "centromere":
        return True, 'feature type is exactly "centromere"'

    # Gene/protein annotations are not removed just because their product or
    # gene name mentions centromere, centrosome, or kinetochore.
    if feature_lower in GENE_LIKE_FEATURES:
        return False, "gene/protein feature type is not a centromere region"

    attributes_by_key = parse_attributes(attributes)
    for key, values in attributes_by_key.items():
        key_lower = key.lower()
        if key_lower in {"product", "gene", "gene_id", "transcript_id", "protein_id", "locus_tag"}:
            continue

        for value in values:
            value_lower = value.lower()
            if value_lower in {"centromere", "centromere region", "centromeric region"}:
                return True, f'attribute {key} marks a centromere region'
            if re.search(r"\bcentromeric\s+region\b", value_lower):
                return True, f'attribute {key} marks a centromeric region'
            if re.search(r"\bcentromere\s+region\b", value_lower):
                return True, f'attribute {key} marks a centromere region'
            if key_lower in {"gbkey", "rpt_type", "region_name"} and re.search(r"\bcentromere\b", value_lower):
                return True, f'attribute {key} marks a centromere region'

    return False, "not a true centromere-region annotation"


def parse_attributes(attributes: str) -> dict[str, list[str]]:
    """Parse simple GTF key \"value\" and GFF key=value attributes."""
    parsed: dict[str, list[str]] = {}
    for item in attributes.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
        else:
            parts = item.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts

        key = key.strip()
        value = value.strip().strip('"')
        parsed.setdefault(key, []).append(value)

    return parsed


def merge_and_clip_intervals(
    intervals_by_seqid: dict[str, list[tuple[int, int]]],
    original_records: dict[str, FastaRecord],
) -> dict[str, list[tuple[int, int]]]:
    """Merge selected intervals and clip them to known original contig lengths."""
    merged_by_seqid: dict[str, list[tuple[int, int]]] = {}

    for seqid, intervals in intervals_by_seqid.items():
        original_record = original_records.get(seqid)
        if original_record is None:
            merged_by_seqid[seqid] = []
            continue

        clipped_intervals = []
        sequence_length = len(original_record.sequence)
        for start, end in intervals:
            clipped_start = max(0, min(start, sequence_length))
            clipped_end = max(0, min(end, sequence_length))
            if clipped_end > clipped_start:
                clipped_intervals.append((clipped_start, clipped_end))

        if not clipped_intervals:
            merged_by_seqid[seqid] = []
            continue

        clipped_intervals.sort()
        merged = [clipped_intervals[0]]
        for start, end in clipped_intervals[1:]:
            last_start, last_end = merged[-1]
            if start <= last_end:
                merged[-1] = (last_start, max(last_end, end))
            else:
                merged.append((start, end))

        merged_by_seqid[seqid] = merged

    return merged_by_seqid


def subtract_intervals(
    sequence_length: int,
    removal_intervals: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Return half-open intervals left after removing centromere intervals."""
    if not removal_intervals:
        return [(0, sequence_length)]

    kept_intervals = []
    current_start = 0
    for remove_start, remove_end in removal_intervals:
        if current_start < remove_start:
            kept_intervals.append((current_start, remove_start))
        current_start = max(current_start, remove_end)

    if current_start < sequence_length:
        kept_intervals.append((current_start, sequence_length))

    return kept_intervals


def parse_split_id(record_id: str) -> SplitId | None:
    """Parse IDs like AOUM02000103.1:2470-24374."""
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
        contains_n = "N" in record.sequence.upper()
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
            "expected_length": split_id.end - split_id.start,
            "sequence_matches_original": False,
            "coordinate_convention": "",
        }

    zero_based = coordinate_candidate(
        original_record.sequence,
        split_id.start,
        split_id.end,
        chunk_sequence,
        "0-based half-open",
    )
    one_based = coordinate_candidate(
        original_record.sequence,
        split_id.start,
        split_id.end,
        chunk_sequence,
        "1-based inclusive",
    )

    for candidate in [zero_based, one_based]:
        if candidate["sequence_matches_original"]:
            return candidate

    if zero_based["coordinates_valid"]:
        return zero_based
    if one_based["coordinates_valid"]:
        return one_based
    return zero_based


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
    annotation_result: AnnotationResult | None,
) -> dict[str, FastaRecord]:
    """Split original records into continuous A/C/G/T-only chunks."""
    recreated_records: dict[str, FastaRecord] = {}
    removal_intervals_by_seqid = {}
    if annotation_result is not None:
        removal_intervals_by_seqid = annotation_result.intervals_by_seqid

    # Dict insertion order preserves the original FASTA record order in Python 3.7+.
    for original_id, original_record in original_records.items():
        sequence = original_record.sequence
        removal_intervals = removal_intervals_by_seqid.get(original_id, [])
        kept_intervals = subtract_intervals(len(sequence), removal_intervals)

        for kept_start, kept_end in kept_intervals:
            add_acgt_chunks_from_interval(
                recreated_records,
                original_id,
                sequence,
                kept_start,
                kept_end,
                min_chunk_length,
            )

    return recreated_records


def add_acgt_chunks_from_interval(
    recreated_records: dict[str, FastaRecord],
    original_id: str,
    sequence: str,
    interval_start: int,
    interval_end: int,
    min_chunk_length: int,
) -> None:
    """Split one kept interval at every non-ACGT base."""
    chunk_start_index = None

    for index in range(interval_start, interval_end):
        if sequence[index].upper() in ACGT:
            if chunk_start_index is None:
                chunk_start_index = index
        elif chunk_start_index is not None:
            add_recreated_chunk(
                recreated_records,
                original_id,
                sequence,
                chunk_start_index,
                index,
                min_chunk_length,
            )
            chunk_start_index = None

    if chunk_start_index is not None:
        add_recreated_chunk(
            recreated_records,
            original_id,
            sequence,
            chunk_start_index,
            interval_end,
            min_chunk_length,
        )


def add_recreated_chunk(
    recreated_records: dict[str, FastaRecord],
    original_id: str,
    sequence: str,
    start_index: int,
    end_index: int,
    min_chunk_length: int,
) -> None:
    """Add one half-open recreated chunk if it is long enough."""
    chunk_length = end_index - start_index
    if chunk_length < min_chunk_length:
        return

    record_id = f"{original_id}:{start_index}-{end_index}"
    chunk_sequence = sequence[start_index:end_index]

    recreated_records[record_id] = FastaRecord(
        record_id=record_id,
        header=record_id,
        sequence=chunk_sequence,
    )


def write_fasta(records: dict[str, FastaRecord], path: Path) -> None:
    """Write records to a FASTA file with one sequence line per record."""
    with path.open("w", encoding="utf-8", newline="\n") as output_file:
        for record in records.values():
            output_file.write(f">{record.header}\n")
            output_file.write(record.sequence + "\n")


def write_comparison_report(
    john_records: dict[str, FastaRecord] | None,
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
            "recreated_length",
            "opengenome_length",
            "sequences_match",
            "note",
        ]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        if john_records is None:
            writer.writerow(
                {
                    "comparison": "john_vs_opengenome",
                    "record_id": "N/A",
                    "status": "not_applicable",
                    "recreated_length": "",
                    "opengenome_length": "",
                    "sequences_match": "",
                    "note": "John-filtered FASTA not provided",
                }
            )
        else:
            write_pairwise_rows(
                writer,
                "john_vs_opengenome",
                john_records,
                opengenome_records,
                "John",
            )

        write_pairwise_rows(
            writer,
            "recreated_vs_opengenome",
            recreated_records,
            opengenome_records,
            "recreated",
        )


def write_pairwise_rows(
    writer: csv.DictWriter,
    comparison_name: str,
    left_records: dict[str, FastaRecord],
    opengenome_records: dict[str, FastaRecord],
    left_label: str,
) -> None:
    """Write rows for one pairwise FASTA comparison."""
    all_ids = sorted(set(left_records) | set(opengenome_records))

    for record_id in all_ids:
        left_record = left_records.get(record_id)
        opengenome_record = opengenome_records.get(record_id)

        if left_record is None:
            status = "only_in_opengenome"
            sequences_match = ""
            note = f"missing from {left_label}"
        elif opengenome_record is None:
            status = f"only_in_{left_label.lower()}"
            sequences_match = ""
            note = "missing from OpenGenome/Evo2"
        else:
            sequences_match = left_record.sequence == opengenome_record.sequence
            status = "shared_id_matching_sequence" if sequences_match else "shared_id_sequence_mismatch"
            note = ""

        writer.writerow(
            {
                "comparison": comparison_name,
                "record_id": record_id,
                "status": status,
                "recreated_length": len(left_record.sequence) if left_record else "",
                "opengenome_length": len(opengenome_record.sequence) if opengenome_record else "",
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


def write_annotation_filter_report(annotation_records: list[AnnotationRecord], path: Path) -> None:
    """Write the annotation selection decision for every parsed row."""
    fieldnames = [
        "seqid",
        "source",
        "feature",
        "start",
        "end",
        "attributes",
        "selected_for_removal",
        "reason",
        "bp_removed",
    ]

    with path.open("w", encoding="utf-8", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for record in annotation_records:
            writer.writerow(
                {
                    "seqid": record.seqid,
                    "source": record.source,
                    "feature": record.feature,
                    "start": record.start,
                    "end": record.end,
                    "attributes": record.attributes,
                    "selected_for_removal": record.selected_for_removal,
                    "reason": record.reason,
                    "bp_removed": record.bp_removed,
                }
            )


def write_summary(
    original_path: Path,
    john_path: Path | None,
    opengenome_path: Path,
    annotation_path: Path | None,
    annotation_result: AnnotationResult | None,
    original_stats: FastaStats,
    john_stats: FastaStats | None,
    opengenome_stats: FastaStats,
    recreated_stats: FastaStats,
    john_vs_opengenome: dict[str, object] | None,
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
    recreated_matches_actual = records_match_exactly(recreated_vs_opengenome)

    with path.open("w", encoding="utf-8", newline="\n") as summary_file:
        summary_file.write("Evo2/OpenGenome Filtering Reproducibility Summary\n")
        summary_file.write("================================================\n\n")

        summary_file.write("Input files\n")
        summary_file.write(f"- Original NCBI FASTA: {original_path}\n")
        if john_path is None:
            summary_file.write("- John's filtered FASTA: N/A (not provided)\n")
        else:
            summary_file.write(f"- John's filtered FASTA: {john_path}\n")
        summary_file.write(f"- OpenGenome/Evo2 filtered FASTA: {opengenome_path}\n")
        if annotation_path is None:
            summary_file.write("- Annotation GTF/GFF: N/A (not provided)\n")
            summary_file.write("- Annotation filtering provided: NO\n\n")
        else:
            summary_file.write(f"- Annotation GTF/GFF: {annotation_path}\n")
            summary_file.write("- Annotation filtering provided: YES\n")
            summary_file.write(
                f"- Total GTF/GFF records parsed: {annotation_result.total_records}\n"
            )
            summary_file.write(
                f"- True centromere regions found: {annotation_result.true_centromere_regions}\n"
            )
            summary_file.write(
                f"- Total centromere bp removed: {annotation_result.total_centromere_bp_removed}\n\n"
            )

        summary_file.write("FASTA summary statistics\n")
        write_stats_line(summary_file, "Original NCBI", original_stats)
        if john_stats is None:
            summary_file.write("- John filtered: N/A (not provided)\n")
        else:
            write_stats_line(summary_file, "John filtered", john_stats)
        write_stats_line(summary_file, "Official OpenGenome2", opengenome_stats)
        write_stats_line(summary_file, "Recreated OpenGenome-like", recreated_stats)
        summary_file.write("\n")

        summary_file.write("Strict sequence checks\n")
        summary_file.write(boolean_line("Official OpenGenome2 contains N", opengenome_stats.total_n_count > 0))
        summary_file.write(
            boolean_line(
                "Official OpenGenome2 contains non-ACGT characters",
                opengenome_stats.total_non_acgt_count > 0,
            )
        )
        summary_file.write(boolean_line("Recreated output contains N", recreated_stats.total_n_count > 0))
        summary_file.write(
            boolean_line(
                "Recreated output contains non-ACGT characters",
                recreated_stats.total_non_acgt_count > 0,
            )
        )
        summary_file.write("\n")

        summary_file.write("John vs OpenGenome/Evo2\n")
        if john_vs_opengenome is None:
            summary_file.write("- Comparison status: N/A (John-filtered FASTA not provided)\n\n")
        else:
            summary_file.write(f"- Records only in John: {len(john_vs_opengenome['only_left'])}\n")
            summary_file.write(f"- Records only in OpenGenome/Evo2: {len(john_vs_opengenome['only_right'])}\n")
            summary_file.write(f"- Records present in both: {len(john_vs_opengenome['shared_ids'])}\n")
            summary_file.write(
                "- Total bp difference (John minus OpenGenome/Evo2): "
                f"{john_vs_opengenome['left_total_bp'] - john_vs_opengenome['right_total_bp']}\n\n"
            )

        summary_file.write("OpenGenome/Evo2 chunk validation\n")
        summary_file.write(f"- OpenGenome-style split records detected: {split_records}\n")
        summary_file.write(
            "- Split records matching the original NCBI subsequence exactly: "
            f"{split_records_matching_original} of {split_records}\n"
        )
        summary_file.write(f"- Split records containing non-ACGT characters: {split_records_with_non_acgt}\n")
        summary_file.write(f"- Coordinate convention observed: {coordinate_conventions}\n\n")

        summary_file.write("Recreated OpenGenome-like comparison\n")
        summary_file.write(f"- Exact matching record IDs: {len(recreated_vs_opengenome['shared_ids'])}\n")
        summary_file.write(
            "- Exact matching sequences: "
            f"{len(recreated_vs_opengenome['matching_sequences'])}\n"
        )
        summary_file.write(f"- Missing records from recreated output: {recreated_missing}\n")
        summary_file.write(f"- Extra records in recreated output: {recreated_extra}\n")
        summary_file.write(f"- Sequence mismatches: {recreated_mismatches}\n")
        summary_file.write(
            boolean_line(
                "Recreated output matched official OpenGenome2 exactly",
                recreated_matches_actual,
            )
        )

        if recreated_matches_actual:
            summary_file.write(
                "- SUCCESS: Recreated OpenGenome/Evo2 output exactly at the record and sequence level.\n"
            )
        else:
            summary_file.write(
                "- NOT EXACT: Recreated output differs from OpenGenome/Evo2 output.\n"
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


def records_match_exactly(comparison: dict[str, object]) -> bool:
    """Return true when there are no missing, extra, or mismatched records."""
    return (
        len(comparison["only_left"]) == 0
        and len(comparison["only_right"]) == 0
        and len(comparison["sequence_mismatches"]) == 0
    )


def require_existing_file(path: Path, label: str) -> None:
    """Give a clear error if an expected input file is missing."""
    if not path.exists():
        raise FileNotFoundError(f"{label} file was not found: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"{label} path is not a file: {path}")


def build_arg_parser(
    default_original_path: Path,
    default_john_path: Path | None,
    default_opengenome_path: Path,
    default_annotation_path: Path | None,
    default_outdir: Path,
    description: str,
) -> argparse.ArgumentParser:
    """Create the command-line parser."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--original",
        type=Path,
        default=default_original_path,
        help=f"Original NCBI input FASTA. Default: {default_original_path}",
    )
    parser.add_argument(
        "--john",
        type=optional_path,
        default=None,
        help=(
            "Optional John's filtered output FASTA. "
            f"Default for this reproduction set: {default_john_path or 'N/A'}"
        ),
    )
    parser.add_argument(
        "--opengenome",
        type=Path,
        default=default_opengenome_path,
        help=f"OpenGenome/Evo2 filtered output FASTA. Default: {default_opengenome_path}",
    )
    parser.add_argument(
        "--annotation",
        type=optional_path,
        default=default_annotation_path,
        help=(
            "Optional GTF/GFF annotation file for conservative centromere removal. "
            f"Default: {default_annotation_path or 'N/A'}"
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=default_outdir,
        help=f"Directory where analysis outputs will be written. Default: {default_outdir}",
    )
    return parser



def run_workflow(
    default_original_path: Path,
    default_john_path: Path | None,
    default_opengenome_path: Path,
    default_outdir: Path,
    description: str,
    default_annotation_path: Path | None = None,
) -> None:
    """Run the FASTA comparison workflow."""
    args = build_arg_parser(
        default_original_path,
        default_john_path,
        default_opengenome_path,
        default_annotation_path,
        default_outdir,
        description,
    ).parse_args()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    if (
        args.john is None
        and default_john_path is not None
        and args.original == default_original_path
        and args.opengenome == default_opengenome_path
    ):
        args.john = default_john_path

    require_existing_file(args.original, "Original NCBI FASTA")
    if args.john is not None:
        require_existing_file(args.john, "John filtered FASTA")
    require_existing_file(args.opengenome, "OpenGenome/Evo2 filtered FASTA")
    annotation_path = args.annotation
    if annotation_path is not None and not annotation_path.is_file():
        print(f"Annotation file not found; skipping annotation filtering: {annotation_path}")
        annotation_path = None

    original_records = parse_fasta(args.original)
    john_records = parse_fasta(args.john) if args.john is not None else None
    opengenome_records = parse_fasta(args.opengenome)

    annotation_result = None
    if annotation_path is not None:
        annotation_result = parse_annotation_file(annotation_path, original_records)

    original_stats = calculate_stats(original_records)
    john_stats = calculate_stats(john_records) if john_records is not None else None
    opengenome_stats = calculate_stats(opengenome_records)

    john_vs_opengenome = (
        compare_record_sets(john_records, opengenome_records)
        if john_records is not None
        else None
    )
    validation_rows = validate_opengenome_chunks(original_records, opengenome_records)

    recreated_records = recreate_opengenome_like_records(
        original_records,
        min_chunk_length=MIN_CHUNK_LENGTH,
        annotation_result=annotation_result,
    )
    recreated_stats = calculate_stats(recreated_records)
    recreated_vs_opengenome = compare_record_sets(recreated_records, opengenome_records)

    recreated_fasta_path = outdir / "recreated_opengenome_like.fasta"
    summary_path = outdir / "summary.txt"
    comparison_report_path = outdir / "comparison_report.csv"
    chunk_validation_path = outdir / "opengenome_chunk_validation.csv"
    annotation_report_path = outdir / "annotation_filter_report.csv"

    write_fasta(recreated_records, recreated_fasta_path)
    write_comparison_report(
        john_records,
        opengenome_records,
        recreated_records,
        comparison_report_path,
    )
    write_chunk_validation_csv(validation_rows, chunk_validation_path)
    if annotation_result is not None:
        write_annotation_filter_report(annotation_result.records, annotation_report_path)

    write_summary(
        args.original,
        args.john,
        args.opengenome,
        annotation_path,
        annotation_result,
        original_stats,
        john_stats,
        opengenome_stats,
        recreated_stats,
        john_vs_opengenome,
        recreated_vs_opengenome,
        validation_rows,
        summary_path,
    )

    if opengenome_stats.total_non_acgt_count > 0:
        print(
            "WARNING: Official OpenGenome2 file contains non-ACGT characters. "
            "This may indicate the provided OpenGenome2 file is inconsistent "
            "with the expected strict filtering rule."
        )

    print(f"Wrote summary: {summary_path}")
    print(f"Wrote comparison report: {comparison_report_path}")
    print(f"Wrote chunk validation: {chunk_validation_path}")
    print(f"Wrote recreated FASTA: {recreated_fasta_path}")
    if annotation_result is not None:
        print(f"Wrote annotation filter report: {annotation_report_path}")

    if records_match_exactly(recreated_vs_opengenome):
        print("SUCCESS: Recreated OpenGenome/Evo2 output exactly at the record and sequence level.")
    else:
        print("NOT EXACT: Recreated output differs from OpenGenome/Evo2 output.")

    if recreated_stats.total_n_count == 0 and recreated_stats.total_non_acgt_count == 0:
        print("STRICT CLEAN OUTPUT: recreated FASTA contains only A/C/G/T bases.")
    else:
        print("STRICT CLEAN OUTPUT FAILED: recreated FASTA contains N or non-ACGT bases.")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATASET_DIR = PROJECT_ROOT / "Datasets" / "Evo2 Data Reproduction 1"
DEFAULT_ORIGINAL_PATH = DEFAULT_DATASET_DIR / "orginial_from_ncbi_GCA_000359685.2.fasta"
DEFAULT_JOHN_PATH = DEFAULT_DATASET_DIR / "john_filtered_GCA_000359685.2.fasta"
DEFAULT_OPENGENOME_PATH = DEFAULT_DATASET_DIR / "open_genome2_filtered_GCA_000359685.2.fasta"
DEFAULT_ANNOTATION_PATH = None
DEFAULT_OUTDIR = PROJECT_ROOT / "Results" / "Evo2 Data Reproduction"


def main() -> None:
    """Run this Evo2/OpenGenome2 comparison workflow."""
    run_workflow(
        default_original_path=DEFAULT_ORIGINAL_PATH,
        default_john_path=DEFAULT_JOHN_PATH,
        default_opengenome_path=DEFAULT_OPENGENOME_PATH,
        default_annotation_path=DEFAULT_ANNOTATION_PATH,
        default_outdir=DEFAULT_OUTDIR,
        description="Compare John and OpenGenome/Evo2 filtered FASTA outputs.",
    )


if __name__ == "__main__":
    main()


