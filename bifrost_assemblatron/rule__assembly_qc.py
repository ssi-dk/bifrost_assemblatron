#!/usr/bin/env python3
"""
Filter contigs from an assembly FASTA by:
  1) minimum contig length
  2) coverage threshold parsed from the contig header

Outputs:
  - passed contigs FASTA (+ stats TSV)
  - optional failed_length FASTA (+ stats TSV)
  - optional failed_coverage FASTA (+ stats TSV)

Assumes headers shares the pattern obtained by spades assembly
i.e contain a substring like: length_<N>_cov_<COV>
(e.g. >contig_1_length_1234_cov_17.3)
"""

import argparse
import re
from typing import Dict, List, Tuple, Optional, TextIO

# Regex to extract length and coverage from headers
PAT = re.compile(r"length_(\d+)_cov_([0-9.]+)")


def n50(lengths: List[int]) -> int:
    """
    Compute N50 from a list of contig lengths.

    N50 is the contig length L such that 50% of the total assembly length
    is contained in contigs of length >= L.
    """
    if not lengths:
        return 0
    total = sum(lengths)
    half = total / 2.0
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= half:
            return length
    return 0


def init_stats() -> Dict[str, object]:
    """
    Create a stats dict for a group of contigs.

      - no_contigs: number of contigs
      - total_len: sum of contig lengths
      - cov_x_len_sum: sum of (coverage * contig_len) across contigs where coverage was parsed
      - lengths: list of all contig lengths (used for N50)
    """
    return {"no_contigs": 0, "total_len": 0, "cov_x_len_sum": 0.0, "lengths": []}


def update_stats(stats: Dict[str, object], contig_len: int, coverage: Optional[float]) -> None:
    """
    Update stats for one contig.
    Shared for ALL buckets (passed, failed_length, failed_coverage).

    If coverage is None (header didn't match), cov_x_len_sum is not updated.
    """
    stats["no_contigs"] = int(stats["no_contigs"]) + 1
    stats["total_len"] = int(stats["total_len"]) + contig_len
    stats["lengths"].append(contig_len)  # type: ignore[attr-defined]

    if coverage is None:
        return

    stats["cov_x_len_sum"] = float(stats["cov_x_len_sum"]) + (coverage * contig_len)


def classify_contig(
    header_line: str,
    contig_len: int,
    cov_threshold: float,
    min_contig_len: int,
    pat: re.Pattern[str] = PAT,
) -> Tuple[str, int, int, Optional[float]]:
    """
    Classify a contig into one of three buckets:
      - "passed"          : length OK and coverage >= threshold
      - "failed_coverage" : length OK but coverage < threshold
      - "failed_length"   : header doesn't parse coverage OR contig too short

    Returns:
      (status, bad_header_delta, filtered_short_delta, coverage_or_none)

    Counters:
      - bad_header_delta: +1 if we couldn't parse coverage from the header
      - filtered_short_delta: +1 if header parsed but contig_len < min_contig_len
    """
    m = pat.search(header_line)
    if not m:
        # Can't parse coverage -> failed_length, coverage unknown
        return ("failed_length", 1, 0, None)

    cov = float(m.group(2))

    # Length filter comes before coverage filter (only long-enough contigs
    # are eligible for passed/failed_coverage)
    if contig_len < min_contig_len:
        return ("failed_length", 0, 1, cov)

    # Coverage threshold check
    if cov < cov_threshold:
        return ("failed_coverage", 0, 0, cov)

    return ("passed", 0, 0, cov)


def write_fasta_record(fh: TextIO, header: str, seq_lines: List[str]) -> None:
    """
    Write a single FASTA record to the given file handle.
    seq_lines are written exactly as read (no re-wrapping).
    """
    fh.write(header + "\n")
    for sline in seq_lines:
        fh.write(sline + "\n")


def fasta_records(path: str) -> Tuple[str, List[str], int]:
    """
    Generator yielding FASTA records as:
      (header, seq_lines, contig_len)

    contig_len is computed as the sum of the lengths of sequence lines.
    Blank lines are ignored.
    """
    header: Optional[str] = None
    seq_lines: List[str] = []
    contig_len = 0

    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Emit the previous record when we see a new header
                if header is not None:
                    yield (header, seq_lines, contig_len)

                # Start a new record
                header = line.strip()
                seq_lines = []
                contig_len = 0
            else:
                raw = line.strip()
                if not raw:
                    continue  # ignore blank lines

                seq_lines.append(raw)
                contig_len += len(raw)  # includes N etc. in length

        # Emit the last record at EOF
        if header is not None:
            yield (header, seq_lines, contig_len)


def process_record(
    header: str,
    seq_lines: List[str],
    contig_len: int,
    stats_by_status: Dict[str, Dict[str, object]],
    out_by_status: Dict[str, Optional[TextIO]],
    cov_threshold: float,
    min_contig_len: int,
) -> Tuple[int, int, int]:
    """
    Process one FASTA record:
      1) classify into a bucket
      2) update stats for that bucket (including cov_x_len_sum if coverage is known)
      3) write the record to the corresponding FASTA (if enabled)

    Returns deltas for global counters:
      (bad_header_delta, filtered_short_delta, failed_coverage_delta)
    """
    status, d_bad, d_short, cov = classify_contig(
        header,
        contig_len,
        cov_threshold=cov_threshold,
        min_contig_len=min_contig_len,
    )

    # One shared stats updater for ALL buckets
    update_stats(stats_by_status[status], contig_len, cov)

    # Write record only if that output bucket was requested (fh != None)
    out_fh = out_by_status.get(status)
    if out_fh is not None:
        write_fasta_record(out_fh, header, seq_lines)

    failed_cov_delta = 1 if status == "failed_coverage" else 0
    return d_bad, d_short, failed_cov_delta


def filter_assembly(
    assembly_path: str,
    stats_by_status: Dict[str, Dict[str, object]],
    cov_threshold: float,
    min_contig_len: int,
    passed_fasta: str,
    failed_length_fasta: str | None = None,
    failed_coverage_fasta: str | None = None,
) -> Tuple[int, int, int]:
    """
    Stream the FASTA file once and route each contig to:
      - passed
      - failed_length
      - failed_coverage

    Returns:
      (bad_headers, filtered_short, failed_coverage_count)
    """
    bad_headers = 0
    filtered_short = 0
    failed_coverage_count = 0

    # Output FASTAs (passed is always written; the others are optional)
    passed_fh = open(passed_fasta, "w")
    failed_len_fh = open(failed_length_fasta, "w") if failed_length_fasta else None
    failed_cov_fh = open(failed_coverage_fasta, "w") if failed_coverage_fasta else None

    # Map classification bucket -> output file handle (or None to skip writing)
    out_by_status: Dict[str, Optional[TextIO]] = {
        "passed": passed_fh,
        "failed_length": failed_len_fh,
        "failed_coverage": failed_cov_fh,
    }

    try:
        for header, seq_lines, contig_len in fasta_records(assembly_path):
            d_bad, d_short, d_cov = process_record(
                header=header,
                seq_lines=seq_lines,
                contig_len=contig_len,
                stats_by_status=stats_by_status,
                out_by_status=out_by_status,
                cov_threshold=cov_threshold,
                min_contig_len=min_contig_len,
            )
            # Accumulate global counters for QA/logging
            bad_headers += d_bad
            filtered_short += d_short
            failed_coverage_count += d_cov
    finally:
        # Always close outputs
        passed_fh.close()
        if failed_len_fh:
            failed_len_fh.close()
        if failed_cov_fh:
            failed_cov_fh.close()

    return bad_headers, filtered_short, failed_coverage_count


def calculate_stats(
    output_path: str,
    label: str,
    stats: Dict[str, object],
    cov_threshold: float,
    min_contig_len: int,
    stdout: bool = False,
) -> None:
    """
    Write a 1-line TSV stats file for a given bucket.

    Columns:
      group, cov_threshold, min_contig_len, no_contigs, sum_len, n50, cov_x_len_sum, mean_cov

    mean_cov is the length-weighted mean coverage:
      mean_cov = cov_x_len_sum / sum_len
    """
    no_contigs = int(stats["no_contigs"])
    total_len = int(stats["total_len"])
    lengths = stats["lengths"]  # type: ignore[assignment]
    cov_x_len_sum = float(stats["cov_x_len_sum"])

    mean_cov = (cov_x_len_sum / total_len) if total_len > 0 else 0.0

    with open(output_path, "w") as out:
        out.write(
            "group\tcov_threshold\tmin_contig_len\tno_contigs\tsum_len\tn50\tmean_cov\n"
        )
        out.write(
            f"{label}\t{cov_threshold:g}\t{min_contig_len}\t{no_contigs}\t{total_len}\t{n50(lengths)}\t{mean_cov:.2f}\n"
        )

    if stdout:
        print(
            f"{label}\tcontigs={no_contigs}\tsum_len={total_len}\tN50={n50(lengths)}\tmean_cov={mean_cov:.2f}"
        )

def parser():
    """
    CLI parser.
    Output arguments are PREFIXES; we write:
      {prefix}.fasta
      {prefix}_stat.tsv
    """
    p = argparse.ArgumentParser(
        description="Filter contigs by min length and coverage threshold; write FASTA outputs and per-output stats TSVs."
    )
    p.add_argument("--assembly", required=True, help="Input assembly FASTA.")
    p.add_argument(
        "--cov-threshold",
        type=float,
        default=10.0,
        help="Coverage threshold for passing (>=thr). Default: 10x",
    )
    p.add_argument(
        "--min-contig-len",
        type=int,
        default=500,
        help="Minimum contig length to keep. Default: 500",
    )

    # Prefix outputs (code appends .fasta and _stat.tsv)
    p.add_argument(
        "--passed",
        required=True,
        nargs="?",
        const="assembly_pass",
        help="PREFIX for passing contigs. Writes {prefix}.fasta and {prefix}_stat.tsv. "
        "If provided without a value, defaults to: assembly_pass",
    )
    p.add_argument(
        "--failed-length",
        nargs="?",
        const="assembly_length_fail",
        default=None,
        help="Optional PREFIX for contigs that fail due to bad header and/or too short. "
        "Writes {prefix}.fasta and {prefix}_stat.tsv. "
        "If provided without a value, defaults to: assembly_length_fail",
    )
    p.add_argument(
        "--failed-coverage",
        nargs="?",
        const="assembly_cov_fail",
        default=None,
        help="Optional PREFIX for contigs that pass length but fail coverage (<thr). "
        "Writes {prefix}.fasta and {prefix}_stat.tsv. "
        "If provided without a value, defaults to: assembly_cov_fail",
    )

    p.add_argument("--stdout", action="store_true", help="Also print summaries to stdout.")
    return p.parse_args()


def main() -> int:
    args = parser()
    cov_threshold = float(args.cov_threshold)
    min_contig_len = int(args.min_contig_len)

    # Build output paths from user-provided prefixes
    passed_prefix = args.passed
    passed_fasta = f"{passed_prefix}.fasta"
    passed_stat = f"{passed_prefix}_stat.tsv"

    failed_len_prefix = args.failed_length
    failed_len_fasta = f"{failed_len_prefix}.fasta" if failed_len_prefix else None
    failed_len_stat = f"{failed_len_prefix}_stat.tsv" if failed_len_prefix else None

    failed_cov_prefix = args.failed_coverage
    failed_cov_fasta = f"{failed_cov_prefix}.fasta" if failed_cov_prefix else None
    failed_cov_stat = f"{failed_cov_prefix}_stat.tsv" if failed_cov_prefix else None

    # Stats dicts for each bucket
    stats_by_status: Dict[str, Dict[str, object]] = {
        "passed": init_stats(),
        "failed_length": init_stats(),
        "failed_coverage": init_stats(),
    }

    # Stream input once and write FASTA outputs
    bad_headers, filtered_short, failed_coverage_count = filter_assembly(
        assembly_path=args.assembly,
        stats_by_status=stats_by_status,
        cov_threshold=cov_threshold,
        min_contig_len=min_contig_len,
        passed_fasta=passed_fasta,
        failed_length_fasta=failed_len_fasta,
        failed_coverage_fasta=failed_cov_fasta,
    )

    # Write stats TSV(s) corresponding to outputs created
    calculate_stats(
        output_path=passed_stat,
        label="passed",
        stats=stats_by_status["passed"],
        cov_threshold=cov_threshold,
        min_contig_len=min_contig_len,
        stdout=args.stdout,
    )

    if failed_len_stat is not None:
        calculate_stats(
            output_path=failed_len_stat,
            label="failed_length",
            stats=stats_by_status["failed_length"],
            cov_threshold=cov_threshold,
            min_contig_len=min_contig_len,
            stdout=args.stdout,
        )

    if failed_cov_stat is not None:
        calculate_stats(
            output_path=failed_cov_stat,
            label="failed_coverage",
            stats=stats_by_status["failed_coverage"],
            cov_threshold=cov_threshold,
            min_contig_len=min_contig_len,
            stdout=args.stdout,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
