#!/usr/bin/env python3
import argparse
import re
import sys
from typing import Dict, List, Tuple

import numpy as np

PAT = re.compile(r"length_(\d+)_cov_([0-9.]+)")


def n50(lengths: List[int]) -> int:
    if not lengths:
        return 0
    total = sum(lengths)
    half = total / 2.0
    cummulative = 0
    for Length in sorted(lengths, reverse=True):
        cummulative += Length
        if cummulative >= half:
            return Length
    return 0


def contig_header_filter(
    header_line: str,
    contig_len: int,
    gc: int,
    acgt: int,
    groups: Dict[str, Dict[str, object]],
    cov_threshold: float,
    min_contig_len: int,
    pat: re.Pattern[str] = PAT,
) -> Tuple[str, int, int]:
    """
    Returns:
      (status, bad_header_delta, filtered_short_delta)

    status:
      - "passed"          : header OK, length >= min_contig_len, coverage >= cov_threshold
      - "failed_coverage" : header OK, length >= min_contig_len, coverage <  cov_threshold
      - "failed_length"   : bad header OR too short

    Side-effect:
      Updates `groups` for contigs with header OK and length >= min_contig_len (both cov bins).
      - groups["cov_lower"]   == length-pass + cov <  cov_threshold  (i.e. failed_coverage)
      - groups["cov_greater"] == length-pass + cov >= cov_threshold  (i.e. passed)
    """
    m = pat.search(header_line)
    if not m:
        return ("failed_length", 1, 0)

    if contig_len < min_contig_len:
        return ("failed_length", 0, 1)

    cov = float(m.group(2))

    key = "cov_lower" if cov < cov_threshold else "cov_greater"
    g = groups[key]
    g["n"] = int(g["n"]) + 1
    g["total_len"] = int(g["total_len"]) + contig_len
    g["total_gc"] = int(g["total_gc"]) + gc
    g["total_acgt"] = int(g["total_acgt"]) + acgt
    g["lengths"].append(contig_len)  # type: ignore[attr-defined]

    return ("failed_coverage" if cov < cov_threshold else "passed", 0, 0)


def filter_assembly(
    assembly_path: str,
    groups: Dict[str, Dict[str, object]],
    cov_threshold: float,
    min_contig_len: int,
    passed_fasta: str,
    failed_length_fasta: str | None = None,
    failed_coverage_fasta: str | None = None,
) -> Tuple[int, int, int, Dict[str, object]]:
    """
    Stream the assembly FASTA, classify contigs, update `groups`, and write FASTA outputs.

    Returns:
      (bad_headers, filtered_short, failed_coverage_count, failed_length_stats)

    Notes:
      - passed_fasta is always written
      - failed_length_fasta / failed_coverage_fasta are only written if provided
    """
    bad_headers = 0
    filtered_short = 0
    failed_coverage_count = 0

    failed_length_stats: Dict[str, object] = {"n": 0, "total_len": 0, "total_gc": 0, "total_acgt": 0, "lengths": []}

    passed_fh = open(passed_fasta, "w")
    failed_len_fh = open(failed_length_fasta, "w") if failed_length_fasta else None
    failed_cov_fh = open(failed_coverage_fasta, "w") if failed_coverage_fasta else None

    try:
        header: str | None = None
        seq_lines: List[str] = []
        contig_len = 0
        gc = 0
        acgt = 0

        def update_failed_length_stats() -> None:
            failed_length_stats["n"] = int(failed_length_stats["n"]) + 1
            failed_length_stats["total_len"] = int(failed_length_stats["total_len"]) + contig_len
            failed_length_stats["total_gc"] = int(failed_length_stats["total_gc"]) + gc
            failed_length_stats["total_acgt"] = int(failed_length_stats["total_acgt"]) + acgt
            failed_length_stats["lengths"].append(contig_len)  # type: ignore[attr-defined]

        with open(assembly_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    # finalize previous contig ONLY if we already have one
                    if header is not None:
                        status, d_bad, d_short = contig_header_filter(
                            header,
                            contig_len,
                            gc,
                            acgt,
                            groups,
                            cov_threshold=cov_threshold,
                            min_contig_len=min_contig_len,
                        )
                        bad_headers += d_bad
                        filtered_short += d_short
                        if status == "failed_coverage":
                            failed_coverage_count += 1
                        if status == "failed_length":
                            update_failed_length_stats()

                        # write FASTA record (only if that output was requested)
                        if status == "passed":
                            out_fh = passed_fh
                        elif status == "failed_coverage":
                            out_fh = failed_cov_fh
                        else:  # failed_length
                            out_fh = failed_len_fh

                        if out_fh is not None:
                            out_fh.write(header + "\n")
                            for sline in seq_lines:
                                out_fh.write(sline + "\n")

                    # start new contig
                    header = line.strip()
                    seq_lines = []
                    contig_len = 0
                    gc = 0
                    acgt = 0
                else:
                    raw = line.strip()
                    if not raw:
                        continue

                    seq_lines.append(raw)

                    s = raw.upper()
                    contig_len += len(s)  # includes N etc in contig length
                    a = s.count("A")
                    c = s.count("C")
                    g_ = s.count("G")
                    t = s.count("T")
                    gc += c + g_
                    acgt += a + c + g_ + t  # excludes N etc from denominator

            # finalize last contig (if any)
            if header is not None:
                status, d_bad, d_short = contig_header_filter(
                    header,
                    contig_len,
                    gc,
                    acgt,
                    groups,
                    cov_threshold=cov_threshold,
                    min_contig_len=min_contig_len,
                )
                bad_headers += d_bad
                filtered_short += d_short
                if status == "failed_coverage":
                    failed_coverage_count += 1
                if status == "failed_length":
                    update_failed_length_stats()

                if status == "passed":
                    out_fh = passed_fh
                elif status == "failed_coverage":
                    out_fh = failed_cov_fh
                else:
                    out_fh = failed_len_fh

                if out_fh is not None:
                    out_fh.write(header + "\n")
                    for sline in seq_lines:
                        out_fh.write(sline + "\n")

    finally:
        passed_fh.close()
        if failed_len_fh:
            failed_len_fh.close()
        if failed_cov_fh:
            failed_cov_fh.close()

    return bad_headers, filtered_short, failed_coverage_count, failed_length_stats


def calculate_stats(
    output_path: str,
    label: str,
    stats: Dict[str, object],
    cov_threshold: float,
    min_contig_len: int,
    stdout: bool = False,
) -> None:
    n = int(stats["n"])
    total_len = int(stats["total_len"])
    total_gc = int(stats["total_gc"])
    total_acgt = int(stats["total_acgt"])
    lengths = stats["lengths"]  # type: ignore[assignment]

    gc_frac = np.nan if total_acgt == 0 else (total_gc / total_acgt)
    gc_frac_str = "NA" if np.isnan(gc_frac) else f"{gc_frac:.6f}"
    gc_pct_str = "NA" if np.isnan(gc_frac) else f"{(gc_frac * 100.0):.3f}"

    with open(output_path, "w") as out:
        out.write("group\tcov_threshold\tmin_contig_len\tcontigs\tsum_len\tgc_fraction\tgc_percent\tn50\n")
        out.write(
            f"{label}\t{cov_threshold:g}\t{min_contig_len}\t{n}\t{total_len}\t{gc_frac_str}\t{gc_pct_str}\t{n50(lengths)}\n"
        )

    if stdout:
        print(f"{label}\tcontigs={n}\tsum_len={total_len}\tGC={gc_pct_str}%\tN50={n50(lengths)}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Filter contigs by min length and coverage threshold; write FASTA outputs and per-output stats TSVs."
    )
    parser.add_argument("--assembly", required=True, help="Input assembly FASTA.")
    parser.add_argument(
        "--cov-threshold",
        type=float,
        default=10.0,
        help="Coverage threshold for passing (>=thr). Default: 10x",
    )
    parser.add_argument(
        "--min-contig-len",
        type=int,
        default=500,
        help="Minimum contig length to keep. Default: 500",
    )

    # Prefix outputs (code appends .fasta and _stat.tsv)
    parser.add_argument(
        "--passed",
        required=True,
        nargs="?",
        const="assembly_pass",
        help="PREFIX for passing contigs. Writes {prefix}.fasta and {prefix}_stat.tsv. "
             "If provided without a value, defaults to: assembly_pass",
    )
    parser.add_argument(
        "--failed-length",
        nargs="?",
        const="assembly_length_fail",
        default=None,
        help="Optional PREFIX for contigs that fail due to bad header and/or too short. "
             "Writes {prefix}.fasta and {prefix}_stat.tsv. "
             "If provided without a value, defaults to: assembly_length_fail",
    )
    parser.add_argument(
        "--failed-coverage",
        nargs="?",
        const="assembly_cov_fail",
        default=None,
        help="Optional PREFIX for contigs that pass length but fail coverage (<thr). "
             "Writes {prefix}.fasta and {prefix}_stat.tsv. "
             "If provided without a value, defaults to: assembly_cov_fail",
    )

    parser.add_argument("--stdout", action="store_true", help="Also print summaries to stdout.")
    args = parser.parse_args()

    cov_threshold = float(args.cov_threshold)
    min_contig_len = int(args.min_contig_len)

    # Build output paths from prefixes
    passed_prefix = args.passed
    passed_fasta = f"{passed_prefix}.fasta"
    passed_stat = f"{passed_prefix}_stat.tsv"

    failed_len_prefix = args.failed_length
    failed_len_fasta = f"{failed_len_prefix}.fasta" if failed_len_prefix else None
    failed_len_stat = f"{failed_len_prefix}_stat.tsv" if failed_len_prefix else None

    failed_cov_prefix = args.failed_coverage
    failed_cov_fasta = f"{failed_cov_prefix}.fasta" if failed_cov_prefix else None
    failed_cov_stat = f"{failed_cov_prefix}_stat.tsv" if failed_cov_prefix else None

    # Coverage-bin stats for length-pass contigs (used to make passed/failed_coverage stats)
    groups: Dict[str, Dict[str, object]] = {
        "cov_lower": {"n": 0, "total_len": 0, "total_gc": 0, "total_acgt": 0, "lengths": []},
        "cov_greater": {"n": 0, "total_len": 0, "total_gc": 0, "total_acgt": 0, "lengths": []},
    }

    bad_headers, filtered_short, failed_coverage_count, failed_length_stats = filter_assembly(
        assembly_path=args.assembly,
        groups=groups,
        cov_threshold=cov_threshold,
        min_contig_len=min_contig_len,
        passed_fasta=passed_fasta,
        failed_length_fasta=failed_len_fasta,
        failed_coverage_fasta=failed_cov_fasta,
    )

    # Stats per created output (same prefix -> {prefix}_stat.tsv)
    calculate_stats(
        output_path=passed_stat,
        label="passed",
        stats=groups["cov_greater"],
        cov_threshold=cov_threshold,
        min_contig_len=min_contig_len,
        stdout=args.stdout,
    )

    if failed_len_stat is not None:
        calculate_stats(
            output_path=failed_len_stat,
            label="failed_length",
            stats=failed_length_stats,
            cov_threshold=cov_threshold,
            min_contig_len=min_contig_len,
            stdout=args.stdout,
        )

    if failed_cov_stat is not None:
        calculate_stats(
            output_path=failed_cov_stat,
            label="failed_coverage",
            stats=groups["cov_lower"],
            cov_threshold=cov_threshold,
            min_contig_len=min_contig_len,
            stdout=args.stdout,
        )

    # Messages to stderr
    if filtered_short:
        print(f"INFO: filtered out {filtered_short} contig(s) with length < {min_contig_len}.", file=sys.stderr)
    if bad_headers:
        print(f"WARNING: {bad_headers} header(s) did not match 'length_..._cov_...' pattern.", file=sys.stderr)
    if failed_cov_fasta and failed_coverage_count:
        print(
            f"INFO: wrote {failed_coverage_count} contig(s) failing coverage (<{cov_threshold:g}) to {failed_cov_fasta}.",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
