#!/usr/bin/env python3
"""
Assembly QC with overlapping buckets:

A: len_lt_500        (length < 500)
B: cov_lt_1          (length >= 500 AND cov < 1)
C: cov_ge_1          (length >= 500 AND cov >= 1)
D: cov_ge_10         (length >= 500 AND cov >= 10)

A contig may belong to multiple buckets (except B vs C).
"""

import argparse
import re
from typing import Dict, List, Tuple, Optional, TextIO

PAT = re.compile(r"length_(\d+)_cov_([0-9.]+)")


def n50(lengths: List[int]) -> int:
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
    return {"no_contigs": 0, "total_len": 0, "cov_x_len_sum": 0.0, "lengths": []}


def update_stats(stats: Dict[str, object], contig_len: int, coverage: Optional[float]) -> None:
    stats["no_contigs"] += 1
    stats["total_len"] += contig_len
    stats["lengths"].append(contig_len)
    if coverage is not None:
        stats["cov_x_len_sum"] += coverage * contig_len


def fasta_records(path: str):
    header = None
    seq_lines = []
    contig_len = 0

    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield (header, seq_lines, contig_len)
                header = line.strip()
                seq_lines = []
                contig_len = 0
            else:
                raw = line.strip()
                if raw:
                    seq_lines.append(raw)
                    contig_len += len(raw)

        if header is not None:
            yield (header, seq_lines, contig_len)


def classify_groups(
    header_line: str,
    contig_len: int,
    cov_threshold: float,
    min_contig_len: int,
) -> Tuple[List[str], Optional[float]]:
    """
    Returns a list of groups the contig belongs to.
    Groups are overlapping.
    """
    m = PAT.search(header_line)
    if not m:
        # If header cannot be parsed, treat as length-only
        if contig_len < min_contig_len:
            return ["len_lt_500"], None
        return [], None

    cov = float(m.group(2))
    groups = []

    if contig_len < min_contig_len:
        groups.append("len_lt_500")
        return groups, cov

    # length >= 500
    if cov < 1:
        groups.append("cov_lt_1")
    if cov >= 1:
        groups.append("cov_ge_1")
    if cov >= cov_threshold:
        groups.append("cov_ge_10")

    return groups, cov


def process_record(
    header: str,
    seq_lines: List[str],
    contig_len: int,
    stats_by_group: Dict[str, Dict[str, object]],
    out_by_group: Dict[str, Optional[TextIO]],
    cov_threshold: float,
    min_contig_len: int,
) -> None:

    groups, cov = classify_groups(
        header_line=header,
        contig_len=contig_len,
        cov_threshold=cov_threshold,
        min_contig_len=min_contig_len,
    )

    for g in groups:
        update_stats(stats_by_group[g], contig_len, cov)
        fh = out_by_group.get(g)
        if fh is not None:
            fh.write(header + "\n")
            for s in seq_lines:
                fh.write(s + "\n")


def filter_assembly(
    assembly_path: str,
    stats_by_group: Dict[str, Dict[str, object]],
    cov_threshold: float,
    min_contig_len: int,
    out_by_group: Dict[str, Optional[TextIO]],
) -> None:

    for header, seq_lines, contig_len in fasta_records(assembly_path):
        process_record(
            header=header,
            seq_lines=seq_lines,
            contig_len=contig_len,
            stats_by_group=stats_by_group,
            out_by_group=out_by_group,
            cov_threshold=cov_threshold,
            min_contig_len=min_contig_len,
        )


def calculate_stats(
    output_path: str,
    label: str,
    stats: Dict[str, object],
    cov_threshold: float,
    min_contig_len: int,
    stdout: bool = False,
) -> None:

    no_contigs = stats["no_contigs"]
    total_len = stats["total_len"]
    lengths = stats["lengths"]
    cov_x_len_sum = stats["cov_x_len_sum"]

    mean_cov = cov_x_len_sum / total_len if total_len > 0 else 0.0

    with open(output_path, "w") as out:
        out.write("group\tcov_threshold\tmin_contig_len\tno_contigs\tsum_len\tn50\tmean_cov\n")
        out.write(
            f"{label}\t{cov_threshold:g}\t{min_contig_len}\t{no_contigs}\t"
            f"{total_len}\t{n50(lengths)}\t{mean_cov:.2f}\n"
        )

    if stdout:
        print(
            f"{label}\tcontigs={no_contigs}\tsum_len={total_len}\t"
            f"N50={n50(lengths)}\tmean_cov={mean_cov:.2f}"
        )


def parser():
    p = argparse.ArgumentParser()
    p.add_argument("--assembly", required=True)
    p.add_argument("--cov-threshold", type=float, default=10.0)
    p.add_argument("--min-contig-len", type=int, default=500)

    p.add_argument("--passed", required=True)
    p.add_argument("--failed-length")
    p.add_argument("--failed-coverage")
    p.add_argument("--low-cov-above-1x")
    p.add_argument("--stdout", action="store_true")
    return p.parse_args()


def main():
    args = parser()

    cov_threshold = args.cov_threshold
    min_contig_len = args.min_contig_len

    # Map CLI prefixes to group names
    out_by_group = {
        "cov_ge_10": open(args.passed + ".fasta", "w"),
        "cov_ge_1": open(args.low_cov_above_1x + ".fasta", "w") if args.low_cov_above_1x else None,
        "cov_lt_1": open(args.failed_coverage + ".fasta", "w") if args.failed_coverage else None,
        "len_lt_500": open(args.failed_length + ".fasta", "w") if args.failed_length else None,
    }

    stats_by_group = {
        "cov_ge_10": init_stats(),
        "cov_ge_1": init_stats(),
        "cov_lt_1": init_stats(),
        "len_lt_500": init_stats(),
    }

    filter_assembly(
        assembly_path=args.assembly,
        stats_by_group=stats_by_group,
        cov_threshold=cov_threshold,
        min_contig_len=min_contig_len,
        out_by_group=out_by_group,
    )

    # Close FASTA files
    for fh in out_by_group.values():
        if fh:
            fh.close()

    # Write stats
    if args.passed:
        calculate_stats(args.passed + "_stat.tsv", "cov_ge_10",
                        stats_by_group["cov_ge_10"], cov_threshold, min_contig_len, args.stdout)

    if args.low_cov_above_1x:
        calculate_stats(args.low_cov_above_1x + "_stat.tsv", "cov_ge_1",
                        stats_by_group["cov_ge_1"], cov_threshold, min_contig_len, args.stdout)

    if args.failed_coverage:
        calculate_stats(args.failed_coverage + "_stat.tsv", "cov_lt_1",
                        stats_by_group["cov_lt_1"], cov_threshold, min_contig_len, args.stdout)

    if args.failed_length:
        calculate_stats(args.failed_length + "_stat.tsv", "len_lt_500",
                        stats_by_group["len_lt_500"], cov_threshold, min_contig_len, args.stdout)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

