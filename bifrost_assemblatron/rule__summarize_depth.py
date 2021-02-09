# script for use with snakemake
import traceback
from typing import Dict
from bifrostlib import common

def rule__summarize_depth(input: object, output: object, component_json: Dict, log: object) -> None:
    try:
        # Variables being used
        """
        Break this into 2?
        """
        depth_dict = {}
        with open(input.coverage, "r") as input_file:
            for line in input_file:
                if line[0] != "#":
                    contig = line.split("\t")[0]
                    depth = int(line.split("\t")[2].strip())
                    if contig not in depth_dict:
                        depth_dict[contig] = {}
                    if depth in depth_dict[contig]:
                        depth_dict[contig][depth] += 1
                    else:
                        depth_dict[contig][depth] = 1

        contig_depth_summary_dict = {}
        for contig in depth_dict:
            total_depth = 0
            total_length = 0
            for depth in depth_dict[contig]:
                length = depth_dict[contig][depth]
                total_depth = total_depth + (length * depth)
                total_length = total_length + length
            contig_depth_summary_dict[contig] = {
                "coverage": total_depth / total_length,
                "total_depth": total_depth,
                "total_length": total_length
            }

        # Removing as it's also looped in vcf so no need to check twice
        # dict is made now to cycle over on range and on contigs
        # binned_depth_summary_dict = {}
        binned_depth = [0] * 100
        # depth_limits = config["serum"]["summarize"]["depth_range"]
        # depth_range = list(range(depth_limits[0], depth_limits[1]))
        # for bound in depth_range:
        #     binned_depth_summary_dict[bound] = 0

        for contig in depth_dict:
            for depth in depth_dict[contig]:
                for i in range(1, 100):
                    if depth >= i:
                        binned_depth[i - 1] += depth_dict[contig][depth]

        common.save_yaml({"contig_depth": contig_depth_summary_dict}, output._file)
        common.save_yaml({"binned_depth": binned_depth}, output._file2)
    except Exception:
        with open(log.err_file, "w+") as fh:
            fh.write(traceback.format_exc())


rule__greater_than_min_reads_check(
    snakemake.input,
    snakemake.output,
    snakemake.params.component_json,
    snakemake.log)
