from bifrostlib import common
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Category
from typing import Dict
import os


def extract_contigs_sum_cov(denovo_assembly: Category, mapping_qc: Category, results: Dict, component_name: str) -> None:
    file_name = "contigs.sum.cov"
    file_key = common.json_key_cleaner(file_name)
    file_path = os.path.join(component_name, file_name)

    yaml = common.get_yaml(file_path)
    contig_summary_yaml = yaml["contig_depth"]

    # For x1
    number_contigs = 0
    length_contigs = 0
    depth_contigs = 0
    number_contigs_500 = 0
    for contig in contig_summary_yaml:
        if contig_summary_yaml[contig]['total_length'] > 500:
            number_contigs_500 += 1
        number_contigs += 1
        length_contigs += contig_summary_yaml[contig]["total_length"] 
        depth_contigs += contig_summary_yaml[contig]["total_depth"]
    denovo_assembly["summary"]["contigs"] = number_contigs
    denovo_assembly["summary"]["length"] = length_contigs
    denovo_assembly["summary"]["depth"] = float(depth_contigs/length_contigs)
    denovo_assembly['summary']['contigs_500'] = number_contigs_500

    # For x10
    number_contigs = 0
    length_contigs = 0
    depth_contigs = 0
    for contig in contig_summary_yaml:
        if contig_summary_yaml[contig]["coverage"] >= float(10):
            number_contigs += 1
            length_contigs += contig_summary_yaml[contig]["total_length"] 
            depth_contigs += contig_summary_yaml[contig]["total_depth"]
    mapping_qc["summary"]['values_at_floor_of_depth'] = {
        'x10': {
            'contigs': number_contigs,
            'length': length_contigs,
            'depth': float(depth_contigs/length_contigs)
        }
    }


def extract_bbuk_log(denovo_assembly: Category, results: Dict, component_name: str) -> None:
    file_name = "log/setup__filter_reads_with_bbduk.err.log"
    file_key = common.json_key_cleaner(file_name)
    file_path = os.path.join(component_name, file_name)
    reads_in = common.get_group_from_file("readsIn:\s*([0-9]+),", file_path)
    reads_in = int(reads_in) if reads_in != None else 0
    reads_removed = common.get_group_from_file("readsRemoved:\s*([0-9]+),", file_path)
    reads_removed = int(reads_removed) if reads_removed != None else 0
    denovo_assembly["summary"]["number_of_reads"] = reads_in
    denovo_assembly["summary"]["number_of_filtered_reads"] = reads_in - reads_removed 


def extract_quast_report(denovo_assembly: Category, results: Dict, component_name: str) -> None:
    file_name = "quast/report.tsv"
    file_key = common.json_key_cleaner(file_name)
    file_path = os.path.join(component_name, file_name)
    results[file_key] = {
        "N90": int(common.get_group_from_file("N90\t([0-9]+)", file_path)),
        "L50": int(common.get_group_from_file("L50\t([0-9]+)", file_path)),
        "L90": int(common.get_group_from_file("L90\t([0-9]+)", file_path))
    }
    denovo_assembly['summary']["GC"] = float(common.get_group_from_file("GC \(%\)\t([0-9]+[\.]?[0-9]*)", file_path))
    denovo_assembly['summary']["N50"] = int(common.get_group_from_file("N50\t([0-9]+)", file_path))


def extract_contig_variants(mapping_qc: Category, results: Dict, component_name: str) -> None:
    file_name = "contigs.variants"
    file_key = common.json_key_cleaner(file_name)
    file_path = os.path.join(component_name, file_name)
    yaml = common.get_yaml(file_path)
    mapping_qc["summary"]["snps"] = {
        'x10_10%':
            {
                'snps': yaml["variant_table"][9][9],
                'indels': yaml["indels"],
                'deletions': yaml["deletions"]
            }
    }


def extract_contig_stats(mapping_qc: Category, results: Dict, component_name: str) -> None:
    file_name = "contigs.stats"
    file_key = common.json_key_cleaner(file_name)
    file_path = os.path.join(component_name, file_name)
    results[file_key] = {}
    with open(file_path, "r") as fh:
        buffer = fh.readlines()
    for line in buffer:
        if line.startswith("SN"):
            temp = line.replace("SN\t","")
            key = temp.split(":")[0].strip()
            key = key.replace(" ","_")
            value = temp.split(":")[1].strip().split("\t")[0]
            results[file_key][key] = value
    mapping_qc["summary"]["mapped"] = {
        'reads_mapped': int(results[file_key]["reads_mapped"]),
        'reads_unmapped': int(results[file_key]["reads_unmapped"]),
        'insert_size_average': float(results[file_key]["insert_size_average"]),
        'insert_size_standard_deviation': float(results[file_key]["insert_size_standard_deviation"]),
    }


def save_contigs_locations(contigs: Category, results: Dict, component_name: str) -> None:
    file_name = "contigs.fasta"
    file_path = os.path.join(component_name, file_name)
    contigs["summary"]["data"] = file_path

def datadump(samplecomponent_ref_json: Dict):
    samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    sample = Sample.load(samplecomponent.sample)
    component = Component.load(samplecomponent.component)
    denovo_assembly = samplecomponent.get_category("denovo_assembly")
    if denovo_assembly is None:
        denovo_assembly = Category(value={
            "name": "denovo_assembly",
            "component": {"id": samplecomponent["component"]["_id"], "name": samplecomponent["component"]["name"]},
            "summary": {},
            "report": {}
        }
        )
    mapping_qc = samplecomponent.get_category("mapping_qc")
    if mapping_qc is None:
        mapping_qc = Category(value={
            "name": "mapping_qc",
            "component": {"id": samplecomponent["component"]["_id"], "name": samplecomponent["component"]["name"]},
            "summary": {},
            "report": {}
        }
        )
    contigs = samplecomponent.get_category("contigs")
    if contigs is None:
        contigs = Category(value={
            "name": "contigs",
            "component": {"id": samplecomponent["component"]["_id"], "name": samplecomponent["component"]["name"]},
            "summary": {},
            "report": {}
        }
        )
    extract_contigs_sum_cov(denovo_assembly, mapping_qc, samplecomponent["results"], samplecomponent["component"]["name"])
    extract_bbuk_log(denovo_assembly, samplecomponent["results"], samplecomponent["component"]["name"])
    extract_quast_report(denovo_assembly, samplecomponent["results"], samplecomponent["component"]["name"])
    extract_contig_variants(mapping_qc, samplecomponent["results"], samplecomponent["component"]["name"])
    extract_contig_stats(mapping_qc, samplecomponent["results"], samplecomponent["component"]["name"])
    save_contigs_locations(contigs, samplecomponent["results"], samplecomponent["component"]["name"])
    samplecomponent.set_category(denovo_assembly)
    samplecomponent.set_category(mapping_qc)
    samplecomponent.set_category(contigs)
    sample.set_category(denovo_assembly)
    sample.set_category(mapping_qc)
    sample.set_category(contigs)
    samplecomponent.save_files()
    common.set_status_and_save(sample, samplecomponent, "Success")
    with open(os.path.join(samplecomponent["component"]["name"], "datadump_complete"), "w+") as fh:
        fh.write("done")


datadump(
    snakemake.params.samplecomponent_ref_json,
)
