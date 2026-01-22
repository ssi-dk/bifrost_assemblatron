from bifrostlib import common
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Category
from bifrostlib.database_interface import save_file
from typing import Dict
import os


# def extract_fastp_log(denovo_assembly: Category, results: Dict, component_name: str) -> None:
#     file_name = "log/setup__filter_reads_with_fastp.err.log"
#     file_key = common.json_key_cleaner(file_name)
#     file_path = os.path.join(component_name, file_name)
#     reads_in = common.get_group_from_file("Input:\s*([0-9]+) reads", file_path)
#     reads_in = int(reads_in) if reads_in != None else 0
#     reads_removed = common.get_group_from_file("Total Removed:\s*([0-9]+) reads", file_path)
#     reads_removed = int(reads_removed) if reads_removed != None else 0
#     denovo_assembly["summary"]["number_of_reads"] = reads_in
#     denovo_assembly["summary"]["number_of_filtered_reads"] = reads_in - reads_removed 

def save_contigs_location(contigs: Category, component_name: str, sample_name: str) -> None:
    file_name = f"{sample_name}.fasta"
    file_path = os.path.join(os.getcwd(), component_name, file_name)
    contigs["summary"]["data"] = file_path

def save_trimmed_reads_location(paired_reads: Category, component_name: str, sample_name: str) -> None:

    file_paths = [os.path.join(os.getcwd, f"{sample_name}.R1.trim.fastq.gz"),
                  os.path.join(os.getcwd, f"{sample_name}.R2.trim.fastq.gz")]
    paired_reads["summary"]["trimmed"] = file_paths

def save_contigs(contigs: Category, component_name: str, sample_name: str) -> None:
    file_path = contigs["summary"]["data"]
    _id = None
    _type = None
    file_id = save_file(_id, sample_name, _type, file_path)
    if file_id is not None:
        contigs["summary"]["file_id"] = file_id


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
    contigs = samplecomponent.get_category("contigs")
    if contigs is None:
        contigs = Category(value={
            "name": "contigs",
            "component": {"id": samplecomponent["component"]["_id"], "name": samplecomponent["component"]["name"]},
            "summary": {},
            "report": {}
        }
        )
    paired_reads = samplecomponent.get_category("paired_reads")
    if paired_reads is None:
        paired_reads = Category(value={
            "name": "paired_reads",
            "component": {"id": samplecomponent["component"]["_id"], "name": samplecomponent["component"]["name"]},
            "summary": {},
            "report": {}
        })
    save_trimmed_reads_location(paired_reads, samplecomponent)
    # extract_fastp_log(denovo_assembly, samplecomponent["results"], samplecomponent["component"]["name"])
    save_contigs_location(contigs, samplecomponent["component"]["name"], samplecomponent["sample"]["name"])
    save_contigs(contigs, samplecomponent["component"]["name"], samplecomponent["sample"]["name"])
    samplecomponent.set_category(paired_reads)
    samplecomponent.set_category(denovo_assembly)
    samplecomponent.set_category(contigs)
    sample.set_category(paired_reads)
    sample.set_category(denovo_assembly)
    sample.set_category(contigs)
    samplecomponent.save_files()
    common.set_status_and_save(sample, samplecomponent, "Success")
    with open(os.path.join(samplecomponent["component"]["name"], "datadump_complete"), "w+") as fh:
        fh.write("done")


datadump(
    snakemake.params.samplecomponent_ref_json,
)
