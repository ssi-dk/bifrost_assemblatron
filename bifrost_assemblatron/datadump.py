from bifrostlib import common
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Category
from bifrostlib.database_interface import save_file
from typing import Dict
import os


def save_contigs_location(contigs: Category, component_name: str, sample_name: str) -> None:
    file_name = f"{sample_name}_trimmed.fasta"
    file_path = os.path.join(os.getcwd(), component_name, file_name)
    contigs["summary"]["data"] = file_path
    contigs["summary"]["type"] = "assembly"


def save_contigs(contigs: Category, component_name: str, sample_name: str) -> None:
    file_path = contigs["summary"]["data"]
    _id = None
    _type = None
    file_id = save_file(_id, sample_name, _type, file_path)
    if file_id is not None:
        contigs["summary"]["file_id"] = file_id

def extract_assembly_statistics(contigs: Category, component_name: str, sample_name: str) -> None:
    # Passed contigs 10x< and 500bp<
    file_path = os.path.join(component_name, f"{sample_name}_above10x_stat.tsv")
    with open(file_path) as fh:
        header = next(fh)
        (group, cov_threshold, min_contig_len, contig_count, sum_len, N50, avg_cov) = fh.readline().strip().split()
        contigs["summary"]["contigs_10x"] = int(contig_count)
        contigs["summary"]["length_10x"] = int(sum_len)
        contigs["summary"]["N50_10x"] = int(N50)
        contigs["summary"]["coverage_10x"] = float(avg_cov)

    # low coverage 1x<10x and 500bp
    file_path = os.path.join(component_name, f"{sample_name}_above1x_stat.tsv")
    with open(file_path) as fh:
        header = next(fh)
        (group, cov_threshold, min_contig_len, contig_count, sum_len, N50, avg_cov) = fh.readline().strip().split()
        contigs["summary"]["contigs_1x"] = int(contig_count)
        contigs["summary"]["length_1x"] = int(sum_len)
        contigs["summary"]["N50_1x"] = int(N50)
        contigs["summary"]["coverage_1x"] = float(avg_cov)

    # Low coverage contigs <1x
    file_path = os.path.join(component_name, f"{sample_name}_below1x_stat.tsv")
    with open(file_path) as fh:
        header = next(fh)
        (group, cov_threshold, min_contig_len, contig_count, sum_len, N50, avg_cov) = fh.readline().strip().split()
        contigs["summary"]["contigs_0x"] = int(contig_count)
        contigs["summary"]["length_0x"] = int(sum_len)
        contigs["summary"]["N50_0x"] = int(N50)
        contigs["summary"]["coverage_0x"] = float(avg_cov)   
    
def datadump(samplecomponent_id: str):
    #samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
    samplecomponent_ref = SampleComponentReference(_id=samplecomponent_id)
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    sample = Sample.load(samplecomponent.sample)
    component = Component.load(samplecomponent.component)
    contigs = samplecomponent.get_category("contigs")
    if contigs is None:
        contigs = Category(value={
            "name": "contigs",
            "component": {"id": samplecomponent["component"]["_id"], "name": samplecomponent["component"]["name"]},
            "summary": {},
            "report": {}
        }
        )
    save_contigs_location(contigs, samplecomponent["component"]["name"], samplecomponent["sample"]["name"])
    extract_assembly_statistics(contigs, samplecomponent["component"]["name"], samplecomponent["sample"]["name"])
    save_contigs(contigs, samplecomponent["component"]["name"], samplecomponent["sample"]["name"])
    samplecomponent.set_category(contigs)
    sample.set_category(contigs)
    samplecomponent.save_files()
    common.set_status_and_save(sample, samplecomponent, "Success")
    with open(os.path.join(samplecomponent["component"]["name"], "datadump_complete"), "w+") as fh:
        fh.write("done")


datadump(
    snakemake.params.samplecomponent_id
)
