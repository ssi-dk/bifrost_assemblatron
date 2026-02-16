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
    # Passed contigs
    file_path = os.path.join(component_name, f"{sample_name}_trimmed_stat.tsv")
    with open(file_path) as fh:
        header = next(fh)
        (group, cov_threshold, min_contig_len, contig_count, sum_len, N50, avg_cov) = fh.readline().strip().split()
        contigs["summary"]["contigs"] = int(contig_count)
        contigs["summary"]["sum_len"] = int(sum_len)
        contigs["summary"]["N50"] = int(N50)
        contigs["summary"]["avg_cov"] = float(avg_cov)

    # Low coverage contigs
    file_path = os.path.join(component_name, f"{sample_name}_cov_fail_stat.tsv")
    with open(file_path) as fh:
        header = next(fh)
        (group, cov_threshold, min_contig_len, contig_count, sum_len, N50, avg_cov) = fh.readline().strip().split()
        contigs["summary"]["low_cov_len"] = int(sum_len)
        contigs["summary"]["low_cov_N50"] = int(N50)
        contigs["summary"]["low_cov_contigs"] = int(contig_count)
        contigs["summary"]["low_cov_avg_cov"] = float(avg_cov)   

def datadump(samplecomponent_ref_json: Dict):
    samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
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
    snakemake.params.samplecomponent_ref_json,
)
