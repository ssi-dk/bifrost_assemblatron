from bifrostlib import common
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Category
from typing import Dict
import os


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

def save_contigs_locations(contigs: Category, results: Dict, component_name: str, sample_name: str) -> None:
    file_name = f"{sample_name}.fasta"
    file_path = os.path.join(os.getcwd(), component_name, file_name)
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
    contigs = samplecomponent.get_category("contigs")
    if contigs is None:
        contigs = Category(value={
            "name": "contigs",
            "component": {"id": samplecomponent["component"]["_id"], "name": samplecomponent["component"]["name"]},
            "summary": {},
            "report": {}
        }
        )
    extract_bbuk_log(denovo_assembly, samplecomponent["results"], samplecomponent["component"]["name"])
    save_contigs_locations(contigs, samplecomponent["results"], samplecomponent["component"]["name"], samplecomponent["sample"]["name"])
    samplecomponent.set_category(denovo_assembly)
    samplecomponent.set_category(contigs)
    sample.set_category(denovo_assembly)
    sample.set_category(contigs)
    samplecomponent.save_files()
    common.set_status_and_save(sample, samplecomponent, "Success")
    with open(os.path.join(samplecomponent["component"]["name"], "datadump_complete"), "w+") as fh:
        fh.write("done")


datadump(
    snakemake.params.samplecomponent_ref_json,
)
