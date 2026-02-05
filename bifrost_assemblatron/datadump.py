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
    save_contigs_location(contigs, samplecomponent["component"]["name"], samplecomponent["sample"]["name"])
    save_contigs(contigs, samplecomponent["component"]["name"], samplecomponent["sample"]["name"])
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
