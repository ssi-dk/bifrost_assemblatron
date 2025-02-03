from bifrostlib import common
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Category
from typing import Dict
import os
from Bio import SeqIO
import hashlib
from datetime import datetime

def extract_bbuk_log(denovo_assembly: Category, results: Dict, component_name: str) -> None:
    file_name = "log/setup__filter_reads_with_bbduk.err.log"
    file_key = common.json_key_cleaner(file_name)
    file_path = os.path.join(component_name, file_name)
    reads_in = common.get_group_from_file("Input:\s*([0-9]+) reads", file_path)
    reads_in = int(reads_in) if reads_in != None else 0
    reads_removed = common.get_group_from_file("Total Removed:\s*([0-9]+) reads", file_path)
    reads_removed = int(reads_removed) if reads_removed != None else 0
    denovo_assembly["summary"]["number_of_reads"] = reads_in
    denovo_assembly["summary"]["number_of_filtered_reads"] = reads_in - reads_removed

def save_contigs_locations(contigs: Category, results: Dict, component_name: str, sample_name: str) -> None:
    file_name = f"{sample_name}.fasta"
    file_path = os.path.join(os.getcwd(), component_name, file_name)
    contigs["summary"]["data"] = file_path

def calculate_md5(sequence: str) -> str:
    md5_hash = hashlib.md5(sequence.encode('utf-8')).hexdigest()
    return md5_hash

def save_contigs_data(contigs: Category, results: Dict, component_name: str, sample_name: str, verbose: int = 0) -> None:
    file_path = os.path.join(os.getcwd(), component_name, f"{sample_name}.fasta")
    print(f"inside save_contigs_data")
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"FASTA file not found at: {file_path}")

    print(f"file path is {file_path}")

    contig_data = {}
    contig_lengths = []
    gc_contents = []
    date = datetime.now().strftime('%Y-%m-%d')

    with open(file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            contig_name = record.id  # Contig header
            print(f"{contig_name}")
            contig_seq = str(record.seq).replace("\n", "")  # Sequence without newlines
            contig_data[contig_name] = contig_seq
            contig_length = len(contig_seq)
            gc_content = round((sum(contig_seq.count(x) for x in "GCgc") / contig_length) * 100, 2)
            contig_lengths.append(contig_length)
            gc_contents.append(gc_content)

            if verbose:
                print(f"Saving contig data for sample: {sample_name}, component: {component_name}, date {date}")
                print(f"Contig Name: {contig_name}")
                print(f"Length: {contig_length}")
                print(f"GC Content: {gc_content}%")
                #print(f"Sequence: {contig_seq[:10]}...")
                print("-" * 50)

    fasta_md5 = calculate_md5("".join(contig_data.values()))

    contigs["summary"]["md5"] = fasta_md5
    contigs["summary"]["num_contigs"] = len(contig_data.keys())
    contigs["summary"]["total_length"] = contig_lengths
    contigs["summary"]["gc_contents"] = gc_contents
    contigs["summary"]["date_added"] = date

    """
    contigs["summary"] = {
        "data": [{
            "path": file_path,
            "fasta_md5": fasta_md5,
            "num_contigs": len(contig_data),
            "total_length": contig_lengths,
            "gc_content": gc_contents,
            "date_added": date
        }]
    }
    """

def datadump(samplecomponent_ref_json: Dict):
    print("INSIDE DATADUMP FUNCTION I HOPE THIS WORK")
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
    save_contigs_data(contigs, samplecomponent["results"], samplecomponent["component"]["name"], samplecomponent["sample"]["name"],0)

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
