import os
import shutil
import sys
import tempfile
from argparse import Namespace
from pathlib import Path

#import pymongo
import pytest
from bifrost_assemblatron import launcher
from bifrostlib import common, database_interface, datahandling
from bifrostlib.datahandling import (Component, ComponentReference, Run,
                                     RunReference, Sample, SampleReference)


bifrost_install_dir = os.environ['BIFROST_INSTALL_DIR']
bifrost_config_and_data_path = Path(f"{bifrost_install_dir}/bifrost/test_data")

@pytest.fixture
def sample(use_collection):
    json_entries = [
        {
            "_id": {"$oid": "000000000000000000000001"},
            "display_name": "S1",
            "name": "S1",
            "components": [],
            "categories": {
                "paired_reads": {
                    "summary": {
                        "data": [f"{bifrost_config_and_data_path}/samples/S1_R1.fastq.gz",
                                 f"{bifrost_config_and_data_path}/samples/S1_R2.fastq.gz"]
                    }
                }
            }
        }
    ]
    bson_entries = [database_interface.json_to_bson(i) for i in json_entries]
    use_collection("samples")
    sample = Sample(value=json_entries[0])
    sample.save()
    return sample


class TestBifrostAssemblatron:
    component_name = "assemblatron__v2.3.2"
    current_dir = os.getcwd()
    json_entries = [
        {
            "_id": {"$oid": "000000000000000000000001"},
            "display_name": "S1",
            "name": "S1",
            "components": [],
            "categories": {
                "paired_reads": {
                    "summary": {
                        "data": [f"{bifrost_config_and_data_path}/samples/S1_R1.fastq.gz",
                                 f"{bifrost_config_and_data_path}/samples/S1_R2.fastq.gz"]
                    }
                }
            }
        }
    ]
    bson_entries = [database_interface.json_to_bson(i) for i in json_entries]

    def test_info(self):
        launcher.run_pipeline(["--info"])

    def test_help(self):
        launcher.run_pipeline(["--help"])

    def test_pipeline(self, sample, use_collection, clean_dir):
        samples = use_collection("samples")
        sample_data = list(samples.find({}))
        test_args = [
            "--sample_name", "S1",
            "--outdir", clean_dir
        ]
        launcher.main(args=test_args)
        assert os.path.isfile(f"{clean_dir}/{self.component_name}/datadump_complete")
        sample_data = list(samples.find({}))
        db_contigs_path = sample_data[0]['categories']['contigs']['summary']['data']
        assert len(sample_data) == 1
        assert db_contigs_path == clean_dir + '/S1.fasta'
        assert os.path.isfile(db_contigs_path)
