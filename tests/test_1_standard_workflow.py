import os
import shutil
import pymongo
import pytest
import argparse
import yaml
import subprocess
from bson.objectid import ObjectId
from bifrost_run_launcher import launcher as brl
from bifrost_assemblatron import launcher as bal
# Test data from bifrost_test_data git repo

@pytest.fixture
def mydb():
    if os.getenv("BIFROST_DB_KEY", None) is not None:
        db_connection = pymongo.MongoClient(os.getenv("BIFROST_DB_KEY"))
        return db_connection.get_database()
    else:
        raise ValueError("BIFROST_DB_KEY not set")

def test_db_connection(mydb):
    mydb.list_collection_names()

def test_clear_db(mydb):
    os.chdir('/bifrost_test_data/')
    col_components = mydb["components"]
    col_samples = mydb["samples"]
    col_runs = mydb["runs"]
    col_sample_components = mydb["sample_components"]
    col_components.drop()
    col_samples.drop()
    col_runs.drop()
    col_sample_components.drop()


def test_install_component(mydb):
    test_clear_db(mydb)

    args: argparse.Namespace = brl.parser(["--install"])
    brl.run_program(args)

    args: argparse.Namespace = bal.parser(["--install"])
    bal.run_program(args)
    

def test_setup_of_data(mydb):
    test_install_component(mydb)
    if not (os.path.isfile("/bifrost_test_data/read_data/S1_R1.fastq.gz") and os.path.isfile("/bifrost_test_data/read_data/S1_R2.fastq.gz")):
        process = subprocess.Popen(
                "bash download_S1.sh",
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True,
                cwd="/bifrost_test_data/read_data/")
        process_out, process_err = process.communicate()
    if os.path.isdir("/bifrost_test_data/output/setup"):
        shutil.rmtree("/bifrost_test_data/output/setup")
    os.makedirs("/bifrost_test_data/output/setup")
    args = brl.parser([
        "-pre", "/bifrost_test_data/pre.sh",
        "-per", "/bifrost_test_data/per_sample.sh",
        "-post", "/bifrost_test_data/post.sh",
        "-meta", "/bifrost_test_data/run_metadata.tsv",
        "-reads", "/bifrost_test_data/read_data",
        "-name", "test_run",
        "-type", "test",
        "-out", "/bifrost_test_data/output/setup"
    ])
    brl.run_program(args)
    assert os.path.isfile("run.yaml")
    with open("run.yaml", "r") as file_handle:
        run_doc = yaml.safe_load(file_handle)
    print(run_doc["samples"])
    assert(run_doc["samples"][0]["name"] == "S1")
    assert os.path.isfile("/bifrost_test_data/output/setup/samples.yaml")
    assert os.path.isfile("/bifrost_test_data/output/setup/run_script.sh")

def test_pipeline(mydb):
    test_setup_of_data(mydb)
    with open("run.yaml", "r") as file_handle:
        run_doc = yaml.safe_load(file_handle)
    sample_id = run_doc["samples"][0]['_id'].strip("ObjectId('").strip("')")
    print(sample_id, str(sample_id))
    if os.path.isdir("/bifrost_test_data/output/test1"):
        shutil.rmtree("/bifrost_test_data/output/test1")
    os.makedirs("/bifrost_test_data/output/test1")
    args = bal.parser(["-id", sample_id,"-out", "/bifrost_test_data/output/test1"])
    bal.run_program(args)
    assert os.path.isfile("/bifrost_test_data/output/test1/assemblatron/datadump_complete")