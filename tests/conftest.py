import os
import shutil
import tempfile
from pathlib import Path

import pymongo
import pytest
from bifrost_assemblatron import launcher
from bifrostlib import common, database_interface, datahandling
from bifrostlib.datahandling import (Component, ComponentReference, Run,
                                     RunReference, Sample, SampleReference)

bifrost_install_dir = os.environ['BIFROST_INSTALL_DIR']
bifrost_config_and_data_path = Path(f"{bifrost_install_dir}/bifrost/test_data")

assert datahandling.has_a_database_connection()
assert "TEST" in os.environ['BIFROST_DB_KEY'].upper()  # A very basic piece of protection ensuring the word test is in the DB


@pytest.fixture
def db():
    client = database_interface.get_connection()
    assert datahandling.has_a_database_connection()
    yield client.get_database()
    client.close()

@pytest.fixture
def initialized_launcher():
    launcher.initialize()

@pytest.fixture
def use_collection(db):
    '''Factory fixture returning a function to clean collections used in the test, and take care of their clean up after use.'''
    created_collections = []
    def _collection(name: str):
        created_collections.append(name)
        #db.drop_collection(name)
        return db[name]
    yield _collection
    # Clean up the used collections
    for name in created_collections:
        db.drop_collection(name)

@pytest.fixture(scope="module")
def clean_dir(request):
    #parent_dir = getattr(request.module, "test_dir", None)
    newpath = tempfile.mkdtemp(dir=bifrost_config_and_data_path/"output")
    yield newpath
    shutil.rmtree(newpath)
