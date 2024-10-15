import os
import shutil

import pytest
from contextlib import contextmanager
import importlib

import tests.resources.xtc
from looseversion import LooseVersion
import sys

from loguru import logger

if LooseVersion(sys.version) >= LooseVersion("3.9"):
    from importlib.resources import files
    def _get_path(parent, filename: str):
        return files(parent).joinpath(filename)
else:
    import inspect
    import pkg_resources
    def _get_path(parent, filename: str):
        parent_string = inspect.getmodule(parent).__name__
        return pkg_resources.resource_filename(parent_string, filename)


######################
# Handle Logging #####
######################

@pytest.fixture(autouse=True)
def disable_loguru_logging():
    # Remove standard logger if present
    try:
        logger.remove(None)
    except ValueError:
        pass
    try:
        logger.remove(0)
    except ValueError:
        pass

@contextmanager
def inject_pytest_logger():
    """
    Injects a loguru logger redirecting to stderr for capturing during pytest executions.
    Needs to be used after the test is initialized.
    """
    logger.enable('TRAMbio')
    ref = logger.add(sys.stderr, format="{level} | {message}", filter="TRAMbio")
    yield
    try:
        logger.remove(ref)
    except ValueError:
        pass
    logger.disable('TRAMbio')


##################
# Add marker #####
##################

def pytest_addoption(parser):
    parser.addoption(
        "--run-slow", action="store_true", default=False, help="run slow tests"
    )
    parser.addoption(
        "--skip-internet", action="store_true", default=False, help="disable test that use internet connection"
    )


_PACKAGE_REQUIREMENTS = [
    "MDAnalysis",
    "mdtraj",
    "lxml"
]
"""List of optional packages required by some tests."""

def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")
    config.addinivalue_line("markers", "internet: mark test as requiring internet connection")
    for package in _PACKAGE_REQUIREMENTS:
        config.addinivalue_line("markers", f"requires_{package}: mark test as requiring the package {package}")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--run-slow"):
        # --runslow not given in cli: skip slow tests
        skip_slow = pytest.mark.skip(reason="need --run-slow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)

    if config.getoption("--skip-internet"):
        skip_internet = pytest.mark.skip(reason="requires internet connection")
        for item in items:
            if "internet" in item.keywords:
                item.add_marker(skip_internet)

    unavailable_packages = []
    for package in _PACKAGE_REQUIREMENTS:
        try:
            importlib.import_module(package)
        except ModuleNotFoundError:
            unavailable_packages.append((f"requires_{package}", pytest.mark.skip(reason=f"requires package {package}")))
        except ImportError:
            unavailable_packages.append((f"requires_{package}", pytest.mark.skip(reason=f"requires package {package}")))

    if len(unavailable_packages) > 0:
        for item in items:
            for marker, skip_flag in unavailable_packages:
                if marker in item.keywords:
                    item.add_marker(skip_flag)


######################
# Resource Paths #####
######################

# PDB and XTC paths

@pytest.fixture
def path_pdb_test_empty():
    yield str(_get_path(tests.resources.pdb, "test_empty.pdb"))

@pytest.fixture
def path_pdb_test_single_model():
    yield str(_get_path(tests.resources.pdb, "test_model_1.pdb"))

@pytest.fixture
def path_pdb_test_multiple_models():
    yield str(_get_path(tests.resources.pdb, "test_model_2.pdb"))

@pytest.fixture
def path_pdb_test_hetatm():
    yield str(_get_path(tests.resources.pdb, "test_hetatm.pdb"))

@pytest.fixture
def path_pdb_sample():
    yield str(_get_path(tests.resources.xtc, "B2AR_sample_chain_e_f.pdb"))

@pytest.fixture
def path_trajectory_sample(tmp_path):
    pdb_name = "B2AR_sample_chain_e_f.pdb"
    xtc_name = "B2AR_sample_chain_e_f.xtc"
    xtc_resource = str(_get_path(tests.resources.xtc, xtc_name))

    tmp_resource = str(tmp_path / xtc_name)
    shutil.copyfile(xtc_resource, tmp_resource)

    yield str(_get_path(tests.resources.xtc, pdb_name)), tmp_resource

    os.remove(tmp_resource)

    # currently requires teardown code for MDAnalysis
    for file in [
        str(tmp_path / f".{xtc_name}_offsets.npz"),
        str(tmp_path / f".{xtc_name}_offsets.lock"),
    ]:
        if os.path.exists(file):
            os.remove(file)

@pytest.fixture
def path_trajectory_sample_spliced(tmp_path):
    pdb_name = "B2AR_sample_chain_e_f.pdb"
    xtc_name = "B2AR_sample_spliced_chain_e_f.xtc"
    xtc_resource = str(_get_path(tests.resources.xtc, xtc_name))

    tmp_resource = str(tmp_path / xtc_name)
    shutil.copyfile(xtc_resource, tmp_resource)

    yield str(_get_path(tests.resources.xtc, pdb_name)), tmp_resource

    os.remove(tmp_resource)

    # currently requires teardown code for MDAnalysis
    for file in [
        str(tmp_path / f".{xtc_name}_offsets.npz"),
        str(tmp_path / f".{xtc_name}_offsets.lock"),
    ]:
        if os.path.exists(file):
            os.remove(file)


# XML paths

@pytest.fixture
def path_xml_components_1():
    yield str(_get_path(tests.resources.xml, "components1.xml"))

@pytest.fixture
def path_xml_components_2():
    yield str(_get_path(tests.resources.xml, "components2.xml"))

@pytest.fixture
def path_xml_broken():
    yield str(_get_path(tests.resources.xml, "broken.xml"))

# extended XML paths

@pytest.fixture
def path_xml_components_for_pdb_test_single_model(path_pdb_test_single_model):
    yield (str(_get_path(tests.resources.xml, "test_model_1_components.xml")),
           path_pdb_test_single_model)

@pytest.fixture
def path_xml_components_for_trajectory_sample_spliced(path_trajectory_sample_spliced):
    yield (str(_get_path(tests.resources.xml, "sample_spliced_components.xml")),
           *path_trajectory_sample_spliced)

# BND paths

@pytest.fixture
def path_bnd_protein():
    yield str(_get_path(tests.resources.bnd, "example_protein.bnd"))

@pytest.fixture
def path_bnd_trajectory():
    yield str(_get_path(tests.resources.bnd, "example_trajectory.bnd"))

# GraphML paths

@pytest.fixture
def path_graphml_simple():
    yield str(_get_path(tests.resources.graphml, "simple.graphml"))


@pytest.fixture
def path_graphml_broken():
    yield str(_get_path(tests.resources.graphml, "broken.graphml"))


@pytest.fixture
def path_graphml_directed():
    yield str(_get_path(tests.resources.graphml, "directed.graphml"))
