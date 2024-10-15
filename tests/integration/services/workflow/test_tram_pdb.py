import importlib
import os
import re

import pytest
import sys

from TRAMbio.services import IOServiceRegistry
from TRAMbio.tram_pdb import main
from TRAMbio.services.parameter import (
    AromaticInteractionParameter,
    HydrophobicInteractionParameter,
    HydrogenBondParameter,
    DisulphideBridgeParameter,
    CationPiInteractionParameter,
    PdbEntryInteractionParameter
)
import TRAMbio.services.interactions.registry as tram_interaction_registry
import TRAMbio.services.structure.registry as tram_structure_registry
import TRAMbio.services.workflow.registry as tram_workflow_registry
from tests.conftest import inject_pytest_logger


#############################
# Environment Variables #####
#############################

class TestEnvironmentVariables:

    @pytest.fixture
    def env_vars_tram_pdb_default(self, monkeypatch):
        monkeypatch.setenv(HydrogenBondParameter.INCLUDE.value, "true")
        monkeypatch.setenv(HydrogenBondParameter.ENERGY_THRESHOLD.value, "-0.1")
        monkeypatch.setenv(HydrogenBondParameter.CUTOFF_DISTANCE.value, "3.0")
        monkeypatch.setenv(HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value, "0.0")
        monkeypatch.setenv(HydrogenBondParameter.BAR_COUNT.value, "5")
        monkeypatch.setenv(HydrophobicInteractionParameter.INCLUDE.value, "true")
        monkeypatch.setenv(DisulphideBridgeParameter.INCLUDE.value, "true")
        monkeypatch.setenv(AromaticInteractionParameter.INCLUDE.value, "true")
        monkeypatch.setenv(CationPiInteractionParameter.INCLUDE.value, "true")
        monkeypatch.setenv(PdbEntryInteractionParameter.LINK_INCLUDE.value, "true")
        monkeypatch.setenv(PdbEntryInteractionParameter.SSBOND_INCLUDE.value, "true")
        monkeypatch.setenv(PdbEntryInteractionParameter.CONECT_INCLUDE.value, "true")

        # force reload registry to load new environment variables
        importlib.reload(tram_interaction_registry)
        importlib.reload(tram_structure_registry)
        importlib.reload(tram_workflow_registry)

    @pytest.fixture
    def env_vars_tram_pdb_no_interactions(self, monkeypatch):
        monkeypatch.setenv(HydrogenBondParameter.INCLUDE.value, "false")
        monkeypatch.setenv(HydrophobicInteractionParameter.INCLUDE.value, "false")
        monkeypatch.setenv(DisulphideBridgeParameter.INCLUDE.value, "false")
        monkeypatch.setenv(AromaticInteractionParameter.INCLUDE.value, "false")
        monkeypatch.setenv(CationPiInteractionParameter.INCLUDE.value, "false")
        monkeypatch.setenv(PdbEntryInteractionParameter.LINK_INCLUDE.value, "false")
        monkeypatch.setenv(PdbEntryInteractionParameter.SSBOND_INCLUDE.value, "false")
        monkeypatch.setenv(PdbEntryInteractionParameter.CONECT_INCLUDE.value, "false")

        # force reload registry to load new environment variables
        importlib.reload(tram_interaction_registry)
        importlib.reload(tram_structure_registry)
        importlib.reload(tram_workflow_registry)


#################
# Arguments #####
#################

class TestArguments:

    @pytest.fixture
    def args_sample_with_edges(self, tmp_path, path_pdb_sample):
        orig_argv = list(sys.argv)
        d = tmp_path / "pdb_sample"
        sys.argv = [
            "tram-pdb",
            "--pdb", str(path_pdb_sample),
            "--out-dir", str(d),
            "--edges",
            "--name", "sample",
            "--log-level", "NONE"
        ]

        yield d

        sys.argv = orig_argv

        # clean-up files
        for file in os.listdir(d):
            os.remove(d / file)
        os.rmdir(d)

    @pytest.fixture
    def args_1l2y(self, tmp_path):
        orig_argv = list(sys.argv)
        d = tmp_path / "pdb_sample"
        sys.argv = [
            "tram-pdb",
            "--pdb", "1L2Y",
            "--out-dir", str(d),
            "--name", "sample",
            "--log-level", "NONE"
        ]

        yield d

        sys.argv = orig_argv

        # clean-up files
        for file in os.listdir(d):
            os.remove(d / file)
        os.rmdir(d)

    # arguments for logging
    @pytest.fixture
    def args_empty_pdb_info(self, tmp_path, path_pdb_test_empty):
        orig_argv = list(sys.argv)
        d = tmp_path / "pdb_sample"
        sys.argv = [
            "tram-pdb",
            "--pdb", str(path_pdb_test_empty),
            "--out-dir", str(d),
            "--name", "sample",
            "--log-level", "INFO"
        ]

        yield d

        sys.argv = orig_argv

        # clean-up files
        for file in os.listdir(d):
            os.remove(d / file)
        os.rmdir(d)

    @pytest.fixture
    def args_empty_pdb_trace(self, tmp_path, path_pdb_test_empty):
        orig_argv = list(sys.argv)
        d = tmp_path / "pdb_sample"
        sys.argv = [
            "tram-pdb",
            "--pdb", str(path_pdb_test_empty),
            "--out-dir", str(d),
            "--name", "sample",
            "--log-level", "TRACE"
        ]

        yield d

        sys.argv = orig_argv

        # clean-up files
        for file in os.listdir(d):
            os.remove(d / file)
        os.rmdir(d)

    @pytest.fixture
    def args_sample_pdb_trace(self, tmp_path, path_pdb_sample):
        orig_argv = list(sys.argv)
        d = tmp_path / "pdb_sample"
        sys.argv = [
            "tram-pdb",
            "--pdb", str(path_pdb_sample),
            "--out-dir", str(d),
            "--name", "sample",
            "--log-level", "TRACE"
        ]

        yield d

        sys.argv = orig_argv

        # clean-up files
        for file in os.listdir(d):
            os.remove(d / file)
        os.rmdir(d)


#####################
# General Tests #####
#####################

@pytest.mark.slow
@pytest.mark.requires_lxml
class TestTramPdb(TestArguments, TestEnvironmentVariables):

    def test_tram_pdb_sample(self, args_sample_with_edges, env_vars_tram_pdb_default):
        main()
        assert len(list(args_sample_with_edges.iterdir())) == 2  # components XML & edges BND

        xml_f = args_sample_with_edges / "sample_components.xml"
        assert IOServiceRegistry.XML.single_service().validate_xml(str(xml_f))

    def test_tram_pdb_sample_no_bonds(self, args_sample_with_edges, env_vars_tram_pdb_no_interactions):
        main()
        assert len(list(args_sample_with_edges.iterdir())) == 2  # components XML & edges BND

        xml_f = args_sample_with_edges / "sample_components.xml"
        assert IOServiceRegistry.XML.single_service().validate_xml(str(xml_f))

    def test_tram_pdb_1l2y(self, args_1l2y, env_vars_tram_pdb_default):
        main()
        assert len(list(args_1l2y.iterdir())) == 1  # components XML

        xml_f = args_1l2y / "sample_components.xml"
        assert IOServiceRegistry.XML.single_service().validate_xml(str(xml_f))


#####################################
# Test Error Handling & Logging #####
#####################################

@pytest.mark.slow
class TestTramPdbLogging(TestArguments, TestEnvironmentVariables):

    # evaluate no error without TRACE option (on faulty input)
    def test_no_error_without_trace(self, args_empty_pdb_info, capsys, env_vars_tram_pdb_default):
        with inject_pytest_logger():
            main()

        captured = capsys.readouterr()[0]

        critical_line = next(filter(lambda x: "CRITICAL" in x, str(captured).split("\n")), None)
        assert critical_line is not None
        assert re.search(r"No hydrogen atoms", critical_line)

    # evaluate error with TRACE option (on faulty input)
    def test_error_with_trace(self, args_empty_pdb_trace, env_vars_tram_pdb_default):
        with pytest.raises(ValueError, match="No hydrogen atoms"):
            main()

    # evaluate no error with TRACE (on correct input)
    def test_no_error_with_trace(self, args_sample_pdb_trace, env_vars_tram_pdb_default):
        main()
