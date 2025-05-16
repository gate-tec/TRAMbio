import importlib
import os
import pytest
import sys

import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ParseError

from TRAMbio.services import IOServiceRegistry
from TRAMbio.tram_xtc import main
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


#############################
# Environment Variables #####
#############################

class TestEnvironmentVariables:

    @pytest.fixture
    def env_vars_tram_xtc_default(self, monkeypatch):
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
    def env_vars_tram_xtc_no_interactions(self, monkeypatch):
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
    def args_sample_mdanalysis(self, tmp_path, path_trajectory_sample):
        pdb_path, xtc_path = path_trajectory_sample

        orig_argv = list(sys.argv)
        d = tmp_path / "xtc_sample"
        sys.argv = [
            "tram-xtc",
            "--pdb", pdb_path,
            "--xtc", xtc_path,
            "--out-dir", str(d),
            "--stride", "1",
            "--edges",
            "--module", "MDAnalysis",
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
    def args_sample_spliced_mdanalysis(self, tmp_path, path_trajectory_sample_spliced):
        pdb_path, xtc_path = path_trajectory_sample_spliced

        orig_argv = list(sys.argv)
        d = tmp_path / "xtc_sample"
        sys.argv = [
            "tram-xtc",
            "--pdb", pdb_path,
            "--xtc", xtc_path,
            "--out-dir", str(d),
            "--stride", "1",
            "--module", "MDAnalysis",
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
    def args_sample_mdtraj(self, tmp_path, path_trajectory_sample):
        pdb_path, xtc_path = path_trajectory_sample

        orig_argv = list(sys.argv)
        d = tmp_path / "xtc_sample"
        sys.argv = [
            "tram-xtc",
            "--pdb", pdb_path,
            "--xtc", xtc_path,
            "--out-dir", str(d),
            "--stride", "1",
            "--edges",
            "--module", "mdtraj",
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
    def args_sample_spliced_mdtraj(self, tmp_path, path_trajectory_sample_spliced):
        pdb_path, xtc_path = path_trajectory_sample_spliced

        orig_argv = list(sys.argv)
        d = tmp_path / "xtc_sample"
        sys.argv = [
            "tram-xtc",
            "--pdb", pdb_path,
            "--xtc", xtc_path,
            "--out-dir", str(d),
            "--stride", "1",
            "--module", "mdtraj",
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
class TestTramXtc(TestArguments, TestEnvironmentVariables):

    @pytest.mark.requires_MDAnalysis
    def test_tram_xtc_sample_mdanalysis(self, args_sample_mdanalysis, env_vars_tram_xtc_default):
        main()
        assert len(list(args_sample_mdanalysis.iterdir())) == 2  # components XML & edges BND

        xml_f = str(args_sample_mdanalysis / "sample_components.xml")
        assert os.path.exists(xml_f)
        assert IOServiceRegistry.XML.single_service().validate_xml(xml_f)

        with open(xml_f) as file:
            xml_lines = file.readlines()

        # verify XML structure
        try:
            components_root = ET.fromstringlist(xml_lines)
        except ParseError:
            pytest.fail(reason="Unable to parse constructed XML file.")

        # verify copying of states from temp file
        parsed_states = components_root.find("{tram:components}states")  # noqa components_root is always set
        assert len(parsed_states) == 22

    @pytest.mark.requires_MDAnalysis
    def test_tram_xtc_sample_no_bonds_mdanalysis(self, args_sample_mdanalysis, env_vars_tram_xtc_no_interactions):
        main()
        assert len(list(args_sample_mdanalysis.iterdir())) == 2  # components XML & edges BND

        xml_f = str(args_sample_mdanalysis / "sample_components.xml")
        assert os.path.exists(xml_f)
        assert IOServiceRegistry.XML.single_service().validate_xml(xml_f)

        with open(xml_f) as file:
            xml_lines = file.readlines()

        # verify XML structure
        try:
            components_root = ET.fromstringlist(xml_lines)
        except ParseError:
            pytest.fail(reason="Unable to parse constructed XML file.")

        # verify copying of states from temp file
        parsed_states = components_root.find("{tram:components}states")  # noqa components_root is always set
        assert len(parsed_states) == 1

    @pytest.mark.requires_MDAnalysis
    def test_tram_xtc_sample_spliced_mdanalysis(self, args_sample_spliced_mdanalysis, env_vars_tram_xtc_default):
        main()
        assert len(list(args_sample_spliced_mdanalysis.iterdir())) == 1  # components XML

        xml_f = str(args_sample_spliced_mdanalysis / "sample_components.xml")
        assert os.path.exists(xml_f)
        assert IOServiceRegistry.XML.single_service().validate_xml(xml_f)

        with open(xml_f) as file:
            xml_lines = file.readlines()

        # verify XML structure
        try:
            components_root = ET.fromstringlist(xml_lines)
        except ParseError:
            pytest.fail(reason="Unable to parse constructed XML file.")

        # verify copying of states from temp file
        parsed_states = components_root.find("{tram:components}states")  # noqa components_root is always set
        assert len(parsed_states) == 3

    @pytest.mark.requires_mdtraj
    def test_tram_xtc_sample_mdtraj(self, args_sample_mdtraj, env_vars_tram_xtc_default):
        main()
        assert len(list(args_sample_mdtraj.iterdir())) == 2  # components XML & edges BND

        xml_f = str(args_sample_mdtraj / "sample_components.xml")
        assert os.path.exists(xml_f)
        assert IOServiceRegistry.XML.single_service().validate_xml(xml_f)

        with open(xml_f) as file:
            xml_lines = file.readlines()

        # verify XML structure
        try:
            components_root = ET.fromstringlist(xml_lines)
        except ParseError:
            pytest.fail(reason="Unable to parse constructed XML file.")

        # verify copying of states from temp file
        parsed_states = components_root.find("{tram:components}states")  # noqa components_root is always set
        assert len(parsed_states) == 22

    @pytest.mark.requires_mdtraj
    def test_tram_xtc_sample_no_bonds_mdtraj(self, args_sample_mdtraj, env_vars_tram_xtc_no_interactions):
        main()
        assert len(list(args_sample_mdtraj.iterdir())) == 2  # components XML & edges BND

        xml_f = str(args_sample_mdtraj / "sample_components.xml")
        assert os.path.exists(xml_f)
        assert IOServiceRegistry.XML.single_service().validate_xml(xml_f)

        with open(xml_f) as file:
            xml_lines = file.readlines()

        # verify XML structure
        try:
            components_root = ET.fromstringlist(xml_lines)
        except ParseError:
            pytest.fail(reason="Unable to parse constructed XML file.")

        # verify copying of states from temp file
        parsed_states = components_root.find("{tram:components}states")  # noqa components_root is always set
        assert len(parsed_states) == 1

    @pytest.mark.requires_mdtraj
    def test_tram_xtc_sample_spliced_mdtraj(self, args_sample_spliced_mdtraj, env_vars_tram_xtc_default):
        main()
        assert len(list(args_sample_spliced_mdtraj.iterdir())) == 1  # components XML

        xml_f = str(args_sample_spliced_mdtraj / "sample_components.xml")
        assert os.path.exists(xml_f)
        assert IOServiceRegistry.XML.single_service().validate_xml(xml_f)

        with open(xml_f) as file:
            xml_lines = file.readlines()

        # verify XML structure
        try:
            components_root = ET.fromstringlist(xml_lines)
        except ParseError:
            pytest.fail(reason="Unable to parse constructed XML file.")

        # verify copying of states from temp file
        parsed_states = components_root.find("{tram:components}states")  # noqa components_root is always set
        assert len(parsed_states) == 3
