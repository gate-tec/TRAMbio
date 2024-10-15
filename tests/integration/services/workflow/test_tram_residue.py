import os

import pytest
import sys

import json

from TRAMbio.tram_residue import main


#################
# Arguments #####
#################

class TestArguments:

    # TODO: test options "key" and "max_states"

    @pytest.fixture
    def args_test_components_1(self, tmp_path, path_xml_components_1):
        orig_argv = list(sys.argv)
        d = tmp_path / "pdb_sample"
        sys.argv = [
            "tram-residue",
            "--xml", str(path_xml_components_1),
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
    def args_test_model_1(self, tmp_path, path_xml_components_for_pdb_test_single_model):
        xml_path, _ = path_xml_components_for_pdb_test_single_model
        orig_argv = list(sys.argv)
        d = tmp_path / "pdb_sample"
        sys.argv = [
            "tram-residue",
            "--xml", str(xml_path),
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
    def args_sample_sliced(self, tmp_path, path_xml_components_for_trajectory_sample_spliced):
        xml_path, _, _ = path_xml_components_for_trajectory_sample_spliced

        orig_argv = list(sys.argv)
        d = tmp_path / "pdb_sample"
        sys.argv = [
            "tram-residue",
            "--xml", str(xml_path),
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
    def args_sample_sliced_with_pdb(self, tmp_path, path_xml_components_for_trajectory_sample_spliced):
        xml_path, pdb_path, _ = path_xml_components_for_trajectory_sample_spliced

        orig_argv = list(sys.argv)
        d = tmp_path / "pdb_sample"
        sys.argv = [
            "tram-residue",
            "--xml", str(xml_path),
            "--pdb", str(pdb_path),
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

class TestTramResidue(TestArguments):

    def test_tram_residue_components_1(self, args_test_components_1):
        main()

        assert len(list(args_test_components_1.iterdir())) == 1  # JSON file

        with open(str(args_test_components_1 / "sample_residue_components.json"), "r") as file:
            data = json.load(file)

        assert len(data) == 1  # reduced to a single state
        assert "1" in data.keys()  # first occurrence of state is key
        assert len(data["1"]) == 0  # no rigid CA centers in components

    def test_tram_residue_model_1(self, args_test_model_1):
        main()

        assert len(list(args_test_model_1.iterdir())) == 1  # JSON file

        with open(str(args_test_model_1 / "sample_residue_components.json"), "r") as file:
            data = json.load(file)

        assert len(data) == 1  # reduced to a single state
        assert "-INF" in data.keys()  # first occurrence of state is key
        assert len(data["-INF"]) == 0  # no rigid CA centers in components

    def test_tram_residue_sample_spliced(self, args_sample_sliced):
        main()

        assert len(list(args_sample_sliced.iterdir())) == 1  # JSON file

        with open(str(args_sample_sliced / "sample_residue_components.json"), "r") as file:
            data = json.load(file)

        assert len(data) == 3  # all three states present
        assert "0" in data.keys()
        assert len(data["0"]) == 2  # two main components

    def test_tram_residue_sample_spliced_pdb(self, args_sample_sliced_with_pdb):
        main()

        assert len(list(args_sample_sliced_with_pdb.iterdir())) == 1  # JSON file

        with open(str(args_sample_sliced_with_pdb / "sample_residue_components.json"), "r") as file:
            data = json.load(file)

        assert len(data) == 3  # all three states present
        assert "0" in data.keys()
        assert len(data["0"]) == 2  # two main components