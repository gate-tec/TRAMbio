import os

import pytest
import sys

from TRAMbio.tram_pymol import main


#################
# Arguments #####
#################

class TestArguments:

    # TODO: test options "key" and "max_states"

    @pytest.fixture
    def args_pdb_to_pymol(self, tmp_path, path_xml_components_for_pdb_test_single_model):
        xml_path, pdb_path = path_xml_components_for_pdb_test_single_model

        orig_argv = list(sys.argv)
        d = tmp_path / "pymol_sample"
        sys.argv = [
            "tram-pymol",
            "--pdb", str(pdb_path),
            "--xml", str(xml_path),
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

    @pytest.fixture
    def args_xtc_to_pymol_mdanalysis(self, tmp_path, path_xml_components_for_trajectory_sample_spliced):
        xml_path, pdb_path, xtc_path = path_xml_components_for_trajectory_sample_spliced

        orig_argv = list(sys.argv)
        d = tmp_path / "pymol_sample"
        sys.argv = [
            "tram-pymol",
            "--pdb", str(pdb_path),
            "--xml", str(xml_path),
            "--xtc", str(xtc_path),
            "--out-dir", str(d),
            "--module", "MDAnalysis",
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
    def args_xtc_to_pymol_mdtraj(self, tmp_path, path_xml_components_for_trajectory_sample_spliced):
        xml_path, pdb_path, xtc_path = path_xml_components_for_trajectory_sample_spliced

        orig_argv = list(sys.argv)
        d = tmp_path / "pymol_sample"
        sys.argv = [
            "tram-pymol",
            "--pdb", str(pdb_path),
            "--xml", str(xml_path),
            "--xtc", str(xtc_path),
            "--out-dir", str(d),
            "--module", "mdtraj",
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
    def args_pdb_to_pymol_trace(self, tmp_path, path_xml_components_for_pdb_test_single_model):
        xml_path, pdb_path = path_xml_components_for_pdb_test_single_model

        orig_argv = list(sys.argv)
        d = tmp_path / "pymol_sample"
        sys.argv = [
            "tram-pymol",
            "--pdb", str(pdb_path),
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
    def args_xtc_to_pymol_mdanalysis_trace(self, tmp_path, path_xml_components_for_trajectory_sample_spliced):
        xml_path, pdb_path, xtc_path = path_xml_components_for_trajectory_sample_spliced

        orig_argv = list(sys.argv)
        d = tmp_path / "pymol_sample"
        sys.argv = [
            "tram-pymol",
            "--pdb", str(pdb_path),
            "--xml", str(xml_path),
            "--xtc", str(xtc_path),
            "--out-dir", str(d),
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
    def args_incompatible_components_mdanalysis_trace(self, tmp_path, path_xml_components_for_trajectory_sample_spliced, path_xml_components_1):
        _, pdb_path, xtc_path = path_xml_components_for_trajectory_sample_spliced
        xml_path = path_xml_components_1  # use wrong component path

        orig_argv = list(sys.argv)
        d = tmp_path / "pymol_sample"
        sys.argv = [
            "tram-pymol",
            "--pdb", str(pdb_path),
            "--xml", str(xml_path),
            "--xtc", str(xtc_path),
            "--out-dir", str(d),
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
    def args_xtc_to_pymol_mdtraj_trace(self, tmp_path, path_xml_components_for_trajectory_sample_spliced):
        xml_path, pdb_path, xtc_path = path_xml_components_for_trajectory_sample_spliced

        orig_argv = list(sys.argv)
        d = tmp_path / "pymol_sample"
        sys.argv = [
            "tram-pymol",
            "--pdb", str(pdb_path),
            "--xml", str(xml_path),
            "--xtc", str(xtc_path),
            "--out-dir", str(d),
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


#############
# Tests #####
#############

class TestTramPyMol(TestArguments):

    def test_tram_pymol_pdb(self, args_pdb_to_pymol):
        main()

        assert len(list(args_pdb_to_pymol.iterdir())) == 2  # annotated PDB and .pml file

    @pytest.mark.requires_MDAnalysis
    def test_tram_pymol_xtc_mdanalysis(self, args_xtc_to_pymol_mdanalysis):
        main()

        assert len(list(args_xtc_to_pymol_mdanalysis.iterdir())) == 2  # annotated PDB and .pml file

    @pytest.mark.requires_mdtraj
    def test_tram_pymol_xtc_mdtraj(self, args_xtc_to_pymol_mdtraj):
        main()

        assert len(list(args_xtc_to_pymol_mdtraj.iterdir())) == 2  # annotated PDB and .pml file


#####################################
# Test Error Handling & Logging #####
#####################################

class TestTramPyMolLogging(TestArguments):

    def test_no_error_pdb(self, args_pdb_to_pymol_trace):
        main()

        assert len(list(args_pdb_to_pymol_trace.iterdir())) == 2  # annotated PDB and .pml file

    @pytest.mark.requires_MDAnalysis
    def test_no_error_xtc_mdanalysis(self, args_xtc_to_pymol_mdanalysis_trace):
        main()

        assert len(list(args_xtc_to_pymol_mdanalysis_trace.iterdir())) == 2  # annotated PDB and .pml file

    @pytest.mark.requires_MDAnalysis
    def test_error_with_trace_mdanalysis(self, args_incompatible_components_mdanalysis_trace):
        with pytest.raises(KeyError, match="Unable to find"):
            main()

    @pytest.mark.requires_mdtraj
    def test_no_error_xtc_mdtraj(self, args_xtc_to_pymol_mdtraj_trace):
        main()

        assert len(list(args_xtc_to_pymol_mdtraj_trace.iterdir())) == 2  # annotated PDB and .pml file
