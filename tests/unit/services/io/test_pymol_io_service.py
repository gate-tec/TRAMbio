import os.path
import textwrap

import pytest

from TRAMbio.services import IOServiceRegistry


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "PyMolIOService"


#################################
# No Mock Services required #####
#################################


##################
# Mock Paths #####
##################

class TestMockPaths:

    @pytest.fixture
    def out_pymol_path(self, tmp_path):
        f = str(tmp_path / "commands.pml")

        yield f

        if os.path.exists(f):
            os.remove(f)


#################
# Mock Data #####
#################

class TestMockData:

    @pytest.fixture
    def mock_data_pymol_bonds(self):
        bonds_string = textwrap.dedent("""
        bond state 2 & /DUMMY/*/A/GLY`3/O, state 2 & /DUMMY/*/A/ARG`1/H

        set_bond stick_color, yellow, /DUMMY/*/A/GLY`3/O, /DUMMY/*/A/ARG`1/H
        set_bond stick_radius, 0.1, /DUMMY/*/A/GLY`3/O, /DUMMY/*/A/ARG`1/H
        """)

        yield bonds_string


#############
# Tests #####
#############

# test write_pymol_template
class TestWritePyMolTemplate(TestMockPaths, TestMockData, TestParameters):

    def test_write_template(self, out_pymol_path):
        io_service = IOServiceRegistry.PYMOL.query_service(self.TESTED_SERVICE)

        ref_name = "DUMMY"

        io_service.write_pymol_template(
            pml_path=out_pymol_path,
            out_prefix=ref_name,
            pdb_path="dummy.pdb",
            num_states=None,
            max_color_value=50,
            bond_commands=None
        )

        # Validate output file
        assert os.path.exists(out_pymol_path)
        with open(out_pymol_path) as file:
            lines = file.readlines()

        for line in lines:
            if line.startswith("load"):
                assert "discrete=1" in line
            elif line.startswith("spectrum"):
                parts = [x.strip() for x in line.split(",")]
                assert parts[0] == "spectrum b"
                assert parts[1] == "rainbow"
                assert parts[2] == ref_name
                assert parts[3].startswith("minimum=")
                assert parts[4].startswith("maximum=")


    def test_write_template_for_trajectory(self, out_pymol_path):
        io_service = IOServiceRegistry.PYMOL.query_service(self.TESTED_SERVICE)

        ref_name = "TRAJ"

        io_service.write_pymol_template(
            pml_path=out_pymol_path,
            out_prefix=ref_name,
            pdb_path="traj.pdb",
            num_states=5,
            max_color_value=70,
            bond_commands=None
        )

        # Validate output file
        assert os.path.exists(out_pymol_path)
        with open(out_pymol_path) as file:
            lines = file.readlines()

        first_load = True
        for line in lines:
            if line.startswith("load"):
                if first_load:
                    assert "discrete=1" in line
                    first_load = False
                else:
                    assert "discrete=0" in line
            elif line.startswith("spectrum"):
                parts = [x.strip() for x in line.split(",")]
                assert parts[0] == "spectrum b"
                assert parts[1] == "rainbow"
                assert parts[2] == ref_name
                assert parts[3].startswith("minimum=")
                assert parts[4].startswith("maximum=")

    def test_write_template_with_bonds(self, out_pymol_path, mock_data_pymol_bonds):
        io_service = IOServiceRegistry.PYMOL.query_service(self.TESTED_SERVICE)

        ref_name = "DUMMY"

        io_service.write_pymol_template(
            pml_path=out_pymol_path,
            out_prefix=ref_name,
            pdb_path="dummy.pdb",
            num_states=None,
            max_color_value=10,
            bond_commands=mock_data_pymol_bonds
        )

        # Validate output file
        assert os.path.exists(out_pymol_path)
        with open(out_pymol_path) as file:
            lines = file.readlines()

        for line in lines:
            if line.startswith("load"):
                assert "discrete=1" in line
            elif line.startswith("spectrum"):
                parts = [x.strip() for x in line.split(",")]
                assert parts[0] == "spectrum b"
                assert parts[1] == "rainbow"
                assert parts[2] == ref_name
                assert parts[3].startswith("minimum=")
                assert parts[4].startswith("maximum=")
