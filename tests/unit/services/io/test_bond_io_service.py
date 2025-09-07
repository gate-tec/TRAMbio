import os

import pytest
import pandas as pd
import numpy as np

from TRAMbio.services import IOServiceRegistry


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "BondIOService"


#################################
# No Mock Services required #####
#################################


##################
# Mock Paths #####
##################

class TestMockPaths:

    @pytest.fixture
    def out_path_bonds(self, tmp_path):
        f = str(tmp_path / "out.bnd")

        yield f

        if os.path.exists(f):
            os.remove(f)


#################
# Mock Data #####
#################

class TestMockData:

    @pytest.fixture
    def mock_data_bond_frame(self):
        bond_frame = pd.DataFrame([
            ("A0001-GLY:N", "A0001-GLY:CA", "covalent", ["covalent"], np.nan, []),
            ("A0001-GLY:CA", "A0001-GLY:C", "covalent", ["covalent"], np.nan, []),
            ("A0001-GLY:C", "A0001-GLY:O", "covalent", ["covalent"], np.nan, []),
            ("A0001-GLY:O", "A0004-ASN:H", "hbond", ["hbond"], -1.7, [130.0, 109.0, 150.0]),
        ], columns=["node1", "node2", "type", "type_set", "key", "extra"])

        yield bond_frame


#############
# Tests #####
#############

# test read
class TestRead(TestParameters):

    def test_protein_bonds(self, path_bnd_protein):
        io_service = IOServiceRegistry.BND.query_service(self.TESTED_SERVICE)

        data = io_service.read(path_bnd_protein)

        for _, row in data.iterrows():
            true_row = ~row.isnull()
            assert true_row["node1"] and true_row["node2"]
            assert true_row["type"] and true_row["type_set"]
            if row["type"] == "hbond":
                assert true_row["key"]

    def test_trajectory_bonds(self, path_bnd_trajectory):
        io_service = IOServiceRegistry.BND.query_service(self.TESTED_SERVICE)

        data = io_service.read(path_bnd_trajectory)


        for _, row in data.iterrows():
            true_row = ~row.isnull()
            assert true_row["node1"] and true_row["node2"]
            assert true_row["type"] and true_row["_merge"] and true_row["frame"]
            if row["_merge"] == "right_only":
                assert true_row["type_set"] and true_row["extra"]
            elif row["_merge"] == "left_only":
                assert not true_row["type_set"] and not true_row["key"] and not true_row["extra"]
            else:
                pytest.fail(f"Unexpected _merge identifier {row['_merge']}.")


# test store_bonds
class TestStoreBonds(TestMockPaths, TestMockData, TestParameters):

    def test_write(self, out_path_bonds, mock_data_bond_frame):
        io_service = IOServiceRegistry.BND.query_service(self.TESTED_SERVICE)

        io_service.store_bonds(out_path_bonds, mock_data_bond_frame, mode='w')

        assert os.path.exists(out_path_bonds)
        with open(out_path_bonds) as file:
            # Windows file format reads double line-breaks
            lines = [x for x in file.readlines() if bool(x.strip())]

        assert len(lines) == len(mock_data_bond_frame) + 1
        assert lines[0].startswith("node1")  # header


    def test_append(self, out_path_bonds, mock_data_bond_frame):
        io_service = IOServiceRegistry.BND.query_service(self.TESTED_SERVICE)

        io_service.store_bonds(out_path_bonds, mock_data_bond_frame, mode='w')

        assert os.path.exists(out_path_bonds)
        with open(out_path_bonds) as file:
            # Windows file format reads double line-breaks
            assert len([x for x in file.readlines() if bool(x.strip())]) == len(mock_data_bond_frame) + 1

        io_service.store_bonds(out_path_bonds, mock_data_bond_frame, mode='a')

        with open(out_path_bonds) as file:
            # Windows file format reads double line-breaks
            assert len([x for x in file.readlines() if bool(x.strip())]) == 2 * len(mock_data_bond_frame) + 1


# test get_bonds_for_key
class TestGetBondsForKey(TestParameters):

    def test_bonds_for_protein(self, path_bnd_protein):
        io_service = IOServiceRegistry.BND.query_service(self.TESTED_SERVICE)

        generator = io_service.get_bonds_for_key(path_bnd_protein, all_weighted_bonds=True)

        assert len(generator.send("-INF")) == 0

        assert len(generator.send("-6")) == 0

        set_1 = generator.send("-4").copy()
        assert len(set_1) == 1

        set_2 = generator.send("-0.1")
        assert len(set_2) == 2
        assert len(set_1.difference(set_2)) == 0


    def test_bonds_for_trajectory(self, path_bnd_trajectory):
        io_service = IOServiceRegistry.BND.query_service(self.TESTED_SERVICE)

        generator = io_service.get_bonds_for_key(path_bnd_trajectory, all_weighted_bonds=True)

        set_1 = generator.send("0").copy()
        assert len(set_1) == 1

        set_2 = generator.send("10")
        assert len(set_2) == 1
        assert len(set_1.difference(set_2)) == 1
