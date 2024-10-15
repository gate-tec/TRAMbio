import time

import pytest

from TRAMbio.services import ParameterRegistry, IOServiceRegistry, StructureServiceRegistry
from TRAMbio.services.parameter import HydrogenBondParameter
from TRAMbio.util.wrapper.biopandas.pandas_pdb import CustomPandasPdb
from tests.mock.services.io import MockPdbIOService
from tests.mock.services.structure import MockHydrogenTestPdbStructureService
from tests.util.graphs import graph_methane
from tests.util.protein_graph_utils import convert_to_pdb_lines


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "mdtraj"
    ALLOWED_HASH_COLLISIONS = 1

    @pytest.fixture
    def parameters_xtc_mdanalysis_default(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrogenBondParameter.INCLUDE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_xtc_mdanalysis_no_hydrogen_test(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrogenBondParameter.INCLUDE.value, False)

        yield param_id


#################################
# Required Mock Services:
# - pdb_io_service
# - pdb_structure_service
#################################

class TestMockService:

    @pytest.fixture
    def stored_services(self):
        original_pdb_io_service = IOServiceRegistry.PDB.single_service()
        original_pdb_structure_service = StructureServiceRegistry.PDB.single_service()

        yield

        IOServiceRegistry.PDB.register_service(original_pdb_io_service)
        StructureServiceRegistry.PDB.register_service(original_pdb_structure_service)

    @pytest.fixture
    def mock_services_empty(self, stored_services):
        with pytest.warns(UserWarning, match="No ATOM/HETATM entries"):
            pdb = CustomPandasPdb().read_pdb_from_list([])

        mock_pdb_io_service = MockPdbIOService(pdb, 1336)
        IOServiceRegistry.PDB.register_service(mock_pdb_io_service)
        StructureServiceRegistry.PDB.register_service(MockHydrogenTestPdbStructureService(True))

        yield mock_pdb_io_service.get_result_counter

    @pytest.fixture
    def mock_services_methane(self, stored_services):
        coords, _, _, _, _, _ = graph_methane()
        others = [
            ["CONECT", f"{1:5d}{2:5d}{3:5d}{4:5d}{5:5d}", len(coords) + 1]
        ]

        pdb = CustomPandasPdb().read_pdb_from_list(convert_to_pdb_lines(coords, others))

        mock_pdb_io_service = MockPdbIOService(pdb, 1336)
        IOServiceRegistry.PDB.register_service(mock_pdb_io_service)
        StructureServiceRegistry.PDB.register_service(MockHydrogenTestPdbStructureService(True))

        yield mock_pdb_io_service.get_result_counter

    @pytest.fixture
    def mock_services_empty_no_hydrogen(self, stored_services):
        with pytest.warns(UserWarning, match="No ATOM/HETATM entries"):
            pdb = CustomPandasPdb().read_pdb_from_list([])

        mock_pdb_io_service = MockPdbIOService(pdb, 1336)
        IOServiceRegistry.PDB.register_service(mock_pdb_io_service)
        StructureServiceRegistry.PDB.register_service(MockHydrogenTestPdbStructureService(False))

        yield mock_pdb_io_service.get_result_counter

#############
# Tests #####
#############

class TestRead(TestMockService, TestParameters):

    # test with empty result df
    @pytest.mark.requires_mdtraj
    def test_empty_stride_50(self, mock_services_empty, path_trajectory_sample,
                             parameters_xtc_mdanalysis_no_hydrogen_test):
        xtc_io_service = IOServiceRegistry.XTC.query_service(self.TESTED_SERVICE)
        callback_function = mock_services_empty

        pdb_path, xtc_path = path_trajectory_sample

        num_frames, generator = xtc_io_service.read(
            xtc_path=xtc_path,
            pdb_path=pdb_path,
            stride=50,
            parameter_id=parameters_xtc_mdanalysis_no_hydrogen_test
        )

        total_frames = 1
        unique_frames = 1

        assert num_frames == total_frames
        for i, (frame, raw_df) in enumerate(generator):
            assert i == frame
            assert len(raw_df['ATOM']) == 0
            assert len(raw_df['HETATM']) == 0

        num_calls, num_hashes = callback_function()
        assert num_calls == total_frames
        assert unique_frames - self.ALLOWED_HASH_COLLISIONS <= num_hashes <= unique_frames

    @pytest.mark.requires_mdtraj
    def test_empty_stride_1(self, mock_services_empty, path_trajectory_sample,
                            parameters_xtc_mdanalysis_no_hydrogen_test):
        xtc_io_service = IOServiceRegistry.XTC.query_service(self.TESTED_SERVICE)
        callback_function = mock_services_empty

        pdb_path, xtc_path = path_trajectory_sample

        num_frames, generator = xtc_io_service.read(
            xtc_path=xtc_path,
            pdb_path=pdb_path,
            stride=1,
            parameter_id=parameters_xtc_mdanalysis_no_hydrogen_test
        )

        total_frames = 22
        unique_frames = 22

        assert num_frames == total_frames
        for i, (frame, raw_df) in enumerate(generator):
            assert i == frame
            assert len(raw_df['ATOM']) == 0
            assert len(raw_df['HETATM']) == 0

        num_calls, num_hashes = callback_function()
        assert num_calls == total_frames
        assert unique_frames - self.ALLOWED_HASH_COLLISIONS <= num_hashes <= unique_frames

    @pytest.mark.requires_mdtraj
    def test_methane_stride_1(self, mock_services_methane, path_trajectory_sample,
                              parameters_xtc_mdanalysis_no_hydrogen_test):
        xtc_io_service = IOServiceRegistry.XTC.query_service(self.TESTED_SERVICE)
        callback_function = mock_services_methane

        pdb_path, xtc_path = path_trajectory_sample

        num_frames, generator = xtc_io_service.read(
            xtc_path=xtc_path,
            pdb_path=pdb_path,
            stride=1,
            parameter_id=parameters_xtc_mdanalysis_no_hydrogen_test
        )

        total_frames = 22
        unique_frames = 22

        assert num_frames == total_frames
        for i, (frame, raw_df) in enumerate(generator):
            assert i == frame
            assert len(raw_df['ATOM']) == 0
            assert len(raw_df['HETATM']) == 5

        num_calls, num_hashes = callback_function()
        assert num_calls == total_frames
        assert unique_frames - self.ALLOWED_HASH_COLLISIONS <= num_hashes <= unique_frames

    @pytest.mark.requires_mdtraj
    def test_spliced_stride_1(self, mock_services_empty, path_trajectory_sample_spliced,
                              parameters_xtc_mdanalysis_no_hydrogen_test):
        xtc_io_service = IOServiceRegistry.XTC.query_service(self.TESTED_SERVICE)
        callback_function = mock_services_empty

        pdb_path, xtc_path = path_trajectory_sample_spliced

        num_frames, generator = xtc_io_service.read(
            xtc_path=xtc_path,
            pdb_path=pdb_path,
            stride=1,
            parameter_id=parameters_xtc_mdanalysis_no_hydrogen_test
        )

        total_frames = 15
        unique_frames = 2

        assert num_frames == total_frames
        for i, (frame, raw_df) in enumerate(generator):
            assert i == frame
            assert len(raw_df['ATOM']) == 0
            assert len(raw_df['HETATM']) == 0

        num_calls, num_hashes = callback_function()
        assert num_calls == total_frames
        assert unique_frames - self.ALLOWED_HASH_COLLISIONS <= num_hashes <= unique_frames

    # test for hydrogen detection
    @pytest.mark.requires_mdtraj
    def test_error_on_missing_hydrogen_atoms(self, mock_services_empty_no_hydrogen, path_trajectory_sample_spliced,
                                             parameters_xtc_mdanalysis_default):
        xtc_io_service = IOServiceRegistry.XTC.query_service(self.TESTED_SERVICE)

        pdb_path, xtc_path = path_trajectory_sample_spliced

        num_frames, generator = xtc_io_service.read(
            xtc_path=xtc_path,
            pdb_path=pdb_path,
            stride=1,
            parameter_id=parameters_xtc_mdanalysis_default
        )

        total_frames = 15

        assert num_frames == total_frames

        with pytest.raises(ValueError, match="No hydrogen atoms"):
            next(generator)