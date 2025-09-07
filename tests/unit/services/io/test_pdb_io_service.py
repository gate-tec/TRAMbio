import os.path
import textwrap
from io import StringIO

import pytest

from TRAMbio.services import IOServiceRegistry
from TRAMbio.util.wrapper.base.list import first
from TRAMbio.util.wrapper.biopandas.pandas_pdb import CustomPandasPdb
from tests.util.graphs import graph_gly_donor, graph_asn_donor, graph_phe_aromatic_xy_plane, graph_methane
from tests.util.protein_graph_utils import convert_to_pdb_lines


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "PdbIOService"


#################################
# No Mock Services required #####
#################################


##################
# Mock Paths #####
##################

class TestMockPaths:

    # output paths for pdb_file_context
    @pytest.fixture
    def out_path_pdb(self, tmp_path):
        f = str(tmp_path / "out.pdb")

        yield f

        if os.path.exists(f):
            os.remove(f)


#################
# Mock Data #####
#################

class TestMockData:

    # mock data for read
    @pytest.fixture
    def mock_data_stream(self):
        pdb_stream = StringIO()
        pdb_stream.write(textwrap.dedent("""\
            HEADER    PROTEIN                                 03-FEB-25   XXXX
            TITLE     EXAMPLE PROTEIN
            MODEL        1
            ATOM      1  CA  GLY A   1      -1.561   1.386   0.000  0.00  0.00           C
            ATOM      2  N   GLY A   1      -1.070   0.000   0.000  0.00  0.00           N
            ATOM      3  H   GLY A   1       0.000   0.000   0.000  0.00  0.00           H
            TER       4
            ATOM      5  CG  ASN A   2      -1.915  10.000   0.000  0.00  0.00           C
            ATOM      6  OD1 ASN A   2      -1.915  11.270   0.000  0.00  0.00           O
            ATOM      7  ND2 ASN A   2      -0.535  10.000   0.000  0.00  0.00           N
            ATOM      8 HD21 ASN A   2       0.000  10.927   0.000  0.00  0.00           H
            ATOM      9 HD22 ASN A   2       0.000   9.073   0.000  0.00  0.00           H
            ENDMDL
            END
            """))

        pdb_stream.seek(0)

        yield pdb_stream

        pdb_stream.close()

    @pytest.fixture
    def mock_data_stream_with_charge(self):
        """Requires module :py:mod:`TRAMbio.util.patches.biopandas`"""
        pdb_stream = StringIO()
        pdb_stream.write(textwrap.dedent("""\
        HEADER    PROTEIN                                 03-FEB-25   XXXX
        TITLE     EXAMPLE PROTEIN
        MODEL        1
        ATOM      1  NZ  LYS A   1      -0.357   0.000   0.000  0.00  0.00           N1+
        ATOM      2  HZ1 LYS A   1       0.000   1.009   0.000  0.00  0.00           H
        ATOM      3  HZ2 LYS A   1       0.000  -0.504   0.873  0.00  0.00           H
        ATOM      4  HZ3 LYS A   1       0.000  -0.504  -0.873  0.00  0.00           H
        ATOM      5  CE  LYS A   1      -1.827   0.000   0.000  0.00  0.00           C
        ENDMDL
        """))

        pdb_stream.seek(0)

        yield pdb_stream

        pdb_stream.close()

    # mock data for pdb_file_context
    @pytest.fixture
    def mock_data_pdb_header(self):
        header_stream = StringIO()
        header_stream.write(textwrap.dedent("""\
        HEADER    PROTEIN                                 03-FEB-25   XXXX
        TITLE     EXAMPLE PROTEIN
        """))

        header_stream.seek(0)

        yield header_stream

        header_stream.close()

    @pytest.fixture
    def mock_data_pdb_model_1(self):
        coords = graph_gly_donor(1, 1)[0]
        others = [["TER", f"{len(coords) + 1:5d}", len(coords) + 1]]
        coords += graph_asn_donor(len(coords) + 2, 2, y_offset=10)[0]

        model_lines = convert_to_pdb_lines(coords, others)

        yield CustomPandasPdb().read_pdb_from_list(model_lines)

    @pytest.fixture
    def mock_data_pdb_model_2(self):
        coords = graph_phe_aromatic_xy_plane(1, 1)[0]
        others = []

        model_lines = convert_to_pdb_lines(coords, others)

        yield CustomPandasPdb().read_pdb_from_list(model_lines)

    @pytest.fixture
    def mock_data_pdb_model_3(self):
        coords = graph_methane(1, 1)[0]
        others = [["CONECT", f"{1:5d}{2:5d}{3:5d}{4:5d}{5:5d}", len(coords) + 1]]

        model_lines = convert_to_pdb_lines(coords, others)

        yield CustomPandasPdb().read_pdb_from_list(model_lines)


#############
# Tests #####
#############

# test read
class TestRead(TestMockData, TestParameters):

    # - test reading from PDB file
    def test_pdb_file_empty(self, path_pdb_test_empty):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        data = io_service.read(path_pdb_test_empty)

        assert len(data.get_model(1).df["ATOM"]) == 0

    def test_pdb_file_single_model(self, path_pdb_test_single_model):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        data = io_service.read(path_pdb_test_single_model)

        assert len(data.get_model(1).df["ATOM"]) == 8

    def test_pdb_file_multiple_models(self, path_pdb_test_multiple_models):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        data = io_service.read(path_pdb_test_multiple_models)

        assert len(data.get_model(1).df["ATOM"]) == 8
        assert len(data.get_model(2).df["ATOM"]) == 11

    def test_pdb_file_hetatm_model(self, path_pdb_test_hetatm):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        data = io_service.read(path_pdb_test_hetatm)
        model = data.get_model(0)

        assert len(model.df["ATOM"]) == 0
        assert len(model.df["HETATM"]) == 5
        assert len(model.df["OTHERS"].loc[model.df["OTHERS"]["record_name"] == "CONECT", :]) == 1

    def test_pdb_file_regular_size(self, path_pdb_sample):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        data = io_service.read(path_pdb_sample)

        assert len(data.get_model(1).df["ATOM"]) == 1336

    # - test reading from stream
    def test_stream(self, mock_data_stream):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        data = io_service.read(mock_data_stream)
        model = data.get_model(1)

        assert len(model.df["ATOM"]) == 8
        assert len(model.df["HETATM"]) == 0
        assert len(model.df["OTHERS"].loc[model.df["OTHERS"]["record_name"] == "TER", :]) == 1

    def test_stream_with_charge(self, mock_data_stream_with_charge):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        data = io_service.read(mock_data_stream_with_charge)

        atom_df = data.get_model(1).df["ATOM"]

        assert len(atom_df) == 5
        assert len(atom_df.loc[atom_df.charge == "1+", :]) > 0


    # - test reading from PDB code
    @pytest.mark.internet
    def test_pdb_code(self):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        data = io_service.read("1l2y")

        assert len(data.get_model(1).df["ATOM"]) == 304

    @pytest.mark.internet
    @pytest.mark.slow
    def test_faulty_pdb_code(self):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        with pytest.raises(ValueError, match="Unknown PDB-code"):
            io_service.read("XXXX")


# pdb_file_context
class TestPdbFileContext(TestMockPaths, TestMockData, TestParameters):

    def test_empty_pdb(self, out_path_pdb, mock_data_pdb_header):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        with io_service.pdb_file_context(out_path_pdb, mock_data_pdb_header) as context:
            pass

        # verify file exists
        assert os.path.exists(out_path_pdb)
        with open(out_path_pdb) as file:
            lines = file.readlines()

        assert lines[0].startswith("HEADER")
        assert lines[1].startswith("TITLE")
        assert lines[2].startswith("END")


    def test_single_model(self, out_path_pdb, mock_data_pdb_header, mock_data_pdb_model_1):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        with io_service.pdb_file_context(out_path_pdb, mock_data_pdb_header) as context:
            context.write_model(model=mock_data_pdb_model_1, model_idx=1)

        # verify file exists
        assert os.path.exists(out_path_pdb)
        with open(out_path_pdb) as file:
            lines = file.readlines()

        assert lines[0].startswith("HEADER")
        assert lines[1].startswith("TITLE")

        assert lines[2].startswith("MODEL") and "1" in lines[2]
        assert lines[-2].startswith("ENDMDL")

        assert lines[-1].startswith("END")

    def test_multiple_models(self, out_path_pdb, mock_data_pdb_header, mock_data_pdb_model_1, mock_data_pdb_model_2):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        with io_service.pdb_file_context(out_path_pdb, mock_data_pdb_header) as context:
            context.write_model(model=mock_data_pdb_model_1, model_idx=1)
            context.write_model(model=mock_data_pdb_model_2, model_idx=2)

        # verify file exists
        assert os.path.exists(out_path_pdb)
        with open(out_path_pdb) as file:
            lines = file.readlines()

        assert lines[0].startswith("HEADER")
        assert lines[1].startswith("TITLE")

        # check first model
        assert lines[2].startswith("MODEL") and "1" in lines[2]
        end_first_model = first(idx for idx, line in enumerate(lines) if line.startswith("ENDMDL"))
        assert end_first_model is not None

        # check second model
        assert lines[end_first_model + 1].startswith("MODEL") and "2" in lines[end_first_model + 1]
        assert lines[-2].startswith("ENDMDL")

        assert lines[-1].startswith("END")

    def test_hetatm_model(self, out_path_pdb, mock_data_pdb_header, mock_data_pdb_model_3):
        io_service = IOServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        with io_service.pdb_file_context(out_path_pdb, mock_data_pdb_header) as context:
            context.write_model(model=mock_data_pdb_model_3, model_idx=0)

        # verify file exists
        assert os.path.exists(out_path_pdb)
        with open(out_path_pdb) as file:
            lines = file.readlines()

        assert lines[0].startswith("HEADER")
        assert lines[1].startswith("TITLE")

        # check model
        assert lines[2].startswith("MODEL") and "0" in lines[2]
        for line in lines[3:-3]:
            assert line.startswith("HETATM")
        assert lines[-3].startswith("CONECT")
        assert lines[-2].startswith("ENDMDL")

        assert lines[-1].startswith("END")
