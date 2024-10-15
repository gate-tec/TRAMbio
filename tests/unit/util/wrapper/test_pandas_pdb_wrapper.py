import textwrap

import pytest
from TRAMbio.util.wrapper.biopandas.pandas_pdb import CustomPandasPdb
from biopandas.pdb import PandasPdb

from TRAMbio.util.wrapper.pandas.warning import WarningWrapper


#################
# Mock Data #####
#################

class TestMockData:

    @pytest.fixture
    def example_pdb_data(self):
        raw_pdb = textwrap.dedent("""
        ATOM      1  N   TYR A   1      55.430  44.410  83.560  1.00  0.00           N
        ATOM      2  HT1 TYR A   1      55.770  43.830  82.770  1.00  0.00           H
        ATOM      3  HT2 TYR A   1      55.690  44.060  84.510  1.00  0.00           H
        ATOM      4  HT3 TYR A   1      55.990  45.280  83.460  1.00  0.00           H
        ATOM      5  CA  TYR A   1      53.970  44.790  83.420  1.00  0.00           C
        ATOM      6  HA  TYR A   1      53.930  45.430  82.550  1.00  0.00           H
        ATOM      7  CB  TYR A   1      53.060  43.470  83.280  1.00  0.00           C
        ATOM      8  HB1 TYR A   1      51.980  43.710  83.330  1.00  0.00           H
        ATOM      9  HB2 TYR A   1      53.380  42.850  84.150  1.00  0.00           H
        ATOM     10  CG  TYR A   1      53.310  42.680  82.070  1.00  0.00           C
        ATOM     11  CD1 TYR A   1      53.280  43.290  80.830  1.00  0.00           C
        ATOM     12  HD1 TYR A   1      53.110  44.350  80.780  1.00  0.00           H
        ATOM     13  CE1 TYR A   1      53.490  42.610  79.610  1.00  0.00           C
        ATOM     14  HE1 TYR A   1      53.470  43.180  78.690  1.00  0.00           H
        ATOM     15  CZ  TYR A   1      53.620  41.180  79.680  1.00  0.00           C
        ATOM     16  OH  TYR A   1      53.550  40.380  78.470  1.00  0.00           O
        ATOM     17  HH  TYR A   1      53.760  40.920  77.710  1.00  0.00           H
        ATOM     18  CD2 TYR A   1      53.420  41.270  82.140  1.00  0.00           C
        ATOM     19  HD2 TYR A   1      53.290  40.730  83.070  1.00  0.00           H
        ATOM     20  CE2 TYR A   1      53.550  40.550  80.970  1.00  0.00           C
        ATOM     21  HE2 TYR A   1      53.480  39.470  80.960  1.00  0.00           H
        ATOM     22  C   TYR A   1      53.470  45.570  84.550  1.00  0.00           C
        ATOM     23  O   TYR A   1      54.090  45.650  85.610  1.00  0.00           O
        TER      24
        """)

        yield raw_pdb


#############
# Tests #####
#############

class TestWrapper(TestMockData):

    @pytest.mark.xfail(raises=KeyError, reason="PandasPdb.to_pdb_stream cannot handle OTHERS records correctly")
    def test_convert_to_pdb_stream(self, example_pdb_data):
        pdb_df = PandasPdb().read_pdb_from_list(example_pdb_data.split('\n'))
        with WarningWrapper():
            pdb_df.to_pdb_stream(records=('ATOM', 'HETATM', 'OTHERS'))

    def test_convert_to_pdb_stream_with_custom_wrapper(self, example_pdb_data):
        pdb_df = CustomPandasPdb().read_pdb_from_list(example_pdb_data.split('\n'))
        pdb_df.to_pdb_stream(records=('ATOM', 'HETATM', 'OTHERS'))
