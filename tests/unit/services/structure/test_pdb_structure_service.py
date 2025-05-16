import time

import pandas as pd
import pytest

from TRAMbio.services import StructureServiceRegistry, ParameterRegistry, InteractionServiceRegistry
from TRAMbio.services.parameter import PdbParameter
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import GraphKey, ProteinGraph
from TRAMbio.util.wrapper.biopandas.pandas_pdb import CustomPandasPdb
from tests.mock.services.interactions.mock_aromatic_interaction_service import MockAromaticInteractionService
from tests.util.graphs import construct_protein_graph_base, graph_gly_donor, graph_gly_acceptor, graph_methane, \
    graph_phe_aromatic_xy_plane
from tests.util.protein_graph_utils import *
from tests.mock.services.interactions import *


##################
# Parameters #####
##################

# The only relevant parameters are UNIQUE_BONDS and KEEP_HETS

class TestParameters:

    TESTED_SERVICE = "PdbStructureService"
    LENGTH_EPSILON = 0.05

    @pytest.fixture
    def parameters_pdb_structure_default(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(PdbParameter.KEEP_HETS.value, True)
        registry.set_parameter(PdbParameter.UNIQUE_BONDS.value, False)

        yield param_id

    @pytest.fixture
    def parameters_pdb_structure_exclude_hets(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(PdbParameter.KEEP_HETS.value, False)
        registry.set_parameter(PdbParameter.UNIQUE_BONDS.value, False)

        yield param_id

    @pytest.fixture
    def parameters_pdb_structure_unique_bonds(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(PdbParameter.KEEP_HETS.value, True)
        registry.set_parameter(PdbParameter.UNIQUE_BONDS.value, True)

        yield param_id


#####################
# Mock Services #####
#####################

class TestMockServices:

    @pytest.fixture
    def mock_services_cov_empty(self):
        # store original registry
        original_services = InteractionServiceRegistry.COV.list_services()

        # inject mock services
        InteractionServiceRegistry.COV.register_service(MockEmptyDisulphideBridgeService())
        InteractionServiceRegistry.COV.register_service(MockEmptySSBondEntryService())
        InteractionServiceRegistry.COV.register_service(MockEmptyLinkEntryService())
        InteractionServiceRegistry.COV.register_service(MockEmptyConectEntryService())

        yield

        # restore registry
        for service in original_services:
            InteractionServiceRegistry.COV.register_service(service)

    @pytest.fixture
    def mock_services_non_cov_empty(self):
        # store original registry
        original_services = InteractionServiceRegistry.NON_COV.list_services()

        # inject mock services
        InteractionServiceRegistry.NON_COV.register_service(MockEmptyHydrogenBondService())
        InteractionServiceRegistry.NON_COV.register_service(MockEmptyHydrophobicPotentialInteractionService())
        InteractionServiceRegistry.NON_COV.register_service(MockEmptyHydrophobicSurfaceInteractionService())
        InteractionServiceRegistry.NON_COV.register_service(MockEmptyCationPiInteractionService())
        InteractionServiceRegistry.NON_COV.register_service(MockEmptyAromaticInteractionService())

        yield

        # restore registry
        for service in original_services:
            InteractionServiceRegistry.NON_COV.register_service(service)

    @pytest.fixture
    def mock_services_non_cov_hydrophobic_aromatic(self):
        # store original registry
        original_services = InteractionServiceRegistry.NON_COV.list_services()

        # inject mock services
        InteractionServiceRegistry.NON_COV.register_service(MockEmptyHydrogenBondService())
        InteractionServiceRegistry.NON_COV.register_service(MockEmptyHydrophobicPotentialInteractionService())
        InteractionServiceRegistry.NON_COV.register_service(MockEmptyCationPiInteractionService())
        # non-empty services for aromatic and hydrophobic interaction
        InteractionServiceRegistry.NON_COV.register_service(MockHydrophobicSurfaceInteractionService())
        InteractionServiceRegistry.NON_COV.register_service(MockAromaticInteractionService())

        yield

        # restore registry
        for service in original_services:
            InteractionServiceRegistry.NON_COV.register_service(service)


#################
# Mock Data #####
#################

class TestMockData:

    @pytest.fixture
    def mock_protein_data_empty(self):
        """no coordinates"""
        coords = []
        others = []

        pdb_lines = convert_to_pdb_lines(coords, others)

        with pytest.warns(UserWarning, match="No ATOM/HETATM entries"):
            yield CustomPandasPdb().read_pdb_from_list(pdb_lines).export_first_model(), load_as_atom_df(coords)

    @pytest.fixture
    def mock_protein_data_title_only(self):
        """no coordinates"""
        coords = []
        others = [
            ["TITLE", "    EMPTY DATA", 1]
        ]

        pdb_lines = convert_to_pdb_lines(coords, others)

        with pytest.warns(UserWarning, match="No ATOM/HETATM entries"):
            yield CustomPandasPdb().read_pdb_from_list(pdb_lines).export_first_model(), load_as_atom_df(coords)

    @pytest.fixture
    def mock_protein_data_gly_acceptor(self):
        """single residue ATOM"""
        coords, _, _, _, _, _ = graph_gly_acceptor()
        others = []

        pdb_lines = convert_to_pdb_lines(coords, others)

        yield CustomPandasPdb().read_pdb_from_list(pdb_lines).export_first_model(), load_as_atom_df(coords)

    @pytest.fixture
    def mock_protein_data_gly_gly(self):
        """multi-residue ATOM"""
        coords, _, _, _, _, _ = graph_gly_acceptor(y_offset=3.89)
        coords2, _, _, _, _, _ = graph_gly_donor(len(coords) + 1, 2, rotations=[(180, 'Z')])
        others = []
        coords += coords2

        pdb_lines = convert_to_pdb_lines(coords, others)

        yield CustomPandasPdb().read_pdb_from_list(pdb_lines).export_first_model(), load_as_atom_df(coords)

    @pytest.fixture
    def mock_protein_data_gly_gly_no_peptide(self):
        """multi-residue ATOM"""
        coords, _, _, _, _, _ = graph_gly_acceptor(y_offset=3.89)
        coords2, _, _, _, _, _ = graph_gly_donor(len(coords) + 2, 2, rotations=[(180, 'Z')])
        others = [
            ["TER", f"{len(coords) + 1:5d}", len(coords) + 1]
        ]
        coords += coords2

        pdb_lines = convert_to_pdb_lines(coords, others)

        yield CustomPandasPdb().read_pdb_from_list(pdb_lines).export_first_model(), load_as_atom_df(coords)

    @pytest.fixture
    def mock_protein_data_methane(self):
        """single residue HETATM"""
        coords, _, _, _, _, _ = graph_methane()
        others = [
            ["CONECT", f"{1:5d}{2:5d}{3:5d}{4:5d}{5:5d}", len(coords) + 1]
        ]

        pdb_lines = convert_to_pdb_lines(coords, others)

        yield CustomPandasPdb().read_pdb_from_list(pdb_lines).export_first_model(), load_as_atom_df(coords)

    @pytest.fixture
    def mock_protein_data_gly_gly_methane(self):
        """multi-residue ATOM and single residue HETATM"""
        coords, _, _, _, _, _ = graph_gly_acceptor(x_offset=20, y_offset=3.89)
        coords2, _, _, _, _, _ = graph_gly_donor(len(coords) + 1, 2, x_offset=20, rotations=[(180, 'Z')])
        coords += coords2
        lc = len(coords)

        coords2, _, _, _, _, _ = graph_methane(len(coords) + 2, 3)
        others = [
            ["TER", f"{lc + 1:5d}", lc + 1],
            ["CONECT", f"{lc + 2:5d}{lc + 3:5d}{lc + 4:5d}{lc + 5:5d}{lc + 6:5d}", lc + len(coords2) + 2]
        ]
        coords += coords2

        pdb_lines = convert_to_pdb_lines(coords, others)

        yield CustomPandasPdb().read_pdb_from_list(pdb_lines).export_first_model(), load_as_atom_df(coords)

    @pytest.fixture
    def mock_protein_data_methane_duplicate_atom_names(self):
        """Duplicate atom names in same residue, resulting in duplicate node IDs"""
        d = 0.658  # d = 0.658179 => sqrt(3 * d^2) == 1.14, however, PDB only allows 3 digits
        coords = [
            ["HETATM", 1, "A", 1, "", "CH4", "C", f"A{1:04d}-CH4:C", "C", "",
             0, 0, 0, 1],
            ["HETATM", 2, "A", 1, "", "CH4", "H", f"A{1:04d}-CH4:H", "H", "",
             d, -d, d, 2],
            ["HETATM", 3, "A", 1, "", "CH4", "H", f"A{1:04d}-CH4:H", "H", "",
             d, d, -d, 3],
            ["HETATM", 4, "A", 1, "", "CH4", "H", f"A{1:04d}-CH4:H", "H", "",
             -d, d, d, 4],
            ["HETATM", 5, "A", 1, "", "CH4", "H", f"A{1:04d}-CH4:H", "H", "",
             -d, -d, -d, 5],
        ]
        others = [
            ["CONECT", f"{1:5d}{2:5d}{3:5d}{4:5d}{5:5d}", len(coords) + 1]
        ]

        pdb_lines = convert_to_pdb_lines(coords, others)

        yield CustomPandasPdb().read_pdb_from_list(pdb_lines).export_first_model(), load_as_atom_df(coords)


###################
# Mock Graphs #####
###################

class TestMockGraphs:

    @pytest.fixture
    def mock_protein_graph_gly_acceptor(self):
        """single residue ATOM"""
        coords, hm, cov, pe, pm, pb = graph_gly_acceptor()
        others = []

        protein_graph = construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=hm,
            cov=cov,
            pe=pe,
            pm=pm,
            pb=pb
        )

        yield protein_graph.atom_df.copy(), protein_graph.others_df.copy(), protein_graph

    @pytest.fixture
    def mock_extra_gly_gly_peptide(self):
        yield [
            (f"A{1:04d}-GLY:C", f"A{2:04d}-GLY:N",
             {"kind": {InteractionType.PEPTIDE_BOND.value}, "bond_length": 3.89, "base": True})
        ], [
            (f"A{2:04d}-GLY:N", f"A{1:04d}-GLY:C", {"weight": 0}),
            (f"A{1:04d}-GLY:C", f"A{2:04d}-GLY:N", {"weight": 5})
        ], {
            f"A{1:04d}-GLY:C": 1
        }, [
            (f"A{1:04d}-GLY:C", f"A{2:04d}-GLY:N", 1)
        ]

    @pytest.fixture
    def mock_protein_graph_gly_gly(self, mock_extra_gly_gly_peptide):
        """multi-residue ATOM with peptide bond"""
        coords, hm, cov, pe, pm, pb = graph_gly_acceptor(x_offset=-1.27, y_offset=3.89)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_donor(len(coords) + 1, 2, x_offset=-1.07, rotations=[(180, 'Z')])
        others = []

        extra = mock_extra_gly_gly_peptide
        cov += extra[0]
        pe += extra[1]
        pm = dict(pm, **extra[2])
        pb += extra[3]

        protein_graph = construct_protein_graph_base(
            coords=coords + coords2,
            others=others,
            hm=hm + hm2,
            cov=cov + cov2,
            pe=pe + pe2,
            pm=dict(pm, **pm2),
            pb=pb + pb2
        )

        yield protein_graph.atom_df.copy(), protein_graph.others_df.copy(), protein_graph


    @pytest.fixture
    def mock_protein_graph_gly_gly_no_peptide(self):
        """multi-residue ATOM without peptide bond"""
        coords, hm, cov, pe, pm, pb = graph_gly_acceptor(x_offset=-1.27, y_offset=3.89)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_donor(len(coords) + 1, 2, x_offset=-1.07, rotations=[(180, 'Z')])
        others = [
            ["TER", f"{len(coords) + 1:5d}", len(coords) + 1]
        ]

        protein_graph = construct_protein_graph_base(
            coords=coords + coords2,
            others=others,
            hm=hm + hm2,
            cov=cov + cov2,
            pe=pe + pe2,
            pm=dict(pm, **pm2),
            pb=pb + pb2
        )

        yield protein_graph.atom_df.copy(), protein_graph.others_df.copy(), protein_graph


    @pytest.fixture
    def mock_protein_graph_methane(self):
        """single residue HETATM"""
        coords, hm, cov, pe, pm, pb = graph_methane()
        others = [
            ["CONECT", f"{1:5d}{2:5d}{3:5d}{4:5d}{5:5d}", len(coords) + 1]
        ]

        protein_graph = construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=hm,
            cov=cov,
            pe=pe,
            pm=pm,
            pb=pb
        )

        yield protein_graph.atom_df.copy(), protein_graph.others_df.copy(), protein_graph


    @pytest.fixture
    def mock_protein_graph_gly_gly_methane(self):
        """multi-residue ATOM with peptide bond"""
        coords, hm, cov, pe, pm, pb = graph_gly_acceptor(x_offset=20 - 1.27, y_offset=3.89)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_donor(len(coords) + 1, 2, x_offset=20 - 1.07, rotations=[(180, 'Z')])
        lc = len(coords) + len(coords2) + 2
        coords3, hm3, cov3, pe3, pm3, pb3 = graph_methane(lc, 3)
        others = [
            ["TER", f"{lc + 1:5d}", lc + 1],
            ["CONECT", f"{lc + 2:5d}{lc + 3:5d}{lc + 4:5d}{lc + 5:5d}{lc + 6:5d}", lc + len(coords3) + 2]
        ]

        cov += [
            (f"A{1:04d}-GLY:C", f"A{2:04d}-GLY:N",
             {"kind": {InteractionType.PEPTIDE_BOND.value}, "bond_length": 3.89, "base": True})
        ]
        pe += [
            (f"A{2:04d}-GLY:N", f"A{1:04d}-GLY:C", {"weight": 0}),
            (f"A{1:04d}-GLY:C", f"A{2:04d}-GLY:N", {"weight": 5})
        ]
        pm[f"A{1:04d}-GLY:C"] = 1
        pb += [
            (f"A{1:04d}-GLY:C", f"A{2:04d}-GLY:N", 1)
        ]

        protein_graph = construct_protein_graph_base(
            coords=coords + coords2 + coords3,
            others=others,
            hm=hm + hm2 + hm3,
            cov=cov + cov2 + cov3,
            pe=pe + pe2 + pe3,
            pm=dict(pm, **pm2, **pm3),
            pb=pb + pb2 + pb3
        )

        yield protein_graph.atom_df.copy(), protein_graph.others_df.copy(), protein_graph

    # mock data for copy_graph_for_frame
    @pytest.fixture
    def mock_protein_graph_gly_gly_alt_coords(self, mock_extra_gly_gly_peptide):
        """multi-residue ATOM with peptide bond"""
        coords, hm, cov, pe, pm, pb = graph_gly_acceptor(x_offset=-1.27 + 10, y_offset=3.89, z_offset=5)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_donor(len(coords) + 1, 2, x_offset=-1.07 + 10, z_offset=5, rotations=[(180, 'Z')])
        others = []

        extra = mock_extra_gly_gly_peptide
        cov += extra[0]
        pe += extra[1]
        pm = dict(pm, **extra[2])
        pb += extra[3]

        protein_graph = construct_protein_graph_base(
            coords=coords + coords2,
            others=others,
            hm=hm + hm2,
            cov=cov + cov2,
            pe=pe + pe2,
            pm=dict(pm, **pm2),
            pb=pb + pb2
        )

        yield protein_graph.atom_df.copy(), protein_graph.others_df.copy(), protein_graph

    @pytest.fixture
    def mock_protein_graph_methane_alt_coords(self):
        """single residue HETATM"""
        coords, hm, cov, pe, pm, pb = graph_methane(rotations=[(90, "Z")])
        others = [
            ["CONECT", f"{1:5d}{2:5d}{3:5d}{4:5d}{5:5d}", len(coords) + 1]
        ]

        protein_graph = construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=hm,
            cov=cov,
            pe=pe,
            pm=pm,
            pb=pb
        )

        yield protein_graph.atom_df.copy(), protein_graph.others_df.copy(), protein_graph

    # mock graphs for apply_non_covalent_interactions
    @pytest.fixture
    def mock_protein_graph_aromatic(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_phe_aromatic_xy_plane(1, 1, rotations=[(-90, 'Y')])  # yz-plane
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_phe_aromatic_xy_plane(len(coords1) + 2, 2, 5.0, rotations=[(-90, 'Y')])  # yz-plane
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1 + coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


##################
# Help Class #####
##################

class TestHelpFunctions(TestParameters):

    def compare_graph_data(self, protein_graph: ProteinGraph, target_protein_graph: ProteinGraph):
        for a, b in zip(sorted(protein_graph.graphs['full'].edges(data=True)),
                        sorted(target_protein_graph.graphs['full'].edges(data=True))):
            assert min(a[0], a[1]) == min(b[0], b[1])
            assert max(a[0], a[1]) == max(b[0], b[1])

            assert a[2]["kind"] == b[2]["kind"]
            assert b[2]["bond_length"] - self.LENGTH_EPSILON <= a[2]["bond_length"] <= b[2]["bond_length"] + self.LENGTH_EPSILON

        for a, b in zip(sorted(protein_graph.graphs['pebble'].nodes(data="pebbles")),
                        sorted(target_protein_graph.graphs['pebble'].nodes(data="pebbles"))):
            assert a[0] == b[0]
            assert a[1] == b[1]  # compare number of pebble

        for a, b in zip(sorted(protein_graph.graphs['pebble'].edges(data="weight")),
                        sorted(target_protein_graph.graphs['pebble'].edges(data="weight"))):
            assert a[0] == b[0]
            assert a[1] == b[1]
            assert a[2] == b[2]  # compare edge weight

        for a, b in zip(sorted(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]),
                        sorted(target_protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])):
            assert a[0] == b[0]
            assert a[1] == b[1]
            assert a[2] == b[2]  # compare edge weight

    @staticmethod
    def compare_graph_attributes(protein_graph: ProteinGraph, atom_df: pd.DataFrame):
        for node, d in protein_graph.graphs['full'].nodes(data=True):
            coords = d["coords"]
            assert coords[0] == atom_df.loc[atom_df["node_id"] == node, "x_coord"].values[0]
            assert coords[1] == atom_df.loc[atom_df["node_id"] == node, "y_coord"].values[0]
            assert coords[2] == atom_df.loc[atom_df["node_id"] == node, "z_coord"].values[0]


#############
# Tests #####
#############

# test export_atom_df
class TestExportAtomDf(TestMockData, TestParameters):

    # - base test (no hets, only unique IDs)
    def test_base(self, mock_protein_data_gly_acceptor, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, target_atom_df = mock_protein_data_gly_acceptor

        atom_df = structure_service.export_atom_df(raw_df, False, parameters_pdb_structure_default)

        difference = pd.concat([
            atom_df.loc[:, ["node_id", "x_coord", "y_coord", "z_coord"]],
            target_atom_df.loc[:, ["node_id", "x_coord", "y_coord", "z_coord"]]
        ]).reset_index(drop=True).drop_duplicates(keep=False)

        assert len(difference) == 0

    def test_empty(self, mock_protein_data_title_only, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, target_atom_df = mock_protein_data_title_only

        atom_df = structure_service.export_atom_df(raw_df, False, parameters_pdb_structure_default)

        assert len(atom_df) == 0

    # - test exclude hets
    def test_exclude_hets(self, mock_protein_data_gly_gly_methane, parameters_pdb_structure_exclude_hets):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, target_atom_df = mock_protein_data_gly_gly_methane

        atom_df = structure_service.export_atom_df(raw_df, False, parameters_pdb_structure_exclude_hets)

        difference = pd.concat([
            atom_df.loc[:, ["node_id", "x_coord", "y_coord", "z_coord"]],
            target_atom_df.loc[:, ["node_id", "x_coord", "y_coord", "z_coord"]]
        ]).reset_index(drop=True).drop_duplicates(keep=False)

        assert len(difference) == 5
        for _, row in difference.iterrows():
            assert "CH4" in row["node_id"]

    def test_exclude_hets_empty(self, mock_protein_data_methane, parameters_pdb_structure_exclude_hets):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, target_atom_df = mock_protein_data_methane

        atom_df = structure_service.export_atom_df(raw_df, False, parameters_pdb_structure_exclude_hets)

        assert len(atom_df) == 0

    # - test include hets
    def test_include_hets_only(self, mock_protein_data_methane, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, target_atom_df = mock_protein_data_methane

        atom_df = structure_service.export_atom_df(raw_df, False, parameters_pdb_structure_default)

        difference = pd.concat([
            atom_df.loc[:, ["node_id", "x_coord", "y_coord", "z_coord"]],
            target_atom_df.loc[:, ["node_id", "x_coord", "y_coord", "z_coord"]]
        ]).reset_index(drop=True).drop_duplicates(keep=False)

        assert len(difference) == 0

    def test_include_hets(self, mock_protein_data_gly_gly_methane, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, target_atom_df = mock_protein_data_gly_gly_methane

        atom_df = structure_service.export_atom_df(raw_df, False, parameters_pdb_structure_default)

        difference = pd.concat([
            atom_df.loc[:, ["node_id", "x_coord", "y_coord", "z_coord"]],
            target_atom_df.loc[:, ["node_id", "x_coord", "y_coord", "z_coord"]]
        ]).reset_index(drop=True).drop_duplicates(keep=False)

        assert len(difference) == 0

    # - test detection of multiple node ids
    def test_no_duplicate_node_ids(self, mock_protein_data_gly_gly_methane, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, target_atom_df = mock_protein_data_gly_gly_methane

        atom_df = structure_service.export_atom_df(raw_df, True, parameters_pdb_structure_default)

        difference = pd.concat([
            atom_df.loc[:, ["node_id", "x_coord", "y_coord", "z_coord"]],
            target_atom_df.loc[:, ["node_id", "x_coord", "y_coord", "z_coord"]]
        ]).reset_index(drop=True).drop_duplicates(keep=False)

        assert len(difference) == 0

    def test_detect_duplicate_node_ids(self, mock_protein_data_methane_duplicate_atom_names, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, target_atom_df = mock_protein_data_methane_duplicate_atom_names

        with pytest.raises(KeyError, match="Duplicate"):
            structure_service.export_atom_df(raw_df, True, parameters_pdb_structure_default)


# test has_hydrogen_atoms
class TestHasHydrogenAtoms(TestMockData, TestParameters):

    # - pandas pdb with hydrogen in ATOM records
    def test_detect_hydrogen(self, mock_protein_data_gly_gly, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, _ = mock_protein_data_gly_gly

        assert structure_service.has_hydrogen_atoms(raw_df)

    # - pandas pdb with no hydrogen in ATOM records
    def test_detect_no_hydrogen(self, mock_protein_data_gly_acceptor, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, _ = mock_protein_data_gly_acceptor

        assert not structure_service.has_hydrogen_atoms(raw_df)

    # - dataframe with hydrogen atoms
    def test_detect_hydrogen_in_dataframe(self, mock_protein_data_gly_gly, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        _, atom_df = mock_protein_data_gly_gly

        assert structure_service.has_hydrogen_atoms(atom_df)

    # - dataframe with no hydrogen atoms
    def test_detect_no_hydrogen_in_dataframe(self, mock_protein_data_gly_acceptor, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        _, atom_df = mock_protein_data_gly_acceptor

        assert not structure_service.has_hydrogen_atoms(atom_df)

    # - dataframe with hydrogen atoms from HETATM records
    def test_detect_hydrogen_from_hetatm_in_dataframe(self, mock_protein_data_methane, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        _, atom_df = mock_protein_data_methane

        assert structure_service.has_hydrogen_atoms(atom_df)


# test export_others_df
class TestExportOthersDf(TestMockData, TestParameters):

    # - only test with TER filter
    def test_ter_only_empty(self, mock_protein_data_title_only, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, _ = mock_protein_data_title_only

        others_df = structure_service.export_others_df(raw_df, True, parameters_pdb_structure_default)

        assert len(others_df) == 0

    def test_ter_only(self, mock_protein_data_gly_gly_no_peptide, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, _ = mock_protein_data_gly_gly_no_peptide

        others_df = structure_service.export_others_df(raw_df, True, parameters_pdb_structure_default)

        assert len(others_df) == 1


class TestExportHeaderStream(TestMockData, TestParameters):

    # test export_header_stream
    # - validate export of simple header (and exclusion of "unnecessary" lines)
    def test_title_only(self, mock_protein_data_title_only, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, _ = mock_protein_data_title_only

        header_stream = structure_service.export_header_stream(raw_df, parameter_id=parameters_pdb_structure_default)

        header_stream.seek(0)
        header_lines = header_stream.readlines()

        assert len(header_lines) == 1
        assert header_lines[0].startswith("TITLE")
        assert "EMPTY DATA" in header_lines[0]

    # - test on "empty" header (both (a) fully empty and (b) with only MODEL/END or TER data)
    def test_empty(self, mock_protein_data_empty, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, _ = mock_protein_data_empty

        header_stream = structure_service.export_header_stream(raw_df, parameter_id=parameters_pdb_structure_default)

        header_stream.seek(0)
        header_lines = header_stream.readlines()

        assert len(header_lines) == 1
        assert header_lines[0].startswith("TITLE")
        assert "PROTEIN" in header_lines[0]

    def test_empty_with_custom_name(self, mock_protein_data_empty, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, _ = mock_protein_data_empty

        protein_name = "custom protein name"

        header_stream = structure_service.export_header_stream(raw_df, protein_name, parameters_pdb_structure_default)

        header_stream.seek(0)
        header_lines = header_stream.readlines()

        assert len(header_lines) == 1
        assert header_lines[0].startswith("TITLE")
        assert protein_name.upper() in header_lines[0]

    def test_ter_only(self, mock_protein_data_gly_gly_no_peptide, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        raw_df, _ = mock_protein_data_gly_gly_no_peptide

        header_stream = structure_service.export_header_stream(raw_df, parameter_id=parameters_pdb_structure_default)

        header_stream.seek(0)
        header_lines = header_stream.readlines()

        assert len(header_lines) == 1
        assert header_lines[0].startswith("TITLE")
        assert "PROTEIN" in header_lines[0]


# test create_graph_struct (requires MockServices)
class TestCreateGraphStruct(TestMockServices, TestMockGraphs, TestHelpFunctions):

    # - test with all MockServices empty & check general data in resulting protein graph
    def test_base(self, mock_protein_graph_gly_acceptor, mock_services_cov_empty, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        atom_df, others_df, target_protein_graph = mock_protein_graph_gly_acceptor

        protein_graph = structure_service.create_graph_struct(atom_df, others_df, parameters_pdb_structure_default)

        self.compare_graph_data(
            protein_graph=protein_graph,
            target_protein_graph=target_protein_graph
        )

    def test_peptide(self, mock_protein_graph_gly_gly, mock_services_cov_empty, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        atom_df, others_df, target_protein_graph = mock_protein_graph_gly_gly

        protein_graph = structure_service.create_graph_struct(atom_df, others_df, parameters_pdb_structure_default)

        self.compare_graph_data(
            protein_graph=protein_graph,
            target_protein_graph=target_protein_graph
        )

    def test_no_peptide(self, mock_protein_graph_gly_gly_no_peptide, mock_services_cov_empty, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        atom_df, others_df, target_protein_graph = mock_protein_graph_gly_gly_no_peptide

        protein_graph = structure_service.create_graph_struct(atom_df, others_df, parameters_pdb_structure_default)

        self.compare_graph_data(
            protein_graph=protein_graph,
            target_protein_graph=target_protein_graph
        )

    def test_hetatm(self, mock_protein_graph_methane, mock_services_cov_empty, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        atom_df, others_df, target_protein_graph = mock_protein_graph_methane

        protein_graph = structure_service.create_graph_struct(atom_df, others_df, parameters_pdb_structure_default)

        self.compare_graph_data(
            protein_graph=protein_graph,
            target_protein_graph=target_protein_graph
        )

    def test_peptide_and_hetatm(self, mock_protein_graph_gly_gly_methane, mock_services_cov_empty, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        atom_df, others_df, target_protein_graph = mock_protein_graph_gly_gly_methane

        protein_graph = structure_service.create_graph_struct(atom_df, others_df, parameters_pdb_structure_default)

        self.compare_graph_data(
            protein_graph=protein_graph,
            target_protein_graph=target_protein_graph
        )


# test copy_graph_for_frame
class TestCopyGraphForFrame(TestMockServices, TestMockGraphs, TestHelpFunctions):

    # - test with mock protein graph
    def test_peptide(self, mock_protein_graph_gly_gly, mock_protein_graph_gly_gly_alt_coords, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        _, _, protein_graph = mock_protein_graph_gly_gly

        new_atom_df, new_others_df, target_protein_graph = mock_protein_graph_gly_gly_alt_coords

        copy_graph = structure_service.copy_graph_for_frame(
            atom_df=new_atom_df,
            others_df=new_others_df,
            protein_graph=protein_graph,
            parameter_id=parameters_pdb_structure_default
        )

        # test graph matching
        self.compare_graph_data(
            protein_graph=copy_graph,
            target_protein_graph=target_protein_graph
        )

        # test updated node attributes
        self.compare_graph_attributes(
            protein_graph=copy_graph,
            atom_df=new_atom_df
        )

    def test_methane(self, mock_protein_graph_methane, mock_protein_graph_methane_alt_coords, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        _, _, protein_graph = mock_protein_graph_methane

        new_atom_df, new_others_df, target_protein_graph = mock_protein_graph_methane_alt_coords

        copy_graph = structure_service.copy_graph_for_frame(
            atom_df=new_atom_df,
            others_df=new_others_df,
            protein_graph=protein_graph,
            parameter_id=parameters_pdb_structure_default
        )

        # test graph matching
        self.compare_graph_data(
            protein_graph=copy_graph,
            target_protein_graph=target_protein_graph
        )

        # test updated node attributes
        self.compare_graph_attributes(
            protein_graph=copy_graph,
            atom_df=new_atom_df
        )

    # - test with data created from create_graph_struct (empty MockServices)
    def test_with_struct_peptide(self, mock_services_cov_empty, mock_protein_graph_gly_gly, mock_protein_graph_gly_gly_alt_coords, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        atom_df, others_df, _ = mock_protein_graph_gly_gly

        new_atom_df, new_others_df, target_protein_graph = mock_protein_graph_gly_gly_alt_coords

        # construct protein_graph from service
        protein_graph = structure_service.create_graph_struct(atom_df, others_df, parameters_pdb_structure_default)

        copy_graph = structure_service.copy_graph_for_frame(
            atom_df=new_atom_df,
            others_df=new_others_df,
            protein_graph=protein_graph,
            parameter_id=parameters_pdb_structure_default
        )

        # test graph matching
        self.compare_graph_data(
            protein_graph=copy_graph,
            target_protein_graph=target_protein_graph
        )

        # test updated node attributes
        self.compare_graph_attributes(
            protein_graph=copy_graph,
            atom_df=new_atom_df
        )

    def test_with_struct_methane(self, mock_services_cov_empty, mock_protein_graph_methane, mock_protein_graph_methane_alt_coords, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)
        atom_df, others_df, _ = mock_protein_graph_methane

        new_atom_df, new_others_df, target_protein_graph = mock_protein_graph_methane_alt_coords

        # construct protein_graph from service
        protein_graph = structure_service.create_graph_struct(atom_df, others_df, parameters_pdb_structure_default)

        copy_graph = structure_service.copy_graph_for_frame(
            atom_df=new_atom_df,
            others_df=new_others_df,
            protein_graph=protein_graph,
            parameter_id=parameters_pdb_structure_default
        )

        # test graph matching
        self.compare_graph_data(
            protein_graph=copy_graph,
            target_protein_graph=target_protein_graph
        )

        # test updated node attributes
        self.compare_graph_attributes(
            protein_graph=copy_graph,
            atom_df=new_atom_df
        )


# test apply_non_covalent_interactions (requires MockServices)
class TestApplyNonCovalentInteractions(TestMockServices, TestMockGraphs, TestParameters):

    # - test with all MockServices empty
    def test_base(self, mock_services_non_cov_empty, mock_protein_graph_aromatic, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        protein_graph = mock_protein_graph_aromatic

        old_len = len(protein_graph.graphs["full"].edges)

        structure_service.apply_non_covalent_interactions(protein_graph, parameters_pdb_structure_default)

        assert old_len == len(protein_graph.graphs["full"].edges)

    # - test with multiple bond types on same edge (MockServices)
    def test_multiple_bonds(self, mock_services_non_cov_hydrophobic_aromatic, mock_protein_graph_aromatic, parameters_pdb_structure_default):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        protein_graph = mock_protein_graph_aromatic

        structure_service.apply_non_covalent_interactions(protein_graph, parameters_pdb_structure_default)

        # expect three bonds with double types
        assert len([(n1, n2) for n1, n2, kind in protein_graph.graphs['full'].edges(data='kind') if len(kind) > 1]) == 3


    def test_unique_bonds(self, mock_services_non_cov_hydrophobic_aromatic, mock_protein_graph_aromatic, parameters_pdb_structure_unique_bonds):
        structure_service = StructureServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        protein_graph = mock_protein_graph_aromatic

        structure_service.apply_non_covalent_interactions(protein_graph, parameters_pdb_structure_unique_bonds)

        # only unique bonds
        assert len([(n1, n2) for n1, n2, kind in protein_graph.graphs['full'].edges(data='kind') if len(kind) > 1]) == 0
