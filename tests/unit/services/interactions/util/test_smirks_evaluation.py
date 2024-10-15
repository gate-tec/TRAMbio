import pytest

from TRAMbio.services.interactions.util import PotentialEvaluator
from TRAMbio.util.constants.smirnoff import LENNARD_JONES_12_6
from tests.util.graphs import construct_protein_graph_base, graph_ser_acceptor, graph_hie_imidazole_ring, \
    graph_phe_aromatic_xy_plane


###################
# Mock Graphs #####
###################

class TestMockGraphs:

    @pytest.fixture
    def mock_protein_graph_ser(self):
        coords, hm, cov, pe, pm, pb = graph_ser_acceptor(1, 1)

        yield construct_protein_graph_base(
            coords=coords,
            others=[],
            hm=hm,
            cov=cov,
            pe=pe,
            pm=pm,
            pb=pb
        )

    @pytest.fixture
    def mock_protein_graph_hie(self):
        coords, hm, cov, pe, pm, pb = graph_hie_imidazole_ring(1, 1)

        yield construct_protein_graph_base(
            coords=coords,
            others=[],
            hm=hm,
            cov=cov,
            pe=pe,
            pm=pm,
            pb=pb
        )

    @pytest.fixture
    def mock_protein_graph_phe(self):
        coords, hm, cov, pe, pm, pb = graph_phe_aromatic_xy_plane(1, 1)

        yield construct_protein_graph_base(
            coords=coords,
            others=[],
            hm=hm,
            cov=cov,
            pe=pe,
            pm=pm,
            pb=pb
        )


#############
# Tests #####
#############


class TestEvaluator(TestMockGraphs):

    @staticmethod
    def get_values_for_smirks(smirks: str):
        for entry in LENNARD_JONES_12_6:
            if entry['smirks'] == smirks:
                return entry['epsilon'], entry['rmin_half']
        raise ValueError(f"Unknown smirks {smirks}")

    def test_ser_hg(self, mock_protein_graph_ser):
        evaluator = PotentialEvaluator(
            atom_graph=mock_protein_graph_ser.graphs['full'],
            atom_df=mock_protein_graph_ser.atom_df,
            hydrogen_mapping=mock_protein_graph_ser.hydrogen_mapping
        )

        true_epsilon, true_rmin_half = self.get_values_for_smirks("[#1:1]-[#8]")
        epsilon, rmin_half = evaluator.evaluate_node("A0001-SER:HG")
        assert true_epsilon == epsilon
        assert true_rmin_half == rmin_half


    def test_hie_hd2(self, mock_protein_graph_hie):
        evaluator = PotentialEvaluator(
            atom_graph=mock_protein_graph_hie.graphs['full'],
            atom_df=mock_protein_graph_hie.atom_df,
            hydrogen_mapping=mock_protein_graph_hie.hydrogen_mapping
        )

        true_epsilon, true_rmin_half = self.get_values_for_smirks("[#1:1]-[#6X3]~[#7,#8,#9,#16,#17,#35]")
        epsilon, rmin_half = evaluator.evaluate_node("A0001-HIE:HD2")
        assert true_epsilon == epsilon
        assert true_rmin_half == rmin_half


    def test_hie_he1(self, mock_protein_graph_hie):
        evaluator = PotentialEvaluator(
            atom_graph=mock_protein_graph_hie.graphs['full'],
            atom_df=mock_protein_graph_hie.atom_df,
            hydrogen_mapping=mock_protein_graph_hie.hydrogen_mapping
        )

        true_epsilon, true_rmin_half = self.get_values_for_smirks("[#1:1]-[#6X3](~[#7,#8,#9,#16,#17,#35])~[#7,#8,#9,#16,#17,#35]")
        epsilon, rmin_half = evaluator.evaluate_node("A0001-HIE:HE1")
        assert true_epsilon == epsilon
        assert true_rmin_half == rmin_half


    def test_hie_he2(self, mock_protein_graph_hie):
        evaluator = PotentialEvaluator(
            atom_graph=mock_protein_graph_hie.graphs['full'],
            atom_df=mock_protein_graph_hie.atom_df,
            hydrogen_mapping=mock_protein_graph_hie.hydrogen_mapping
        )

        true_epsilon, true_rmin_half = self.get_values_for_smirks("[#1:1]-[#7]")
        epsilon, rmin_half = evaluator.evaluate_node("A0001-HIE:HE2")
        assert true_epsilon == epsilon
        assert true_rmin_half == rmin_half


    def test_phe_hz(self, mock_protein_graph_phe):
        evaluator = PotentialEvaluator(
            atom_graph=mock_protein_graph_phe.graphs['full'],
            atom_df=mock_protein_graph_phe.atom_df,
            hydrogen_mapping=mock_protein_graph_phe.hydrogen_mapping
        )

        true_epsilon, true_rmin_half = self.get_values_for_smirks("[#1:1]-[#6X3]")
        epsilon, rmin_half = evaluator.evaluate_node("A0001-PHE:HZ")
        assert true_epsilon == epsilon
        assert true_rmin_half == rmin_half
