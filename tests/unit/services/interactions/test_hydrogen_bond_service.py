import numpy as np
import pytest

import time
import re
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import GraphKey
from TRAMbio.services import InteractionServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import HydrogenBondParameter, GeneralWorkflowParameter
from tests.conftest import inject_pytest_logger
from tests.util.graphs import graph_arg_donor, graph_glu_acceptor, graph_gly_donor, graph_gly_acceptor, \
    graph_asn_donor, graph_ser_acceptor, graph_lys_donor, graph_cyh_acceptor
from tests.util.graphs_dna import graph_adenine as graph_dna_a, graph_thymine as graph_dna_t, \
    graph_guanine as graph_dna_g, graph_cytosine as graph_dna_c
from tests.util.graphs_rna import graph_adenine as graph_rna_a, graph_uracil as graph_rna_u, \
    graph_guanine as graph_rna_g, graph_cytosine as graph_rna_c
from tests.util.protein_graph_utils import construct_protein_graph_base


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "HydrogenBondService"
    ENERGY_EPSILON = 0.05  # allowed energy deviation in kcal/mol due to rounding errors in mock data

    @pytest.fixture
    def parameters_no_hbond(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrogenBondParameter.INCLUDE.value, False)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id


    @pytest.fixture
    def parameters_hbond_default(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrogenBondParameter.INCLUDE.value, True)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id


    @pytest.fixture
    def parameters_hbond_unbounded(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrogenBondParameter.INCLUDE.value, True)
        registry.set_parameter(HydrogenBondParameter.MINIMAL_LENGTH.value, 2.0)
        registry.set_parameter(HydrogenBondParameter.ENERGY_THRESHOLD.value, 5000.0)
        registry.set_parameter(HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value, 5000.0)
        registry.set_parameter(HydrogenBondParameter.CUTOFF_DISTANCE.value, 5000.0)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id


    @pytest.fixture
    def parameters_hbond_bars_3(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrogenBondParameter.INCLUDE.value, True)
        registry.set_parameter(HydrogenBondParameter.BAR_COUNT.value, 3)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id


    @pytest.fixture
    def parameters_hbond_scaled_bars_minus5_minus0(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrogenBondParameter.INCLUDE.value, True)
        registry.set_parameter(HydrogenBondParameter.ENERGY_THRESHOLD.value, -0.0)
        registry.set_parameter(HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value, -5.0)
        registry.set_parameter(HydrogenBondParameter.BAR_COUNT.value, 5)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id


    @pytest.fixture
    def parameters_hbond_scaled_bars_minus2_minus0(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrogenBondParameter.INCLUDE.value, True)
        registry.set_parameter(HydrogenBondParameter.ENERGY_THRESHOLD.value, -0.0)
        registry.set_parameter(HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value, -2.0)
        registry.set_parameter(HydrogenBondParameter.BAR_COUNT.value, 5)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id


#################################
# No Mock Services required #####
#################################


###################
# Mock Graphs #####
###################

class TestMockGraphs:

    @pytest.fixture
    def mock_protein_graph_gly_gly(self):
        """sp2 donor and sp2 acceptor"""
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_gly_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_acceptor(len(coords1) + 2, 2, 1.2, 1.5)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1+coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


    @pytest.fixture
    def mock_protein_graph_gly_gly_large_distance(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_gly_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_acceptor(len(coords1) + 2, 2, 5.0, 1.5)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1+coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


    @pytest.fixture
    def mock_protein_graph_gly_gly_small_distance(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_gly_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_acceptor(len(coords1) + 2, 2, 2.2 - 1.07)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1+coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


    @pytest.fixture
    def mock_protein_graph_gly_gly_no_donor_base(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_gly_donor(1, 1, exclude=["CA"])
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_acceptor(len(coords1) + 2, 2, 1.2, 1.5)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1+coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


    @pytest.fixture
    def mock_protein_graph_gly_gly_no_acceptor_base(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_gly_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_acceptor(len(coords1) + 2, 2, 1.2, 1.5, exclude=["CA"])
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1+coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


    @pytest.fixture
    def mock_protein_graph_asn_glu(self):
        """sp2 donor and sp2 acceptor. double hydrogen bond (bifurcated). (Positioning -> Parameters)"""
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_asn_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_glu_acceptor(len(coords1) + 2, 2, 2.0)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1+coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


    @pytest.fixture
    def mock_protein_graph_arg_glu(self):
        """double salt-bridge"""
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_arg_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_glu_acceptor(len(coords1) + 2, 2, 1.6)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1+coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


    @pytest.fixture
    def mock_protein_graph_gly_ser(self):
        """sp2 donor and sp3 acceptor"""
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_gly_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_ser_acceptor(len(coords1) + 2, 2, 1.6)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1+coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


    @pytest.fixture
    def mock_protein_graph_lys_gly(self):
        """sp3 donor and sp2 acceptor"""
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_lys_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_acceptor(len(coords1) + 2, 2, 2.0, 1.7)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1+coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


    @pytest.fixture
    def mock_protein_graph_lys_ser(self):
        """sp3 donor and sp3 acceptor"""
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_lys_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_ser_acceptor(len(coords1) + 2, 2, 2.0, 1.7)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1+coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


    @pytest.fixture
    def mock_protein_graph_lys_glu(self):
        """Salt-bridge"""
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_lys_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_glu_acceptor(len(coords1) + 2, 2, 2.3 + 1.17, rotations=[(-90, 'Z')])
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        yield construct_protein_graph_base(
            coords=coords1+coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )


    @pytest.fixture
    def mock_protein_graph_lys_cyh(self):
        """sp3 donor and sp3 acceptor. sulfur acceptor"""
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_lys_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_cyh_acceptor(len(coords1) + 2, 2, 2.5, 2.0, 1.0)
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


    @pytest.fixture
    def mock_protein_graph_gly_gly_and_arg_glu(self):
        """single hbond and double salt-bridge"""
        coords, hm, cov, pe, pm, pb = graph_gly_donor(1, 1)
        len_coords = len(coords)
        others = [
            ["TER", f"{len_coords + 1:5d}", len_coords + 1]
        ]

        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_acceptor(len_coords + 2, 2, 1.2, 1.5)
        coords, hm, cov, pe, pm, pb = coords + coords2, hm + hm2, cov + cov2, pe + pe2, dict(pm, **pm2), pb + pb2
        len_coords += len(coords2) + 2
        others += [
            ["TER", f"{len_coords + 1:5d}", len_coords + 1]
        ]

        coords2, hm2, cov2, pe2, pm2, pb2 = graph_arg_donor(len_coords + 2, 3, 0.0, 6.0)
        coords, hm, cov, pe, pm, pb = coords + coords2, hm + hm2, cov + cov2, pe + pe2, dict(pm, **pm2), pb + pb2
        len_coords += len(coords2) + 2
        others += [
            ["TER", f"{len_coords + 1:5d}", len_coords + 1]
        ]

        coords2, hm2, cov2, pe2, pm2, pb2 = graph_glu_acceptor(len_coords + 2, 4, 1.6, 6.6)
        coords, hm, cov, pe, pm, pb = coords + coords2, hm + hm2, cov + cov2, pe + pe2, dict(pm, **pm2), pb + pb2

        yield construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=hm,
            cov=cov,
            pe=pe,
            pm=pm,
            pb=pb
        )

    # DNA & RNA

    @pytest.fixture
    def mock_protein_graph_dna_a_t(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_dna_a(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_dna_t(len(coords1) + 2, 2, -0.5, 2.9, -0.25, rotations=[(180, 'Z'), (-175, 'Y')])
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

    @pytest.fixture
    def mock_protein_graph_dna_g_c(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_dna_g(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_dna_c(len(coords1) + 2, 2, 0, 2.9, 0,
                                                        rotations=[(180, 'Z'), (180, 'Y')])
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

    @pytest.fixture
    def mock_protein_graph_rna_a_u(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_rna_a(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_rna_u(len(coords1) + 2, 2, -0.5, 2.9, -0.25,
                                                        rotations=[(180, 'Z'), (-175, 'Y')])
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

    @pytest.fixture
    def mock_protein_graph_rna_g_c(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_rna_g(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_rna_c(len(coords1) + 2, 2, 0, 2.9, 0,
                                                        rotations=[(180, 'Z'), (180, 'Y')])
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


#############
# Tests #####
#############

class TestApplyInteractions(TestMockGraphs, TestParameters):

    # correctly skip detection by parameter HydrogenBondParameter.INCLUDE.value == False
    def test_skip_eval(self, mock_protein_graph_gly_gly, parameters_no_hbond):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly

        interaction_service.apply_interactions(protein_graph, parameters_no_hbond, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 0


    # detect no h-bonds
    def test_detect_no_hbond(self, mock_protein_graph_gly_gly_large_distance, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly_large_distance

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 0


    # detect single h-bond (no energy cut-off, detect correct insertion of H-atom in pebble graph)
    def test_detect_single_hbond(self, mock_protein_graph_gly_gly, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value][0]
        assert "A0001-GLY:H" in bond
        assert "A0002-GLY:O" in bond
        assert bond[2] == 5  # full strength hydrogen bond
        assert bond[3] < 0  # negative energy value

        assert InteractionType.H_BOND.value in protein_graph.graphs['full'].edges['A0001-GLY:H', 'A0002-GLY:O']['kind']


    # detect positive energy value for too small distance (requires unbounded parameters)
    def test_detect_energy_positive_hbond(self, mock_protein_graph_gly_gly_small_distance, parameters_hbond_unbounded):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly_small_distance

        interaction_service.apply_interactions(protein_graph, parameters_hbond_unbounded, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value][0]
        assert "A0001-GLY:H" in bond
        assert "A0002-GLY:O" in bond
        assert bond[2] == 5  # full strength hydrogen bond
        assert bond[3] > 0  # positive energy value

        assert InteractionType.H_BOND.value in protein_graph.graphs['full'].edges['A0001-GLY:H', 'A0002-GLY:O']['kind']


    # detect multiple h-bonds with correct energy values and angles
    def test_hbond_sp2_sp2(self, mock_protein_graph_gly_gly, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value][0]
        assert "A0001-GLY:H" in bond
        assert "A0002-GLY:O" in bond
        assert bond[2] == 5  # full strength hydrogen bond
        assert bond[3] < 0  # negative energy value

        # sp2-sp2 only has angles: theta, phi, gamma
        assert len(protein_graph.graphs['full'].edges['A0001-GLY:H', 'A0002-GLY:O']['extra']) == 3

        # d vector is (2.27, 1.5, 0), |d| is ~2.72
        d = 2.72
        # r vector is (1.2, 1.5, 0), |r| is ~1.92

        # N-H vector is (1.07, 0, 0) | N-H unit vector is (1, 0, 0)
        # O-H vector is -r           | H-O unit vector is (-0.625, -0.781, 0)
        # => theta is arccos(-0.625) ~ 128.68 degrees (2.246 rad)
        theta = 2.246

        # N-H and O-base vectors are parallel
        # => phi = theta
        phi = theta

        # both sp2-centers are in xy-plane
        # => gamma = 180 degrees (3.14 rad)
        gamma = 3.14
        # max(phi, gamma) = gamma

        angular_term = np.cos(theta)**2 * np.exp(-(np.pi - theta)**6) * np.cos(gamma)**2

        # energy = 8 * [5 * (2.80 / d) ** 12 - 6 * (2.80 / d) ** 10] * angular_term
        energy = 8 * (5 * (2.80 / d) ** 12 - 6 * (2.80 / d) ** 10) * angular_term

        assert energy - self.ENERGY_EPSILON <= bond[3] <= energy + self.ENERGY_EPSILON



    def test_hbond_sp2_sp3(self, mock_protein_graph_gly_ser, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_ser

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value][0]
        assert "A0001-GLY:H" in bond
        assert "A0002-SER:OG" in bond
        assert bond[2] == 5  # full strength hydrogen bond
        assert bond[3] < 0  # negative energy value

        # sp2-sp3 only has theta angle
        assert len(protein_graph.graphs['full'].edges['A0001-GLY:H', 'A0002-SER:OG']['extra']) == 1

        # d vector is (2.67, 0, 0), |d| is 2.67
        d = 2.67
        # r vector is (1.6, 0, 0), |r| is 1.6

        # N-H vector is (1.07, 0, 0) | N-H unit vector is (1, 0, 0)
        # O-H vector is -r           | H-O unit vector is (-1, 0, 0)
        # => theta is arccos(-1) = 180 degrees (PI rad)
        theta = np.pi

        angular_term = np.cos(theta)**4 * np.exp(-2 * (np.pi - theta)**6)

        # energy = 8 * [5 * (2.80 / d) ** 12 - 6 * (2.80 / d) ** 10] * angular_term
        energy = 8 * (5 * (2.80 / d) ** 12 - 6 * (2.80 / d) ** 10) * angular_term

        assert energy - self.ENERGY_EPSILON <= bond[3] <= energy + self.ENERGY_EPSILON


    def test_hbond_sp3_sp2(self, mock_protein_graph_lys_gly, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_gly

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value][0]
        assert "A0001-LYS:HZ1" in bond
        assert "A0002-GLY:O" in bond
        assert bond[2] == 5  # full strength hydrogen bond
        assert bond[3] < 0  # negative energy value

        # sp3-sp2 has angles: theta, phi
        assert len(protein_graph.graphs['full'].edges['A0001-LYS:HZ1', 'A0002-GLY:O']['extra']) == 2

        # d vector is (2.357, 1.69, 0), |d| is ~2.9
        d = 2.9
        # r vector is (2, 0.69, 0), |r| is ~2.12

        # N-H vector is ~ (0.357, 1.00, 0) | N-H unit vector is (0.336, 0.941, 0)
        # O-H vector is -r                 | H-O unit vector is (-0.945, -0.326, 0)
        # => theta is arccos(-0.624) ~ 128.6 degrees (2.245 rad)
        theta = 2.245

        # O-base vector is (1.27, 0, 0)    | O-base unit vector is (1, 0, 0)
        # => phi = arccos(-0.945) ~ 160.9 degrees (2.81 rad)
        phi = 2.81

        angular_term = np.cos(theta)**2 * np.exp(-(np.pi - theta)**6) * np.cos(phi)**2

        # energy = 8 * [5 * (2.80 / d) ** 12 - 6 * (2.80 / d) ** 10] * angular_term
        energy = 8 * (5 * (2.80 / d) ** 12 - 6 * (2.80 / d) ** 10) * angular_term

        assert energy - self.ENERGY_EPSILON <= bond[3] <= energy + self.ENERGY_EPSILON


    def test_hbond_sp3_sp3(self, mock_protein_graph_lys_ser, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_ser

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value][0]
        assert "A0001-LYS:HZ1" in bond
        assert "A0002-SER:OG" in bond
        assert bond[2] == 5  # full strength hydrogen bond
        assert bond[3] < 0  # negative energy value

        # sp3-sp3 has angles: theta, phi
        assert len(protein_graph.graphs['full'].edges['A0001-LYS:HZ1', 'A0002-SER:OG']['extra']) == 2

        # d vector is (2.357, 1.69, 0), |d| is ~2.9
        d = 2.9
        # r vector is (2, 0.69, 0), |r| is ~2.12

        # N-H vector is ~ (0.357, 1.00, 0) | N-H unit vector is (0.336, 0.941, 0)
        # O-H vector is -r                 | H-O unit vector is (-0.945, -0.326, 0)
        # => theta is arccos(-0.624) ~ 128.6 degrees (2.245 rad)
        theta = 2.245

        # O-base vector is (1.394, -0.361, 0)    | O-base unit vector is (0.968, -0.251, 0)
        # => phi = arccos(-0.833) ~ 146.4 degrees (2.56 rad)
        phi = 2.56

        angular_term = np.cos(theta)**2 * np.exp(-(np.pi - theta)**6) * np.cos(phi - np.radians(109.5))**2

        # energy = 8 * [5 * (2.80 / d) ** 12 - 6 * (2.80 / d) ** 10] * angular_term
        energy = 8 * (5 * (2.80 / d) ** 12 - 6 * (2.80 / d) ** 10) * angular_term

        assert energy - self.ENERGY_EPSILON <= bond[3] <= energy + self.ENERGY_EPSILON


    # test bifurcated hydrogen bond
    def test_bifurcated_hbond(self, mock_protein_graph_asn_glu, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_asn_glu

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        bonds = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]
        assert len(bonds) == 2

        bond_hd21_oe1 = next(filter(lambda x: 'A0001-ASN:HD21' in x, bonds), None)
        assert bond_hd21_oe1 is not None
        assert 'A0002-GLU:OE1' in bond_hd21_oe1
        assert bond_hd21_oe1[2] == 5
        assert bond_hd21_oe1[3] < 0

        bond_hd22_oe2 = next(filter(lambda x: 'A0001-ASN:HD22' in x, bonds), None)
        assert bond_hd22_oe2 is not None
        assert 'A0002-GLU:OE2' in bond_hd22_oe2
        assert bond_hd22_oe2[2] == 5
        assert bond_hd22_oe2[3] < 0


    # test sulphur participated hydrogen bond
    def test_detect_sulphur_type_hbond(self, mock_protein_graph_lys_cyh, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_cyh

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value][0]
        assert 'A0001-LYS:HZ1' in bond
        assert 'A0002-CYH:SG' in bond
        assert bond[2] == 5
        assert bond[3] < 0

        # sulphur node in hydrogen bond allows for farther bond lengths
        assert 2.6 < protein_graph.graphs['full'].edges['A0001-LYS:HZ1', 'A0002-CYH:SG']['bond_length'] <= 3.0


    # detect single salt-bridge (no energy cut-off)
    def test_detect_single_salt_bridge(self, mock_protein_graph_lys_glu, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_glu

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value][0]
        assert 'A0001-LYS:NZ' in bond
        assert 'A0002-GLU:OE2' in bond
        assert bond[2] == 5
        assert bond[3] < 0
        assert InteractionType.SALT_BRIDGE.value in protein_graph.graphs['full'].edges['A0001-LYS:NZ', 'A0002-GLU:OE2']['kind']


    # detect multiple salt-bridges with correct energy values and angles
    def test_detect_multiple_salt_bridges(self, mock_protein_graph_arg_glu, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_arg_glu

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        bonds = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]
        assert len(bonds) == 4

        for bond in bonds:
            assert 'A0001-ARG:NH1' in bond or 'A0001-ARG:NH2' in bond
            assert 'A0002-GLU:OE1' in bond or 'A0002-GLU:OE2' in bond
            assert bond[2] == 5
            assert bond[3] < 0

        assert InteractionType.SALT_BRIDGE.value in protein_graph.graphs['full'].edges['A0001-ARG:NH1', 'A0002-GLU:OE1']['kind']
        assert InteractionType.SALT_BRIDGE.value in protein_graph.graphs['full'].edges['A0001-ARG:NH2', 'A0002-GLU:OE2']['kind']


    # raise error on missing atoms for calculation (e.g., third atom for calculating gamma angle for donor)
    def test_error_on_missing_donor_base_for_sp2(self, mock_protein_graph_gly_gly_no_donor_base, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly_no_donor_base

        with pytest.raises(ValueError, match="sp2-center at A0001-GLY:N"):
            interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)


    # raise error on missing atoms for calculation (e.g., third atom for calculating gamma angle for acceptor)
    def test_error_on_missing_acceptor_base_for_sp2(self, mock_protein_graph_gly_gly_no_acceptor_base, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly_no_acceptor_base

        with pytest.raises(ValueError, match="sp2-center at A0002-GLY:O"):
            interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)


    # test the correct usage of the BAR_COUNT parameter
    def test_bar_count_parameter(self, mock_protein_graph_gly_gly, parameters_hbond_bars_3):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly

        interaction_service.apply_interactions(protein_graph, parameters_hbond_bars_3, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value][0]
        assert "A0001-GLY:H" in bond
        assert "A0002-GLY:O" in bond
        assert bond[2] == 3


    # correctly scale number of bars by energy (evaluate bar count and strong energy threshold, needs multiple runs)
    def test_scale_bar_count_2(self, mock_protein_graph_gly_gly, parameters_hbond_scaled_bars_minus5_minus0):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly

        interaction_service.apply_interactions(protein_graph, parameters_hbond_scaled_bars_minus5_minus0, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value][0]
        assert "A0001-GLY:H" in bond
        assert "A0002-GLY:O" in bond
        assert -1.75 - self.ENERGY_EPSILON <= bond[3] <= -1.75 + self.ENERGY_EPSILON
        # interval (-5,0] split into 4 parts
        assert bond[2] == 2


    def test_scale_bar_count_4(self, mock_protein_graph_gly_gly, parameters_hbond_scaled_bars_minus2_minus0):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly

        interaction_service.apply_interactions(protein_graph, parameters_hbond_scaled_bars_minus2_minus0, verbose=False)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value][0]
        assert "A0001-GLY:H" in bond
        assert "A0002-GLY:O" in bond
        assert -1.75 - self.ENERGY_EPSILON <= bond[3] <= -1.75 + self.ENERGY_EPSILON
        # interval (-2,0] split into 4 parts
        assert bond[2] == 4


    # evaluate correct logging of (1) no bonds (2) only h-bonds (3) only salt-bridges (4) both
    def test_logging_no_bonds(self,
            mock_protein_graph_gly_gly_large_distance,
            parameters_hbond_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly_large_distance

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_hbond_default)

        captured = capsys.readouterr()[1]

        assert re.search(r"INFO \| .+ 0 hydrogen bonds", captured)
        assert re.search(r"INFO \| .+ 0 salt-bridges", captured)


    def test_logging_single_hbond(self,
            mock_protein_graph_gly_gly,
            parameters_hbond_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_hbond_default)

        captured = capsys.readouterr()[1]

        assert re.search(r"INFO \| .+ 1 hydrogen bonds", captured)
        assert re.search(r"INFO \| .+ 0 salt-bridges", captured)


    def test_logging_single_salt_bridge(self,
            mock_protein_graph_lys_glu,
            parameters_hbond_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_glu

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_hbond_default)

        captured = capsys.readouterr()[1]

        assert re.search(r"INFO \| .+ 0 hydrogen bonds", captured)
        assert re.search(r"INFO \| .+ 1 salt-bridges", captured)


    def test_logging_single_hbond_and_two_salt_bridges(self,
            mock_protein_graph_gly_gly_and_arg_glu,
            parameters_hbond_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_gly_and_arg_glu

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_hbond_default)

        captured = capsys.readouterr()[1]

        assert re.search(r"INFO \| .+ 1 hydrogen bonds", captured)
        assert re.search(r"INFO \| .+ 4 salt-bridges", captured)

    # DNA & RNA

    def test_detect_dna_a_t_hbonds(self, mock_protein_graph_dna_a_t, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_dna_a_t

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        bonds = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]
        assert len(bonds) == 2

        bond_h61_o4 = next(filter(lambda x: 'A0001- DA:H61' in x, bonds), None)
        assert bond_h61_o4 is not None
        assert 'A0002- DT:O4' in bond_h61_o4
        assert bond_h61_o4[2] == 5
        assert bond_h61_o4[3] < 0

        bond_h3_n1 = next(filter(lambda x: 'A0002- DT:H3' in x, bonds), None)
        assert bond_h3_n1 is not None
        assert 'A0001- DA:N1' in bond_h3_n1
        assert bond_h3_n1[2] == 5
        assert bond_h3_n1[3] < 0

    def test_detect_dna_g_c_hbonds(self, mock_protein_graph_dna_g_c, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_dna_g_c

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        bonds = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]
        assert len(bonds) == 3

        bond_h1_n3 = next(filter(lambda x: 'A0001- DG:H1' in x, bonds), None)
        assert bond_h1_n3 is not None
        assert 'A0002- DC:N3' in bond_h1_n3
        assert bond_h1_n3[2] == 5
        assert bond_h1_n3[3] < 0

        bond_h22_o2 = next(filter(lambda x: 'A0001- DG:H22' in x, bonds), None)
        assert bond_h22_o2 is not None
        assert 'A0002- DC:O2' in bond_h22_o2
        assert bond_h22_o2[2] == 5
        assert bond_h22_o2[3] < 0

        bond_h42_o6 = next(filter(lambda x: 'A0002- DC:H42' in x, bonds), None)
        assert bond_h42_o6 is not None
        assert 'A0001- DG:O6' in bond_h42_o6
        assert bond_h42_o6[2] == 5
        assert bond_h42_o6[3] < 0

    def test_detect_rna_a_u_hbonds(self, mock_protein_graph_rna_a_u, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_rna_a_u

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        bonds = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]
        assert len(bonds) == 2

        bond_h61_o4 = next(filter(lambda x: 'A0001-  A:H61' in x, bonds), None)
        assert bond_h61_o4 is not None
        assert 'A0002-  U:O4' in bond_h61_o4
        assert bond_h61_o4[2] == 5
        assert bond_h61_o4[3] < 0

        bond_h3_n1 = next(filter(lambda x: 'A0002-  U:H3' in x, bonds), None)
        assert bond_h3_n1 is not None
        assert 'A0001-  A:N1' in bond_h3_n1
        assert bond_h3_n1[2] == 5
        assert bond_h3_n1[3] < 0

    def test_detect_rna_g_c_hbonds(self, mock_protein_graph_rna_g_c, parameters_hbond_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_rna_g_c

        interaction_service.apply_interactions(protein_graph, parameters_hbond_default, verbose=False)

        bonds = protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]
        assert len(bonds) == 3

        bond_h1_n3 = next(filter(lambda x: 'A0001-  G:H1' in x, bonds), None)
        assert bond_h1_n3 is not None
        assert 'A0002-  C:N3' in bond_h1_n3
        assert bond_h1_n3[2] == 5
        assert bond_h1_n3[3] < 0

        bond_h22_o2 = next(filter(lambda x: 'A0001-  G:H22' in x, bonds), None)
        assert bond_h22_o2 is not None
        assert 'A0002-  C:O2' in bond_h22_o2
        assert bond_h22_o2[2] == 5
        assert bond_h22_o2[3] < 0

        bond_h42_o6 = next(filter(lambda x: 'A0002-  C:H42' in x, bonds), None)
        assert bond_h42_o6 is not None
        assert 'A0001-  G:O6' in bond_h42_o6
        assert bond_h42_o6[2] == 5
        assert bond_h42_o6[3] < 0