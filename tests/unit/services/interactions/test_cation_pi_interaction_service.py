import re
import time

import pytest

from TRAMbio.services import InteractionServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import CationPiInteractionParameter, GeneralWorkflowParameter
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import GraphKey
from tests.conftest import inject_pytest_logger
from tests.util.graphs import construct_protein_graph_base, graph_phe_aromatic_xy_plane, graph_lys_donor


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "CationPiInteractionService"

    @pytest.fixture
    def parameters_no_cation_pi(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(CationPiInteractionParameter.INCLUDE.value, False)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_cation_pi_default(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(CationPiInteractionParameter.INCLUDE.value, True)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_cation_pi_unbounded(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(CationPiInteractionParameter.INCLUDE.value, True)
        registry.set_parameter(CationPiInteractionParameter.CUTOFF_DISTANCE.value, 5000)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_cation_pi_bars_1(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(CationPiInteractionParameter.INCLUDE.value, True)
        registry.set_parameter(CationPiInteractionParameter.BAR_COUNT.value, 1)
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
    def mock_protein_graph_lys_phe(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_lys_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_phe_aromatic_xy_plane(len(coords1) + 2, 2, 4.5, rotations=[(-90, 'Y')])  # yz-plane
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
    def mock_protein_graph_lys_phe_large_distance(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_lys_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_phe_aromatic_xy_plane(len(coords1) + 2, 2, 15.0, rotations=[(-90, 'Y')])  # yz-plane
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
    def mock_protein_graph_lys_phe_wrong_angle(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_lys_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_phe_aromatic_xy_plane(
            len(coords1) + 2, 2, 4.5, rotations=[(-90, 'Y'), (65, 'Z')]
        )  # yz-plane but then rotation along z-axis (more than 60 degrees deviation)
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


#############
# Tests #####
#############

class TestApplyInteractions(TestMockGraphs, TestParameters):

    # correctly skip detection by parameter CationPiInteractionParameter.INCLUDE.value == False
    def test_skip_eval(self, mock_protein_graph_lys_phe, parameters_no_cation_pi):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_phe

        interaction_service.apply_interactions(protein_graph, parameters_no_cation_pi)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 0


    # detect no cation-pi interactions
    def test_detect_no_cation_pi_by_distance(self, mock_protein_graph_lys_phe_large_distance, parameters_cation_pi_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_phe_large_distance

        interaction_service.apply_interactions(protein_graph, parameters_cation_pi_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 0


    def test_detect_no_cation_pi_by_angle(self, mock_protein_graph_lys_phe_wrong_angle, parameters_cation_pi_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_phe_wrong_angle

        interaction_service.apply_interactions(protein_graph, parameters_cation_pi_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 0


    # detect single cation-pi interaction
    def test_detect_single_cation_pi(self, mock_protein_graph_lys_phe, parameters_cation_pi_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_phe

        interaction_service.apply_interactions(protein_graph, parameters_cation_pi_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value][0]
        assert "A0001-LYS:NZ" in bond
        assert "A0002-PHE:CG" in bond
        assert bond[2] == 3

        assert InteractionType.CATION_PI.value in protein_graph.graphs['full'].edges['A0001-LYS:NZ', 'A0002-PHE:CG']['kind']


    # test the correct usage of the BAR_COUNT parameter
    def test_bar_count_parameter(self, mock_protein_graph_lys_phe, parameters_cation_pi_bars_1):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_phe

        interaction_service.apply_interactions(protein_graph, parameters_cation_pi_bars_1)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value][0]
        assert "A0001-LYS:NZ" in bond
        assert "A0002-PHE:CG" in bond
        assert bond[2] == 1


    # test the correct usage of the CUT_OFF parameter
    def test_cut_off_parameter(self, mock_protein_graph_lys_phe_large_distance, parameters_cation_pi_unbounded):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_phe_large_distance

        interaction_service.apply_interactions(protein_graph, parameters_cation_pi_unbounded)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1


    # evaluate correct logging of (1) no bonds (2) detected bonds
    def test_logging_debug_no_cation_pi_interactions(self,
            mock_protein_graph_lys_phe_large_distance,
            parameters_cation_pi_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_phe_large_distance

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_cation_pi_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 0
        debug_line = next(filter(lambda x: x.startswith("DEBUG"), str(captured).split("\n")), None)
        assert debug_line is not None
        assert re.search(r"0 cation-pi", debug_line)


    def test_logging_single_cation_pi_interaction(self,
            mock_protein_graph_lys_phe,
            parameters_cation_pi_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_lys_phe

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_cation_pi_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is not None
        assert re.search(r"1 cation-pi", info_line)