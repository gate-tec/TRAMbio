import re
import time

import pytest

from TRAMbio.services import InteractionServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import AromaticInteractionParameter, GeneralWorkflowParameter
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import GraphKey
from tests.conftest import inject_pytest_logger
from tests.util.graphs import construct_protein_graph_base, graph_phe_aromatic_xy_plane


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "AromaticInteractionService"

    @pytest.fixture
    def parameters_no_aromatic_interactions(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(AromaticInteractionParameter.INCLUDE.value, False)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id


    @pytest.fixture
    def parameters_aromatic_interactions_default(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(AromaticInteractionParameter.INCLUDE.value, True)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id


    @pytest.fixture
    def parameters_aromatic_interactions_unbounded(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(AromaticInteractionParameter.INCLUDE.value, True)
        registry.set_parameter(AromaticInteractionParameter.CUTOFF_DISTANCE_PI.value, 5000)
        registry.set_parameter(AromaticInteractionParameter.CUTOFF_DISTANCE_T.value, 5000)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id


    @pytest.fixture
    def parameters_aromatic_interactions_bars_1(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(AromaticInteractionParameter.INCLUDE.value, True)
        registry.set_parameter(AromaticInteractionParameter.BAR_COUNT.value, 1)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id


#################################
# No Mock Services required #####
#################################

class TestMockGraphs:

    @pytest.fixture
    def mock_protein_graph_yz_plane_pi(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_phe_aromatic_xy_plane(1, 1, rotations=[(-90, 'Y')])  # yz-plane
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_phe_aromatic_xy_plane(len(coords1) + 2, 2,5.0, rotations=[(-90, 'Y')])  # yz-plane
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
    def mock_protein_graph_yz_plane_pi_large_distance(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_phe_aromatic_xy_plane(1, 1, rotations=[(-90, 'Y')])  # yz-plane
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_phe_aromatic_xy_plane(len(coords1) + 2, 2,15.0, rotations=[(-90, 'Y')])  # yz-plane
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
    def mock_protein_graph_yz_plane_t(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_phe_aromatic_xy_plane(1, 1, rotations=[(-90, 'Y')])  # yz-plane
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_phe_aromatic_xy_plane(len(coords1) + 2, 2, 4.0)
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
    def mock_protein_graph_yz_plane_t_large_distance(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_phe_aromatic_xy_plane(1, 1, rotations=[(-90, 'Y')])  # yz-plane
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_phe_aromatic_xy_plane(len(coords1) + 2, 2, 15.0)
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

    # correctly skip detection by parameter AromaticInteractionParameter.INCLUDE.value == False
    def test_skip_eval(self, mock_protein_graph_yz_plane_pi, parameters_no_aromatic_interactions):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_yz_plane_pi

        interaction_service.apply_interactions(protein_graph, parameters_no_aromatic_interactions)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 0


    # detect no aromatic bonds
    def test_detect_no_aromatic_interaction(self, mock_protein_graph_yz_plane_pi_large_distance, parameters_aromatic_interactions_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_yz_plane_pi_large_distance

        interaction_service.apply_interactions(protein_graph, parameters_aromatic_interactions_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 0


    # detect single aromatic bond (Pi-Stacking)
    def test_detect_single_pi_stacking(self, mock_protein_graph_yz_plane_pi, parameters_aromatic_interactions_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_yz_plane_pi

        interaction_service.apply_interactions(protein_graph, parameters_aromatic_interactions_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value][0]
        assert 'A0001-PHE:CG' in bond
        assert 'A0002-PHE:CG' in bond
        assert bond[2] == 3

        assert InteractionType.PI_STACKING.value in protein_graph.graphs['full'].edges['A0001-PHE:CG', 'A0002-PHE:CG']['kind']


    # detect single aromatic bond (T-Stacking)
    def test_detect_single_t_stacking(self, mock_protein_graph_yz_plane_t, parameters_aromatic_interactions_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_yz_plane_t

        interaction_service.apply_interactions(protein_graph, parameters_aromatic_interactions_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value][0]
        assert 'A0001-PHE:CG' in bond
        assert 'A0002-PHE:CG' in bond
        assert bond[2] == 3

        assert InteractionType.T_STACKING.value in protein_graph.graphs['full'].edges['A0001-PHE:CG', 'A0002-PHE:CG']['kind']


    # test the correct usage of the BAR_COUNT parameter
    def test_bar_count_parameter(self, mock_protein_graph_yz_plane_pi, parameters_aromatic_interactions_bars_1):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_yz_plane_pi

        interaction_service.apply_interactions(protein_graph, parameters_aromatic_interactions_bars_1)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value][0]
        assert 'A0001-PHE:CG' in bond
        assert "A0002-PHE:CG" in bond
        assert bond[2] == 1


    # test the correct usage of the CUT_OFF parameter (Pi-Stacking)
    def test_cut_off_parameter_pi(self, mock_protein_graph_yz_plane_pi_large_distance, parameters_aromatic_interactions_unbounded):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_yz_plane_pi_large_distance

        interaction_service.apply_interactions(protein_graph, parameters_aromatic_interactions_unbounded)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1


    # test the correct usage of the CUT_OFF parameter (T-Stacking)
    def test_cut_off_parameter_t(self, mock_protein_graph_yz_plane_t_large_distance, parameters_aromatic_interactions_unbounded):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_yz_plane_t_large_distance

        interaction_service.apply_interactions(protein_graph, parameters_aromatic_interactions_unbounded)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1


    # evaluate correct logging of (1) no bonds (2) detected bonds
    def test_logging_debug_no_aromatic_interactions(self,
            mock_protein_graph_yz_plane_pi_large_distance,
            parameters_aromatic_interactions_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_yz_plane_pi_large_distance

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_aromatic_interactions_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 0
        info_line = next(filter(lambda x: x.startswith("DEBUG"), str(captured).split("\n")), None)
        assert info_line is not None
        assert re.search(r"0 aromatic", info_line)


    def test_logging_single_aromatic_interaction(self,
            mock_protein_graph_yz_plane_pi,
            parameters_aromatic_interactions_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_yz_plane_pi

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_aromatic_interactions_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is not None
        assert re.search(r"1 aromatic", info_line)
