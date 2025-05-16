import pytest

import time
import re
from TRAMbio.services import InteractionServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import DisulphideBridgeParameter, GeneralWorkflowParameter
from TRAMbio.util.structure_library.graph_struct import GraphKey
from tests.conftest import inject_pytest_logger
from tests.util.graphs import construct_protein_graph_base, graph_css, graph_cyh_acceptor


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "DisulphideBridgeService"

    @pytest.fixture
    def parameters_no_disulphide_bridges(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(DisulphideBridgeParameter.INCLUDE.value, False)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_disulphide_bridge_default(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(DisulphideBridgeParameter.INCLUDE.value, True)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_disulphide_bridge_unbounded(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(DisulphideBridgeParameter.INCLUDE.value, True)
        registry.set_parameter(DisulphideBridgeParameter.CUTOFF_DISTANCE.value, 5000)
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
    def mock_protein_graph_css_css(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_css(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_css(len(coords1) + 2, 2, 0.0, 2.04, rotations=[(180, 'Z')])
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
    def mock_protein_graph_cyh_css(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_cyh_acceptor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_css(len(coords1) + 2, 2, 0.0, 2.04, rotations=[(180, 'Z')])
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
    def mock_protein_graph_css_css_non_standard_distance(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_css(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_css(len(coords1) + 2, 2, 0.5, 2.04, rotations=[(180, 'Z')])
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
    def mock_protein_graph_css_css_large_distance(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_css(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_css(len(coords1) + 2, 2, 6.0, 2.04, rotations=[(180, 'Z')])
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

    # correctly skip detection by parameter DisulphideBridgeParameter.INCLUDE.value == False
    def test_skip_eval(self, mock_protein_graph_css_css, parameters_no_disulphide_bridges):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css

        interaction_service.apply_interactions(protein_graph, parameters_no_disulphide_bridges)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value]) == 0


    # detect no disulphide bridge
    def test_detect_no_disulphide_bridge(self, mock_protein_graph_cyh_css, parameters_disulphide_bridge_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_cyh_css

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_disulphide_bridge_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0


    # detect single disulphide bridge
    def test_detect_single_disulphide_bridge(self, mock_protein_graph_css_css, parameters_disulphide_bridge_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_disulphide_bridge_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 1


    # test the correct usage of the CUT_OFF parameter
    def test_cut_off_parameter(self, mock_protein_graph_css_css_large_distance, parameters_disulphide_bridge_unbounded):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css_large_distance

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_disulphide_bridge_unbounded)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 1


    # evaluate correct logging of (1) no bonds (2) detected bonds
    def test_logging_no_disulphide_bridges(self,
            mock_protein_graph_cyh_css,
            parameters_disulphide_bridge_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_cyh_css

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_disulphide_bridge_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is None


    def test_logging_single_disulphide_bridge(self,
            mock_protein_graph_css_css,
            parameters_disulphide_bridge_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_disulphide_bridge_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 1
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is not None
        assert re.search(r"1 disulphide", info_line)


    # evaluate correct debug logging
    def test_logging_debug_unusual_distance(self,
            mock_protein_graph_css_css_non_standard_distance,
            parameters_disulphide_bridge_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css_non_standard_distance

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_disulphide_bridge_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 1
        debug_line = next(filter(lambda x: x.startswith("DEBUG"), str(captured).split("\n")), None)
        assert debug_line is not None
        assert re.search(r"A0001-CSS:SG .* A0002-CSS:SG", debug_line)
        assert re.search(r"doesn't match .* instead of 2.04", debug_line)


    def test_logging_debug_not_enough_suitable_cysteine_residues(self,
            mock_protein_graph_cyh_css,
            parameters_disulphide_bridge_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_cyh_css

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_disulphide_bridge_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0
        debug_line = next(filter(lambda x: x.startswith("DEBUG"), str(captured).split("\n")), None)
        assert debug_line is not None
        assert re.search(r"Cannot add disulphide", debug_line)


    def test_logging_debug_no_valid_interaction(self,
            mock_protein_graph_css_css_large_distance,
            parameters_disulphide_bridge_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css_large_distance

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_disulphide_bridge_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0
        debug_line = next(filter(lambda x: x.startswith("DEBUG"), str(captured).split("\n")), None)
        assert debug_line is not None
        assert re.search(r"No valid disulphide", debug_line)
