import re
import time

import pytest

from TRAMbio.services import InteractionServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import HydrophobicInteractionParameter, GeneralWorkflowParameter
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import GraphKey
from tests.conftest import inject_pytest_logger
from tests.util.graphs import construct_protein_graph_base, graph_phe_aromatic_xy_plane, graph_cyh_acceptor


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "HydrophobicSurfaceInteractionService"

    @pytest.fixture
    def parameters_no_hydrophobics(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrophobicInteractionParameter.INCLUDE.value, False)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_hydrophobics_default(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrophobicInteractionParameter.INCLUDE.value, True)
        registry.set_parameter(HydrophobicInteractionParameter.POTENTIAL.value, False)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_hydrophobics_bars_1(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrophobicInteractionParameter.INCLUDE.value, True)
        registry.set_parameter(HydrophobicInteractionParameter.POTENTIAL.value, False)
        registry.set_parameter(HydrophobicInteractionParameter.BAR_COUNT.value, 1)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_hydrophobics_unbounded(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrophobicInteractionParameter.INCLUDE.value, True)
        registry.set_parameter(HydrophobicInteractionParameter.POTENTIAL.value, False)
        registry.set_parameter(HydrophobicInteractionParameter.SURFACE_CUTOFF_DISTANCE.value, 5000)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_hydrophobics_multiplex(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(HydrophobicInteractionParameter.INCLUDE.value, True)
        registry.set_parameter(HydrophobicInteractionParameter.POTENTIAL.value, False)
        registry.set_parameter(HydrophobicInteractionParameter.MINIMAL_LENGTH.value, False)
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
    def mock_protein_graph_phe(self):
        coords, hm, cov, pe, pm, pb = graph_phe_aromatic_xy_plane(1, 1, rotations=[(-90, 'Y')]) # yz-plane

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
    def mock_protein_graph_cyh_cyh(self):
        # VDW_RADIUS(S) = 1.8
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_cyh_acceptor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_cyh_acceptor(len(coords1) + 2, 2, -3.0, rotations=[(180, 'Z')])
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
    def mock_protein_graph_cyh_cyh_large_distance(self):
        # VDW_RADIUS(S) = 1.8
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_cyh_acceptor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_cyh_acceptor(len(coords1) + 2, 2, -18.0, rotations=[(180, 'Z')])
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
    def mock_protein_graph_phe_cyh(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_phe_aromatic_xy_plane(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_cyh_acceptor(len(coords1) + 2, 2, 4.0)
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
    def mock_protein_graph_phe_cyh_multiplex(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_phe_aromatic_xy_plane(1, 1, rotations=[(-90, 'Y')])
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_cyh_acceptor(len(coords1) + 2, 2, 2.5, z_offset=0.2)
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
    def mock_protein_graph_phe_phe(self):
        # all six pairs of carbon atoms at perfect distance of 3.4 Å = 2 * VDW_RADIUS(C)
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_phe_aromatic_xy_plane(1, 1, rotations=[(-90, 'Y')]) # yz-plane
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_phe_aromatic_xy_plane(len(coords1) + 2, 2, 3.4, rotations=[(-90, 'Y')]) # yz-plane
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

    # correctly skip detection by parameter HydrophobicInteractionParameter.INCLUDE.value == False
    def test_skip_eval(self, mock_protein_graph_phe, parameters_no_hydrophobics):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_phe

        interaction_service.apply_interactions(protein_graph, parameters_no_hydrophobics)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 0


    # detect no hydrophobic interactions
    def test_no_hydrophobics(self, mock_protein_graph_phe, parameters_hydrophobics_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_phe

        interaction_service.apply_interactions(protein_graph, parameters_hydrophobics_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 0


    # detect single hydrophobic interaction
    def test_detect_single_hydrophobic_sulphur_sulphur(self, mock_protein_graph_cyh_cyh, parameters_hydrophobics_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_cyh_cyh

        interaction_service.apply_interactions(protein_graph, parameters_hydrophobics_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value][0]
        assert "A0001-CYH:SG" in bond
        assert "A0002-CYH:SG" in bond
        assert bond[2] == 3

        assert InteractionType.HYDROPHOBIC.value in protein_graph.graphs['full'].edges["A0001-CYH:SG", "A0002-CYH:SG"]["kind"]


    def test_detect_single_hydrophobic_sulphur_carbon(self, mock_protein_graph_phe_cyh, parameters_hydrophobics_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_phe_cyh

        interaction_service.apply_interactions(protein_graph, parameters_hydrophobics_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value][0]
        assert "A0001-PHE:CG" in bond
        assert "A0002-CYH:SG" in bond
        assert bond[2] == 3

        assert InteractionType.HYDROPHOBIC.value in protein_graph.graphs['full'].edges["A0001-PHE:CG", "A0002-CYH:SG"]["kind"]


    def test_detect_single_hydrophobic_sulphur_carbon_minimal_length(self, 
            mock_protein_graph_phe_cyh_multiplex,
            parameters_hydrophobics_default
    ):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_phe_cyh_multiplex

        interaction_service.apply_interactions(protein_graph, parameters_hydrophobics_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value][0]
        assert "A0001-PHE:CG" in bond
        assert "A0002-CYH:SG" in bond
        assert bond[2] == 3

        assert InteractionType.HYDROPHOBIC.value in protein_graph.graphs['full'].edges["A0001-PHE:CG", "A0002-CYH:SG"]["kind"]


    # detect multiple hydrophobic interactions
    def test_detect_multiple_hydrophobics(self, mock_protein_graph_phe_phe, parameters_hydrophobics_default):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_phe_phe

        interaction_service.apply_interactions(protein_graph, parameters_hydrophobics_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 6


    # test the correct usage of the BAR_COUNT parameter
    def test_bar_count_parameter(self, mock_protein_graph_cyh_cyh, parameters_hydrophobics_bars_1):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_cyh_cyh

        interaction_service.apply_interactions(protein_graph, parameters_hydrophobics_bars_1)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        bond = protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value][0]
        assert "A0001-CYH:SG" in bond
        assert "A0002-CYH:SG" in bond
        assert bond[2] == 1


    # test the correct usage of the SURFACE_CUT_OFF parameter
    def test_cut_off_parameter(self, mock_protein_graph_cyh_cyh_large_distance, parameters_hydrophobics_unbounded):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_cyh_cyh_large_distance

        interaction_service.apply_interactions(protein_graph, parameters_hydrophobics_unbounded)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 2


    # test the correct usage of the MINIMAL_LENGTH parameter
    def test_minimal_length_parameter(self, mock_protein_graph_phe_cyh_multiplex, parameters_hydrophobics_multiplex):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_phe_cyh_multiplex

        interaction_service.apply_interactions(protein_graph, parameters_hydrophobics_multiplex)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 6
        for bond in protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]:
            assert "A0002-CYH:SG" in bond


    # evaluate correct logging of (1) no bonds (2) detected bonds
    def test_logging_debug_no_hydrophobic_interactions(self,
            mock_protein_graph_cyh_cyh_large_distance,
            parameters_hydrophobics_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_cyh_cyh_large_distance

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_hydrophobics_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 0
        info_line = next(filter(lambda x: x.startswith("DEBUG"), str(captured).split("\n")), None)
        assert info_line is not None
        assert re.search(r"0 hydrophobic", info_line)


    def test_logging_single_hydrophobic_interaction(self,
            mock_protein_graph_cyh_cyh,
            parameters_hydrophobics_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.NON_COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_cyh_cyh

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_hydrophobics_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.NON_COVALENT_EDGES.value]) == 1
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is not None
        assert re.search(r"1 hydrophobic", info_line)
