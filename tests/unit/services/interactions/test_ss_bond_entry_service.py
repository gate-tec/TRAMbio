import re
import time

import pytest

from TRAMbio.services import InteractionServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import PdbEntryInteractionParameter, GeneralWorkflowParameter
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import GraphKey
from tests.conftest import inject_pytest_logger
from tests.util.graphs import construct_protein_graph_base, graph_css


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "SSBondEntryService"

    @pytest.fixture
    def parameters_no_ss_bonds(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(PdbEntryInteractionParameter.SSBOND_INCLUDE.value, False)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_ss_bond_default(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(PdbEntryInteractionParameter.SSBOND_INCLUDE.value, True)
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
        len_coords = len(coords1) + 1 + len(coords2)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1],
            ["SSBOND", f"{1:4d} CYS A{1:5d}    CYS A{2:5d}{' ' * 26}1555   1555 {2.04:5.2f}", len_coords + 1],
            ["CONECT", f"{1:5d}{len(coords1) + 2:5d}", len_coords + 2]
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
    def mock_protein_graph_css_css_no_record(self):
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
    def mock_protein_graph_css_css_wrong_chain(self):
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_css(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_css(len(coords1) + 2, 2, 0.0, 2.04, rotations=[(180, 'Z')])
        len_coords = len(coords1) + 1 + len(coords2)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1],
            ["SSBOND", f"{1:4d} CYS B{1:5d}    CYS A{2:5d}{' ' * 26}1555   1555 {2.04:5.2f}", len_coords + 1],
            ["CONECT", f"{50:5d}{53:5d}", len_coords + 2]
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

    # correctly skip detection by parameter PdbEntryInteractionParameter.SSBOND_INCLUDE.value == False
    def test_skip_eval(self, mock_protein_graph_css_css, parameters_no_ss_bonds):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_no_ss_bonds)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0


    # insert no ss-bond with missing record
    def test_place_no_ss_bonds(self, mock_protein_graph_css_css_no_record, parameters_ss_bond_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css_no_record

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_ss_bond_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0


    # insert single ss-bond from record
    def test_place_single_ss_bond(self, mock_protein_graph_css_css, parameters_ss_bond_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_ss_bond_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 1
        assert InteractionType.SS_BOND.value in protein_graph.graphs['full'].edges['A0001-CSS:SG', 'A0002-CSS:SG']['kind']


    # correctly handle missing atoms for record
    def test_detect_missing_atoms(self, mock_protein_graph_css_css_wrong_chain, parameters_ss_bond_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css_wrong_chain

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_ss_bond_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0


    # evaluate correct logging of (1) no bonds (2) inserted bonds
    def test_logging_no_ss_bonds(self,
            mock_protein_graph_css_css_no_record,
            parameters_ss_bond_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css_no_record

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_ss_bond_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is None


    def test_logging_single_ss_bond(self,
            mock_protein_graph_css_css,
            parameters_ss_bond_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_ss_bond_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 1
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is not None
        assert re.search(r"Added .* SSBOND", info_line)


    def test_logging_warning_missing_atom_for_record(self,
            mock_protein_graph_css_css_wrong_chain,
            parameters_ss_bond_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_css_css_wrong_chain

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_ss_bond_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0
        warning_line = next(filter(lambda x: x.startswith("WARNING"), str(captured).split("\n")), None)
        assert warning_line is not None
        assert re.search(r"Unable to find", warning_line)
        assert re.search(r"SSBOND", warning_line)
