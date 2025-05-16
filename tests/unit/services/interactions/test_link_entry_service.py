import re
import time

import pytest

from TRAMbio.services import InteractionServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import PdbEntryInteractionParameter, GeneralWorkflowParameter
from TRAMbio.util.structure_library.graph_struct import GraphKey
from tests.conftest import inject_pytest_logger
from tests.util.graphs import construct_protein_graph_base, graph_gly_acceptor


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "LinkEntryService"

    @pytest.fixture
    def parameters_no_link(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(PdbEntryInteractionParameter.LINK_INCLUDE.value, False)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_link_default(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(PdbEntryInteractionParameter.LINK_INCLUDE.value, True)
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
    def mock_protein_graph_gly_na(self):
        coords, hm, cov, pe, pm, pb = graph_gly_acceptor(1, 1)
        others = [
            ["TER", f"{len(coords) + 1:5d}", len(coords) + 1],
            ["LINK", f"       O   GLY A{1:4d}{' ' * 16}NA    NA A{2:4d}     1555   1555 {2.98:5.2f}", len(coords) + 3],
            ["CONECT", f"{3:5d}{len(coords) + 2:5d}", len(coords) + 4]
        ]
        coords += [
            ["HETATM", len(coords) + 2, "A", 2, "", " NA", "NA", f"A{2:04d}- NA:NA", "NA", "",
             2.98, 0.0, 0.0, len(coords) + 2]
        ]

        yield construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=hm,
            cov=cov,
            pe=pe,
            pm=pm,
            pb=pb
        )


    @pytest.fixture
    def mock_protein_graph_gly_na_no_record(self):
        coords, hm, cov, pe, pm, pb = graph_gly_acceptor(1, 1)
        others = [
            ["TER", f"{len(coords) + 1:5d}", len(coords) + 1]
        ]
        coords += [
            ["HETATM", len(coords) + 2, "A", 2, "", " NA", "NA", f"A{2:04d}- NA:NA", "NA", "",
             2.98, 0.0, 0.0, len(coords) + 2]
        ]

        yield construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=hm,
            cov=cov,
            pe=pe,
            pm=pm,
            pb=pb
        )


    @pytest.fixture
    def mock_protein_graph_gly_na_wrong_chain(self):
        coords, hm, cov, pe, pm, pb = graph_gly_acceptor(1, 1)
        others = [
            ["TER", f"{len(coords) + 1:5d}", len(coords) + 1],
            ["LINK", f"       O   GLY B{3:4d}{' ' * 16}NA    NA B{4:4d}     1555   1555 {2.98:5.2f}", len(coords) + 3],
            ["CONECT", f"{3:5d}{len(coords) + 2:5d}", len(coords) + 4],
        ]
        coords += [
            ["HETATM", len(coords) + 2, "A", 2, "", " NA", "NA", f"A{2:04d}- NA:NA", "NA", "",
             2.98, 0.0, 0.0, len(coords) + 2]
        ]

        yield construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=hm,
            cov=cov,
            pe=pe,
            pm=pm,
            pb=pb
        )


#############
# Tests #####
#############

class TestApplyInteractions(TestMockGraphs, TestParameters):

    # correctly skip detection by parameter PdbEntryInteractionParameter.LINK_INCLUDE.value == False
    def test_skip_eval(self, mock_protein_graph_gly_na, parameters_no_link):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_na

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_no_link)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0


    # insert no link-bond with missing record
    def test_place_no_link(self, mock_protein_graph_gly_na_no_record, parameters_link_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_na_no_record

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_link_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0


    # insert single link-bond from record
    def test_place_single_link(self, mock_protein_graph_gly_na, parameters_link_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_na

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_link_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 1


    # correctly handle faulty record
    # TODO

    # correctly handle missing atoms for record
    def test_detect_missing_atoms(self, mock_protein_graph_gly_na_wrong_chain, parameters_link_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_na_wrong_chain

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_link_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0


    # evaluate correct logging of (1) no bonds (2) inserted bonds
    def test_logging_no_link(self,
            mock_protein_graph_gly_na_no_record,
            parameters_link_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_na_no_record

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_link_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is None


    def test_logging_single_link(self,
            mock_protein_graph_gly_na,
            parameters_link_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_na

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_link_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 1
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is not None
        assert re.search(r"Added 1 .* LINK", info_line)


    def test_logging_warning_missing_atom_for_record(self,
            mock_protein_graph_gly_na_wrong_chain,
            parameters_link_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_gly_na_wrong_chain

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_link_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0
        warning_line = next(filter(lambda x: x.startswith("WARNING"), str(captured).split("\n")), None)
        assert warning_line is not None
        assert re.search(r"Unable to find", warning_line)
        assert re.search(r"LINK", warning_line)
