import re
import time

import pytest

from TRAMbio.services import InteractionServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import PdbEntryInteractionParameter, GeneralWorkflowParameter
from TRAMbio.util.structure_library.graph_struct import GraphKey
from tests.conftest import inject_pytest_logger
from tests.util.graphs import construct_protein_graph_base


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "ConectEntryService"

    @pytest.fixture
    def parameters_no_conect(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(PdbEntryInteractionParameter.CONECT_INCLUDE.value, False)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield param_id

    @pytest.fixture
    def parameters_conect_default(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(PdbEntryInteractionParameter.CONECT_INCLUDE.value, True)
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
    def mock_protein_graph_methane(self):
        coords = [
            ["HETATM", 1, "A", 1, "", "CH4", "C", f"A{1:04d}-CH4:C", "C", "",
             0, 0, 0, 0],
            ["HETATM", 2, "A", 1, "", "CH4", "H1", f"A{1:04d}-CH4:H1", "H", "",
             0, 0, 0, 0],
            ["HETATM", 3, "A", 1, "", "CH4", "H2", f"A{1:04d}-CH4:H2", "H", "",
             0, 0, 0, 0],
            ["HETATM", 4, "A", 1, "", "CH4", "H3", f"A{1:04d}-CH4:H3", "H", "",
             0, 0, 0, 0],
            ["HETATM", 5, "A", 1, "", "CH4", "H4", f"A{1:04d}-CH4:H4", "H", "",
             0, 0, 0, 0],
        ]
        others = [
            ["CONECT", f"{1:5d}{2:5d}{3:5d}{4:5d}{5:5d}", len(coords) + 1]
        ]
        pm = {
            f"A{1:04d}-CH4:C": 6,
            f"A{1:04d}-CH4:H1": 6,
            f"A{1:04d}-CH4:H2": 6,
            f"A{1:04d}-CH4:H3": 6,
            f"A{1:04d}-CH4:H4": 6,
        }

        yield construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=[],
            cov=[],
            pe=[],
            pm=pm,
            pb=[]
        )


    @pytest.fixture
    def mock_protein_graph_water_single_bond(self):
        coords = [
            ["HETATM", 1, "A", 1, "", "H2O", "O", f"A{1:04d}-H2O:O", "O", "",
             0, 0, 0, 0],
            ["HETATM", 2, "A", 1, "", "H2O", "H1", f"A{1:04d}-H2O:H1", "H", "",
             0, 0, 0, 0],
            ["HETATM", 3, "A", 1, "", "H2O", "H2", f"A{1:04d}-H2O:H2", "H", "",
             0, 0, 0, 0],
        ]
        others = [
            ["CONECT", f"{1:5d}{2:5d}", len(coords) + 1]
        ]
        pm = {
            f"A{1:04d}-H2O:O": 6,
            f"A{1:04d}-H2O:H1": 6,
            f"A{1:04d}-H2O:H2": 6,
        }

        yield construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=[],
            cov=[],
            pe=[],
            pm=pm,
            pb=[]
        )


    @pytest.fixture
    def mock_protein_graph_water_both_bonds(self):
        coords = [
            ["HETATM", 1, "A", 1, "", "H2O", "O", f"A{1:04d}-H2O:O", "O", "",
             0, 0, 0, 0],
            ["HETATM", 2, "A", 1, "", "H2O", "H1", f"A{1:04d}-H2O:H1", "H", "",
             0, 0, 0, 0],
            ["HETATM", 3, "A", 1, "", "H2O", "H2", f"A{1:04d}-H2O:H2", "H", "",
             0, 0, 0, 0],
        ]
        others = [
            ["CONECT", f"{1:5d}{2:5d}{3:5d}", len(coords) + 1]
        ]
        pm = {
            f"A{1:04d}-H2O:O": 6,
            f"A{1:04d}-H2O:H1": 6,
            f"A{1:04d}-H2O:H2": 6,
        }

        yield construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=[],
            cov=[],
            pe=[],
            pm=pm,
            pb=[]
        )


    @pytest.fixture
    def mock_protein_graph_methane_no_record(self):
        coords = [
            ["HETATM", 1, "A", 1, "", "CH4", "C", f"A{1:04d}-CH4:C", "C", "",
             0, 0, 0, 0],
            ["HETATM", 2, "A", 1, "", "CH4", "H1", f"A{1:04d}-CH4:H1", "H", "",
             0, 0, 0, 0],
            ["HETATM", 3, "A", 1, "", "CH4", "H2", f"A{1:04d}-CH4:H2", "H", "",
             0, 0, 0, 0],
            ["HETATM", 4, "A", 1, "", "CH4", "H3", f"A{1:04d}-CH4:H3", "H", "",
             0, 0, 0, 0],
            ["HETATM", 5, "A", 1, "", "CH4", "H4", f"A{1:04d}-CH4:H4", "H", "",
             0, 0, 0, 0],
        ]
        others = []
        pm = {
            f"A{1:04d}-CH4:C": 6,
            f"A{1:04d}-CH4:H1": 6,
            f"A{1:04d}-CH4:H2": 6,
            f"A{1:04d}-CH4:H3": 6,
            f"A{1:04d}-CH4:H4": 6,
        }

        yield construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=[],
            cov=[],
            pe=[],
            pm=pm,
            pb=[]
        )


    @pytest.fixture
    def mock_protein_graph_methane_wrong_atom_numbers(self):
        coords = [
            ["HETATM", 1, "A", 1, "", "CH4", "C", f"A{1:04d}-CH4:C", "C", "",
             0, 0, 0, 0],
            ["HETATM", 2, "A", 1, "", "CH4", "H1", f"A{1:04d}-CH4:H1", "H", "",
             0, 0, 0, 0],
            ["HETATM", 3, "A", 1, "", "CH4", "H2", f"A{1:04d}-CH4:H2", "H", "",
             0, 0, 0, 0],
            ["HETATM", 4, "A", 1, "", "CH4", "H3", f"A{1:04d}-CH4:H3", "H", "",
             0, 0, 0, 0],
            ["HETATM", 5, "A", 1, "", "CH4", "H4", f"A{1:04d}-CH4:H4", "H", "",
             0, 0, 0, 0],
        ]
        others = [
            ["CONECT", f"{6:5d}{7:5d}{8:5d}{9:5d}{10:5d}", len(coords) + 1]
        ]
        pm = {
            f"A{1:04d}-CH4:C": 6,
            f"A{1:04d}-CH4:H1": 6,
            f"A{1:04d}-CH4:H2": 6,
            f"A{1:04d}-CH4:H3": 6,
            f"A{1:04d}-CH4:H4": 6,
        }

        yield construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=[],
            cov=[],
            pe=[],
            pm=pm,
            pb=[]
        )


#############
# Tests #####
#############

class TestApplyInteractions(TestMockGraphs, TestParameters):

    # correctly skip detection by parameter PdbEntryInteractionParameter.CONECT_INCLUDE.value == False
    def test_skip_eval(self, mock_protein_graph_methane, parameters_no_conect):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_methane

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_no_conect)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0


    # insert no conect-bond with missing record
    def test_place_no_conect(self, mock_protein_graph_methane_no_record, parameters_conect_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_methane_no_record

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_conect_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0

    # insert single conect-bond from record
    def test_place_single_conect_bond(self, mock_protein_graph_water_single_bond, parameters_conect_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_water_single_bond

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_conect_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 1


    # insert multiple conect-bonds from record
    def test_place_multiple_conect_bonds(self, mock_protein_graph_methane, parameters_conect_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_methane

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_conect_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 4


    # correctly handle missing atoms for record
    def test_detect_missing_atoms(self, mock_protein_graph_methane_wrong_atom_numbers, parameters_conect_default):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_methane_wrong_atom_numbers

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        interaction_service.apply_interactions(protein_graph, parameters_conect_default)

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0


    # evaluate correct logging of (1) no bonds (2) inserted bonds
    def test_logging_no_conect(self,
            mock_protein_graph_methane_no_record,
            parameters_conect_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_methane_no_record

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_conect_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is None


    def test_logging_single_conect(self,
            mock_protein_graph_water_single_bond,
            parameters_conect_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_water_single_bond

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_conect_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 1
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is not None
        assert re.search(r"Added 1 .* CONECT", info_line)


    def test_logging_multiple_conects_in_record(self,
            mock_protein_graph_methane,
            parameters_conect_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_methane

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_conect_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 4
        info_line = next(filter(lambda x: x.startswith("INFO"), str(captured).split("\n")), None)
        assert info_line is not None
        assert re.search(r"Added 4 .* CONECT", info_line)


    def test_logging_warning_missing_atom_for_record(self,
            mock_protein_graph_methane_wrong_atom_numbers,
            parameters_conect_default,
            capsys
    ):
        interaction_service = InteractionServiceRegistry.COV.query_service(self.TESTED_SERVICE)
        protein_graph = mock_protein_graph_methane_wrong_atom_numbers

        before = len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value])

        with inject_pytest_logger():
            interaction_service.apply_interactions(protein_graph, parameters_conect_default)

        captured = capsys.readouterr()[1]

        assert len(protein_graph.graphs['pebble'].graph[GraphKey.COVALENT_EDGES.value]) - before == 0
        warning_line = next(filter(lambda x: x.startswith("WARNING"), str(captured).split("\n")), None)
        assert warning_line is not None
        assert re.search(r"Unable to find", warning_line)
        assert re.search(r"CONECT", warning_line)
