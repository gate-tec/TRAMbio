import os

import pytest

from TRAMbio.services import WorkflowServiceRegistry
from TRAMbio.util.structure_library.generator import CustomGenerator
from TRAMbio.util.structure_library.components import PebbleGameResult, IntermediateComponents


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "PdbWorkflowService"

    @staticmethod
    def assert_content_order(lines, content_list):
        iterator = iter(content_list)
        item = None
        for line in lines:
            if item is None:
                try:
                    item = next(iterator)
                except StopIteration:
                    return

            if item in line:
                item = None

        if item is not None:
            pytest.fail(f"Missing item '{item}' in data.")


##################
# Test Paths #####
##################

class TestPaths:

    @pytest.fixture
    def temp_out_path(self, tmp_path):
        f = str(tmp_path / "temp.xml")

        yield f

        if os.path.exists(f):
            os.remove(f)


######################
# Mock Generator #####
######################

class TestMockGenerator:

    @pytest.fixture
    def mock_generator_identical_components(self):
        def generator():
            yield "0", [
                PebbleGameResult(
                    size=3, halo_size=0,
                    nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N"], halo=[]
                )
            ]

            yield "1", [
                PebbleGameResult(
                    size=3, halo_size=0,
                    nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N"], halo=[]
                )
            ]

        yield CustomGenerator(generator())

    @pytest.fixture
    def mock_generator_identical_components_with_halo(self):
        def generator():
            yield "0", [
                PebbleGameResult(
                    size=3, halo_size=1,
                    nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N"], halo=["A0002-GLY:H"]
                )
            ]

            yield "1", [
                PebbleGameResult(
                    size=3, halo_size=1,
                    nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N"], halo=["A0001-GLY:CA"]
                )
            ]

        yield CustomGenerator(generator())

    @pytest.fixture
    def mock_generator_disjunct_components(self):
        def generator():
            yield "0", [
                PebbleGameResult(
                    size=3, halo_size=0,
                    nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N"], halo=[]
                )
            ]

            yield "1", [
                PebbleGameResult(
                    size=4, halo_size=0,
                    nodes=["A0003-GLY:C", "A0003-GLY:O", "A0004-GLY:N", "A0004-GLY:CA"], halo=[]
                )
            ]

        yield CustomGenerator(generator())

    @pytest.fixture
    def mock_generator_subset_components(self):
        def generator():
            yield "0", [
                PebbleGameResult(
                    size=3, halo_size=0,
                    nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N"], halo=[]
                )
            ]

            yield "1", [
                PebbleGameResult(
                    size=4, halo_size=0,
                    nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N", "A0002-GLY:CA"], halo=[]
                )
            ]

        yield CustomGenerator(generator())

    @pytest.fixture
    def mock_generator_supset_components(self):
        def generator():
            yield "0", [
                PebbleGameResult(
                    size=4, halo_size=0,
                    nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N", "A0002-GLY:CA"], halo=[]
                )
            ]

            yield "1", [
                PebbleGameResult(
                    size=3, halo_size=0,
                    nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N"], halo=[]
                )
            ]

        yield CustomGenerator(generator())

    @pytest.fixture
    def mock_generator_overlapping_components(self):
        def generator():
            yield "0", [
                PebbleGameResult(
                    size=4, halo_size=0,
                    nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N", "A0002-GLY:CA"], halo=[]
                )
            ]

            yield "1", [
                PebbleGameResult(
                    size=3, halo_size=0,
                    nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N", "A0002-GLY:H"], halo=[]
                )
            ]

        yield CustomGenerator(generator())


#################
# Mock Data #####
#################

class TestMockData:

    @pytest.fixture
    def mock_data_single_component(self):
        comp_dict = {
            "C-0": IntermediateComponents(
                size=3, nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N"], components=None
            )
        }
        hydrogen_dict = {}

        yield comp_dict, hydrogen_dict

    @pytest.fixture
    def mock_data_single_component_with_halo(self):
        comp_dict = {
            "C-0": IntermediateComponents(
                size=3, nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N"], components=None
            )
        }
        hydrogen_dict = {"A0002-GLY:N": ["A0002-GLY:H"]}

        yield comp_dict, hydrogen_dict

    @pytest.fixture
    def mock_data_overlapping_components(self):
        comp_dict = {
            "C-0": IntermediateComponents(
                size=4, nodes=None, components=["C-1", "C-2"]
            ),
            "C-1": IntermediateComponents(
                size=3, nodes=["A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N"], components=None
            ),
            "C-2": IntermediateComponents(
                size=1, nodes=["A0001-GLY:CA"], components=None
            )
        }
        hydrogen_dict = {}

        yield comp_dict, hydrogen_dict


#############
# Tests #####
#############

# tests for calculate_components_from_generator
class TestCalculateComponentsFromGenerator(TestParameters, TestMockGenerator, TestPaths):

    def test_identical_components(self, temp_out_path, mock_generator_identical_components):
        workflow_service = WorkflowServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        comp_dict, discard_list = workflow_service.calculate_components_from_generator(
            mock_generator_identical_components, temp_out_path
        )

        assert len(comp_dict) == 1
        assert len(discard_list) == 1
        assert "1" in discard_list

        assert "C-0" in comp_dict.keys()
        assert comp_dict["C-0"]["size"] == 3

        with open(temp_out_path, "r") as file:
            self.assert_content_order(file.readlines(), [
                "key=\"0\"", "C-0", "key=\"1\"", "C-0"
            ])

    def test_identical_components_with_halo(self, temp_out_path, mock_generator_identical_components_with_halo):
        workflow_service = WorkflowServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        comp_dict, discard_list = workflow_service.calculate_components_from_generator(
            mock_generator_identical_components_with_halo, temp_out_path
        )

        assert len(comp_dict) == 1
        assert len(discard_list) == 0

        assert "C-0" in comp_dict.keys()
        assert comp_dict["C-0"]["size"] == 3

        with open(temp_out_path, "r") as file:
            self.assert_content_order(file.readlines(), [
                "key=\"0\"", "A0002-GLY:H", "C-0", "key=\"1\"", "A0001-GLY:CA", "C-0"
            ])

    def test_disjunct_components(self, temp_out_path, mock_generator_disjunct_components):
        workflow_service = WorkflowServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        comp_dict, discard_list = workflow_service.calculate_components_from_generator(
            mock_generator_disjunct_components, temp_out_path
        )

        assert len(comp_dict) == 2
        assert len(discard_list) == 0

        assert "C-0" in comp_dict.keys()
        assert comp_dict["C-0"]["size"] == 3

        assert "C-1" in comp_dict.keys()
        assert comp_dict["C-1"]["size"] == 4

        with open(temp_out_path, "r") as file:
            self.assert_content_order(file.readlines(), [
                "key=\"0\"", "C-0", "key=\"1\"", "C-1"
            ])

    def test_subset_components(self, temp_out_path, mock_generator_subset_components):
        workflow_service = WorkflowServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        comp_dict, discard_list = workflow_service.calculate_components_from_generator(
            mock_generator_subset_components, temp_out_path
        )

        assert len(comp_dict) == 2
        assert len(discard_list) == 0

        assert "C-0" in comp_dict.keys()
        assert comp_dict["C-0"]["size"] == 3

        assert "C-1" in comp_dict.keys()
        assert comp_dict["C-1"]["size"] == 1

        with open(temp_out_path, "r") as file:
            self.assert_content_order(file.readlines(), [
                "key=\"0\"", "C-0", "key=\"1\"", "C-1", "C-0"
            ])

    def test_supset_components(self, temp_out_path, mock_generator_supset_components):
        workflow_service = WorkflowServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        comp_dict, discard_list = workflow_service.calculate_components_from_generator(
            mock_generator_supset_components, temp_out_path
        )

        assert len(comp_dict) == 3
        assert len(discard_list) == 0

        assert "C-0" in comp_dict.keys()
        assert comp_dict["C-0"]["size"] == 4
        assert "C-1" in comp_dict["C-0"]["components"]
        assert "C-2" in comp_dict["C-0"]["components"]

        assert "C-1" in comp_dict.keys()
        assert comp_dict["C-1"]["size"] == 3

        assert "C-2" in comp_dict.keys()
        assert comp_dict["C-2"]["size"] == 1

        with open(temp_out_path, "r") as file:
            self.assert_content_order(file.readlines(), [
                "key=\"0\"", "C-0", "key=\"1\"", "C-1"
            ])

    def test_overlapping_components(self, temp_out_path, mock_generator_overlapping_components):
        workflow_service = WorkflowServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        comp_dict, discard_list = workflow_service.calculate_components_from_generator(
            mock_generator_overlapping_components, temp_out_path
        )

        assert len(comp_dict) == 4
        assert len(discard_list) == 0

        assert "C-0" in comp_dict.keys()
        assert comp_dict["C-0"]["size"] == 4
        assert "C-2" in comp_dict["C-0"]["components"]
        assert "C-3" in comp_dict["C-0"]["components"]

        assert "C-1" in comp_dict.keys()
        assert comp_dict["C-1"]["size"] == 1

        assert "C-2" in comp_dict.keys()
        assert comp_dict["C-2"]["size"] == 3

        assert "C-3" in comp_dict.keys()
        assert comp_dict["C-3"]["size"] == 1

        with open(temp_out_path, "r") as file:
            self.assert_content_order(file.readlines(), [
                "key=\"0\"", "C-0", "key=\"1\"", "C-1", "C-2"
            ])


# tests for convert_component_archive_to_mapping
class TestConvertComponentArchiveToMapping(TestParameters, TestMockData):

    def test_single_component(self, mock_data_single_component):
        workflow_service = WorkflowServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        comp_dict, hydrogen_mapping = mock_data_single_component

        base_components, num_base_components, component_mapping = workflow_service.convert_component_archive_to_mapping(
            comp_dict, hydrogen_mapping
        )

        assert num_base_components == 1
        self.assert_content_order(base_components[0].split("\n"), [
            "structure=\"peptide unit\"", "<nodes>", "A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N", "<halo/>"
        ])

        assert len(component_mapping) == 1
        assert "C-0" in component_mapping.keys()
        assert len(component_mapping["C-0"]["out"]) == 1  # single component entry
        assert component_mapping["C-0"]["stack"] is None  # empty stack

    def test_single_component_with_halo(self, mock_data_single_component_with_halo):
        workflow_service = WorkflowServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        comp_dict, hydrogen_mapping = mock_data_single_component_with_halo

        base_components, num_base_components, component_mapping = workflow_service.convert_component_archive_to_mapping(
            comp_dict, hydrogen_mapping
        )

        assert num_base_components == 1
        self.assert_content_order(base_components[0].split("\n"), [
            "structure=\"peptide unit\"", "<nodes>", "A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N", "<halo>", "A0002-GLY:H"
        ])

        assert len(component_mapping) == 1
        assert "C-0" in component_mapping.keys()
        assert len(component_mapping["C-0"]["out"]) == 1  # single component entry
        assert component_mapping["C-0"]["stack"] is None  # empty stack

    def test_overlapping_components(self, mock_data_overlapping_components):
        workflow_service = WorkflowServiceRegistry.PDB.query_service(self.TESTED_SERVICE)

        comp_dict, hydrogen_mapping = mock_data_overlapping_components

        base_components, num_base_components, component_mapping = workflow_service.convert_component_archive_to_mapping(
            comp_dict, hydrogen_mapping
        )

        assert num_base_components == 2
        self.assert_content_order(base_components[0].split("\n"), [
            "structure=\"peptide unit\"", "<nodes>", "A0001-GLY:C", "A0001-GLY:O", "A0002-GLY:N", "<halo/>",
        ])
        self.assert_content_order(base_components[1].split("\n"), [
            "<nodes>", "A0001-GLY:CA", "<halo/>"
        ])

        assert len(component_mapping) == 3
        for comp_id, expected_lines in [("C-0", 2), ("C-1", 1), ("C-2", 1)]:
            assert comp_id in component_mapping.keys()
            assert len(component_mapping[comp_id]["out"]) == expected_lines
            assert component_mapping[comp_id]["stack"] is None  # empty stack

        # verify that component mapping of C-0 is a union of C-1 and C-2
        list_c0 = list(sorted(component_mapping["C-0"]["out"]))
        list_c1_c2 = list(sorted(component_mapping["C-1"]["out"] + component_mapping["C-2"]["out"]))
        assert len(list_c0) == len(list_c1_c2)
        assert all(map(lambda x: x[0] == x[1], zip(list_c0, list_c1_c2)))
