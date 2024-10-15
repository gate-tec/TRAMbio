import textwrap

import pytest

import xml.etree.ElementTree as ET
from TRAMbio.services import StructureServiceRegistry


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "XmlStructureService"


#################################
# No Mock Services required #####
#################################


#################
# Mock Data #####
#################

class TestMockData:

    @pytest.fixture
    def mock_data_base_components_1(self):
        base_component_xml = textwrap.dedent("""
        <?xml version="1.0" encoding="UTF-8"?>
        <components xmlns="tram:components" size="1">
            <component id="1" size="5">
                <nodes>
                    <node>A0001-ARG:O</node>
                    <node>A0001-ARG:C</node>
                    <node>A0002-GLY:N</node>
                </nodes>
                <halo>
                    <node>A0002-GLY:CA</node>
                    <node>A0002-GLY:H</node>
                </halo>
            </component>
        </components>
        """).strip()

        target_dict = {
            "1": {"A0001-ARG:O", "A0001-ARG:C", "A0002-GLY:N", "A0002-GLY:CA", "A0002-GLY:H"},
        }

        yield ET.fromstringlist(base_component_xml.split("\n")), target_dict

    @pytest.fixture
    def mock_data_base_components_1_2_3(self):
        base_component_xml = textwrap.dedent("""
            <?xml version="1.0" encoding="UTF-8"?>
            <components xmlns="tram:components" size="3">
                <component id="1" size="5">
                    <nodes>
                        <node>A0001-ARG:O</node>
                        <node>A0001-ARG:C</node>
                        <node>A0002-GLY:N</node>
                    </nodes>
                    <halo>
                        <node>A0002-GLY:CA</node>
                        <node>A0002-GLY:H</node>
                    </halo>
                </component>
                <component id="2" size="1">
                    <nodes>
                        <node>A0003-ARG:CA</node>
                    </nodes>
                    <halo/>
                </component>
                <component id="3" size="1">
                    <nodes>
                        <node>A0004-PHE:CA</node>
                    </nodes>
                    <halo/>
                </component>
            </components>
            """).strip()

        target_dict = {
            "1": {"A0001-ARG:O", "A0001-ARG:C", "A0002-GLY:N", "A0002-GLY:CA", "A0002-GLY:H"},
            "2": {"A0003-ARG:CA"},
            "3": {"A0004-PHE:CA"}
        }

        yield ET.fromstringlist(base_component_xml.split("\n")), target_dict

    @pytest.fixture
    def mock_data_states_1(self):
        states_xml = textwrap.dedent("""
            <?xml version="1.0" encoding="UTF-8"?>
            <states xmlns="tram:components">
                <state key="1">
                    <component>
                        <halo/>
                        <components>
                            <component id="1"/>
                        </components>
                    </component>
                </state>
            </states>
            """).strip()

        yield ET.fromstringlist(states_xml.split("\n"))

    @pytest.fixture
    def mock_data_states_1_2(self):
        states_xml = textwrap.dedent("""
            <?xml version="1.0" encoding="UTF-8"?>
            <states xmlns="tram:components">
                <state key="1">
                    <component>
                        <halo/>
                        <components>
                            <component id="1"/>
                            <component id="2"/>
                            <component id="3"/>
                        </components>
                    </component>
                </state>
                <state key="2">
                    <component>
                        <halo/>
                        <components>
                            <component id="1"/>
                        </components>
                    </component>
                    <component>
                        <halo/>
                        <components>
                            <component id="2"/>
                            <component id="3"/>
                        </components>
                    </component>
                </state>
            </states>
            """).strip()

        yield ET.fromstringlist(states_xml.split("\n"))

    # mock data for consistent_color_neighbor_states
    @pytest.fixture
    def mock_data_neighbor_states_3_3(self):
        state_1 = [
            {"1", "2"},
            {"3", "4"},
            {"5", "6"}
        ]
        state_2 = [
            {"1", "2"},
            {"5", "6"},
            {"3", "4"},
        ]

        base_lengths = {
            "1": 5,
            "2": 1,
            "3": 2,
            "4": 3,
            "5": 1,
            "6": 3
        }

        color_map = {
            0: 0,
            1: 3,
            2: 4
        }
        next_color = 5

        yield state_1, state_2, base_lengths, color_map, next_color

    @pytest.fixture
    def mock_data_neighbor_states_3_4(self):
        state_1 = [
            {"1", "2"},
            {"3", "4"},
            {"5", "6"}
        ]
        state_2 = [
            {"2"},
            {"1"},
            {"3", "4"},
            {"6"}
        ]

        base_lengths = {
            "1": 5,
            "2": 1,
            "3": 2,
            "4": 3,
            "5": 1,
            "6": 3
        }

        color_map = {
            0: 0,
            1: 3,
            2: 4
        }
        next_color = 5

        yield state_1, state_2, base_lengths, color_map, next_color


#############
# Tests #####
#############

# test resolve_base_mapping
class TestResolveBaseMapping(TestMockData, TestParameters):

    def test_resolve_base_comp_1(self, mock_data_base_components_1):
        structure_service = StructureServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        base_xml, target_dict = mock_data_base_components_1

        base_components, base_lengths = structure_service.resolve_base_mapping(base_xml)

        assert len(base_components) == len(target_dict)
        for k, v in target_dict.items():
            assert k in base_components.keys()
            assert len(base_components[k]) == len(v)
            assert base_lengths[k] == len(v)
            assert base_components[k] == v

    def test_resolve_base_comp_1_2_3(self, mock_data_base_components_1_2_3):
        structure_service = StructureServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        base_xml, target_dict = mock_data_base_components_1_2_3

        base_components, base_lengths = structure_service.resolve_base_mapping(base_xml)

        assert len(base_components) == len(target_dict)
        for k, v in target_dict.items():
            assert k in base_components.keys()
            assert len(base_components[k]) == len(v)
            assert base_lengths[k] == len(v)
            assert base_components[k] == v


# test create_list_from_state
class TestCreateListFromState(TestMockData, TestParameters):

    def test_states_1_state_1(self, mock_data_states_1, mock_data_base_components_1):
        structure_service = StructureServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        _, node_map = mock_data_base_components_1

        comp_list, node_list = structure_service.create_list_from_state(mock_data_states_1[0], node_map)

        assert len(comp_list) == 1
        assert node_list[0] == node_map["1"]

    def test_states_2_state_1(self, mock_data_states_1_2, mock_data_base_components_1_2_3):
        structure_service = StructureServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        _, node_map = mock_data_base_components_1_2_3

        comp_list, node_list = structure_service.create_list_from_state(mock_data_states_1_2[0], node_map)

        assert len(comp_list) == 1
        assert node_list[0] == node_map["1"] | node_map["2"] | node_map["3"]

    def test_states_2_state_2(self, mock_data_states_1_2, mock_data_base_components_1_2_3):
        structure_service = StructureServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        _, node_map = mock_data_base_components_1_2_3

        comp_list, node_list = structure_service.create_list_from_state(mock_data_states_1_2[1], node_map)

        assert len(comp_list) == 2
        assert node_list[0] == node_map["1"]
        assert node_list[1] == node_map["2"] | node_map["3"]


# test consistent_color_neighbor_states
class TestConsistentColorNeighborStates(TestMockData, TestParameters):

    def test_color_states_3_3(self, mock_data_neighbor_states_3_3):
        structure_service = StructureServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        state_1, state_2, base_lengths, color_map, next_color = mock_data_neighbor_states_3_3

        new_color_map, new_next_color = structure_service.consistent_color_neighbor_states(
            state_1, state_2, base_lengths, color_map, next_color
        )

        assert len(new_color_map) == len(state_2)
        assert new_color_map[0] == color_map[0]
        assert new_color_map[1] == color_map[2]
        assert new_color_map[2] == color_map[1]

    def test_color_states_3_4(self, mock_data_neighbor_states_3_4):
        structure_service = StructureServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        state_1, state_2, base_lengths, color_map, next_color = mock_data_neighbor_states_3_4

        new_color_map, new_next_color = structure_service.consistent_color_neighbor_states(
            state_1, state_2, base_lengths, color_map, next_color
        )

        assert len(new_color_map) == len(state_2)
        assert new_color_map[0] == next_color
        assert new_color_map[1] == color_map[0]
        assert new_color_map[2] == color_map[1]
        assert new_color_map[3] == color_map[2]


# test consistent_color_components
class TestConsistentColorComponents(TestMockData, TestParameters):

    def test_states_1(self, mock_data_states_1, mock_data_base_components_1):
        structure_service = StructureServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        base_components, _ = mock_data_base_components_1

        generator = structure_service.consistent_color_components(base_components, mock_data_states_1)

        assert generator is not None

        iterator = iter(generator)
        try:
            key, node_sets, color_map = next(iterator)
        except StopIteration:
            pytest.fail("Generator terminated before expected yield.")

        assert key == "1"  # noqa; state key
        assert len(node_sets) == 1  # noqa; only one component in state
        assert len(node_sets[0]) == 5  # noqa; five nodes in first component
        assert color_map[0] == 1  # noqa; first component has color 1

        try:
            next(iterator)
        except StopIteration:
            pass
        else:
            pytest.fail("Unexpected additional yield in generator.")

        assert generator.stop() == 2

    def test_states_2(self, mock_data_states_1_2, mock_data_base_components_1_2_3):
        structure_service = StructureServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        base_components, _ = mock_data_base_components_1_2_3

        generator = structure_service.consistent_color_components(base_components, mock_data_states_1_2)

        assert generator is not None

        iterator = iter(generator)

        # Validate first state
        try:
            key, node_sets, color_map = next(iterator)
        except StopIteration:
            pytest.fail("Generator terminated before expected yield.")

        assert key == "1"  # noqa; state key
        assert len(node_sets) == 1  # noqa; only one component in state
        assert len(node_sets[0]) == 7  # noqa; seven nodes in first component
        assert color_map[0] == 1  # noqa; first component has color 1

        # Validate second state
        try:
            key, node_sets, color_map = next(iterator)
        except StopIteration:
            pytest.fail("Generator terminated before expected yield.")

        assert key == "2"  # state key
        assert len(node_sets) == 2  # two components in state
        assert len(node_sets[0]) == 5  # five nodes in first component
        assert len(node_sets[1]) == 2  # two nodes in second component
        assert color_map[0] == 1  # first component has color 1
        assert color_map[1] == 2  # second component has color 2

        try:
            next(iterator)
        except StopIteration:
            pass
        else:
            pytest.fail("Unexpected additional yield in generator.")

        assert generator.stop() == 3
