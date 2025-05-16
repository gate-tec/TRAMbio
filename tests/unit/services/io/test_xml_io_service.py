import os
import re
import sys
import time
from io import StringIO

import pytest
import textwrap

import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ParseError

from TRAMbio.services import ParameterRegistry, IOServiceRegistry
from TRAMbio.services.parameter import PebbleGameParameter, XtcParameter
from TRAMbio.util.constants.pebble_game import PebbleGameCategory
from TRAMbio.util.constants.xml import XMLConstants
from TRAMbio.util.structure_library.components import StructureRef
from tests.conftest import inject_pytest_logger


##################
# Parameters #####
##################

class TestParameters:

    TESTED_SERVICE = "XmlIOService"

    @pytest.fixture
    def parameters_xml_io_default(self):
        param_id = f"TEST_PARAMETER_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(param_id)

        registry.set_parameter(PebbleGameParameter.K.value, 2)
        registry.set_parameter(PebbleGameParameter.L.value, 3)
        registry.set_parameter(XtcParameter.STRIDE.value, 50)

        yield param_id


#################################
# No Mock Services required #####
#################################


##################
# Mock Paths #####
##################

class TestMockPaths:

    # paths for convert_temp_to_xml
    @pytest.fixture
    def temp_xml_path_C_0(self, tmp_path):
        f_temp = str(tmp_path / ".temp_test_components.xml")
        f_out = str(tmp_path / "test_components.xml")

        with open(f_temp, "w") as file:
            file.write("        " + textwrap.dedent("""
                    <state key="1">
                        <component>
                            <halo/>
                            <components>
            #C-0
                            </components>
                        </component>
                    </state>
            """).strip() + "\n")

        yield f_temp, f_out

        if os.path.exists(f_temp):
            os.remove(f_temp)
        if os.path.exists(f_out):
            os.remove(f_out)

    @pytest.fixture
    def temp_xml_path_C_1_C_2(self, tmp_path):
        f_temp = str(tmp_path / ".temp_test_components.xml")
        f_out = str(tmp_path / "test_components.xml")

        with open(f_temp, "w") as file:
            file.write("        " + textwrap.dedent("""
                        <state key="1">
                            <component>
                                <halo/>
                                <components>
                #C-1
                #C-2
                                </components>
                            </component>
                        </state>
                        <state key="2">
                            <component>
                                <halo/>
                                <components>
                #C-1
                                </components>
                            </component>
                            <component>
                                <halo/>
                                <components>
                #C-2
                                </components>
                            </component>
                        </state>
                """).strip() + "\n")

        yield f_temp, f_out

        if os.path.exists(f_temp):
            os.remove(f_temp)
        if os.path.exists(f_out):
            os.remove(f_out)

    # paths for write_pebble_game_results
    @pytest.fixture
    def out_path_pebble_game(self, tmp_path):
        f = str(tmp_path / "pebble_game.xml")

        yield f

        if os.path.exists(f):
            os.remove(f)


#################
# Mock Data #####
#################

class TestMockData:

    # mock data for state_context
    @pytest.fixture
    def mock_data_state_context_key_1(self):
        xml_buffer = StringIO()
        structure_key = "1"

        yield xml_buffer, structure_key

        xml_buffer.close()

    @pytest.fixture
    def mock_data_state_context_key_minus_1_5(self):
        xml_buffer = StringIO()
        structure_key = "-1.5"

        yield xml_buffer, structure_key

        xml_buffer.close()

    # mock data for write_temp_xml_fragment
    @pytest.fixture
    def mock_data_temp_fragment(self):
        xml_buffer = StringIO()
        halo = [
            "A0002-GLY:H",
            "A0002-GLY:CA"
        ]
        sub_components = [
            "C-0"
        ]

        yield xml_buffer, halo, sub_components

        xml_buffer.close()

    @pytest.fixture
    def mock_data_temp_fragment_no_halo(self):
        xml_buffer = StringIO()
        halo = None
        sub_components = [
            "C-0"
        ]

        yield xml_buffer, halo, sub_components

        xml_buffer.close()

    @pytest.fixture
    def mock_data_temp_fragment_multiple_sub_components(self):
        xml_buffer = StringIO()
        halo = [
            "A0002-GLY:H",
            "A0002-GLY:CA"
        ]
        sub_components = [
            "C-1",
            "C-2"
        ]

        yield xml_buffer, halo, sub_components

        xml_buffer.close()

    # mock data for convert_temp_to_xml
    @pytest.fixture
    def mock_data_base_components_C_0(self):
        base_components = [
            f"{' ' * 4 * 2}<component id=\"1\" size=\"5\">\n"
            f"{' ' * 4 * 3}<nodes>\n"
            f"{' ' * 4 * 4}<node>A0001-ARG:O</node>\n"
            f"{' ' * 4 * 4}<node>A0001-ARG:C</node>\n"
            f"{' ' * 4 * 4}<node>A0002-GLY:N</node>\n"
            f"{' ' * 4 * 3}</nodes>\n"
            f"{' ' * 4 * 3}<halo>\n"
            f"{' ' * 4 * 4}<node>A0002-GLY:CA</node>\n"
            f"{' ' * 4 * 4}<node>A0002-GLY:H</node>\n"
            f"{' ' * 4 * 3}</halo>\n"
            f"{' ' * 4 * 2}</component>\n"
        ]
        num_base_components = 1
        component_mapping = {
            "C-0": StructureRef(
                out=[
                    f"{' ' * 4 * 5}<component id=\"1\"/>\n"
                ], stack=None)
        }

        yield base_components, num_base_components, component_mapping

    @pytest.fixture
    def mock_data_base_components_C_1_C_2(self):
        base_components = [
            f"{' ' * 4 * 2}<component id=\"1\" size=\"5\">\n"
            f"{' ' * 4 * 3}<nodes>\n"
            f"{' ' * 4 * 4}<node>A0001-ARG:O</node>\n"
            f"{' ' * 4 * 4}<node>A0001-ARG:C</node>\n"
            f"{' ' * 4 * 4}<node>A0002-GLY:N</node>\n"
            f"{' ' * 4 * 3}</nodes>\n"
            f"{' ' * 4 * 3}<halo>\n"
            f"{' ' * 4 * 4}<node>A0002-GLY:CA</node>\n"
            f"{' ' * 4 * 4}<node>A0002-GLY:H</node>\n"
            f"{' ' * 4 * 3}</halo>\n"
            f"{' ' * 4 * 2}</component>\n",

            f"{' ' * 4 * 2}<component id=\"2\" size=\"1\">\n"
            f"{' ' * 4 * 3}<nodes>\n"
            f"{' ' * 4 * 4}<node>A0001-ARG:CA</node>\n"
            f"{' ' * 4 * 3}</nodes>\n"
            f"{' ' * 4 * 3}<halo/>\n"
            f"{' ' * 4 * 2}</component>\n",

            f"{' ' * 4 * 2}<component id=\"3\" size=\"1\">\n"
            f"{' ' * 4 * 3}<nodes>\n"
            f"{' ' * 4 * 4}<node>A0001-ARG:CA</node>\n"
            f"{' ' * 4 * 3}</nodes>\n"
            f"{' ' * 4 * 3}<halo/>\n"
            f"{' ' * 4 * 2}</component>\n"
        ]
        num_base_components = 1
        component_mapping = {
            "C-1": StructureRef(
                out=[
                    f"{' ' * 4 * 5}<component id=\"1\"/>\n"
                ], stack=None),
            "C-2": StructureRef(
                out=[
                    f"{' ' * 4 * 5}<component id=\"2\"/>\n"
                    f"{' ' * 4 * 5}<component id=\"3\"/>\n"
                ], stack=None)
        }

        yield base_components, num_base_components, component_mapping

    # mock data for write_pebble_game_results
    @pytest.fixture
    def mock_data_pebble_game_2_3(self):
        category = PebbleGameCategory.WELL_CONSTRAINED.value
        components = [
            ["n0", "n1", "n2"],
            ["n3", "n4", "n6", "n7"]
        ]

        yield category, components


#############
# Tests #####
#############

# test read
class TestRead(TestParameters):

    def test_read_xml_1(self, path_xml_components_1):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        base_components, states = xml_io_service.read(path_xml_components_1)

        # check components
        assert len(base_components) == 1
        assert base_components[0].get("id", None) == "1"  # verify id value

        # check states
        assert len(states) == 1
        first_comp_list = states[0][0].find("./tram:components", {"tram": "tram:components"})
        assert first_comp_list is not None
        assert len(first_comp_list) == 1
        assert first_comp_list[0].get("id", None) == "1"  # verify id reference

    def test_read_xml_2(self, path_xml_components_2):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        base_components, states = xml_io_service.read(path_xml_components_2)

        # check components
        assert len(base_components) == 3
        assert len({base_components[i].get("id", None) for i in range(3)}) == 3  # verify uniqueness of ids

        # check states
        assert len(states) == 2
        first_comp_list = states[0][0].find("./tram:components", {"tram": "tram:components"})
        assert first_comp_list is not None
        assert len(first_comp_list) == 3
        assert len({first_comp_list[i].get("id", None) for i in range(3)}) == 3  # verify id reference

    def test_read_broken_xml(self, path_xml_broken):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        with pytest.raises(KeyError, match="requires .* <components> .* <states>"):
            xml_io_service.read(path_xml_broken)



# test read_graphml
class TestReadGraphml(TestParameters):

    # - test parsing
    def test_read_simple_graph(self, path_graphml_simple):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        graph = xml_io_service.read_graphml(path_graphml_simple)

        assert graph.number_of_nodes() == 11
        assert len(graph.edges) == 12

    def test_read_broken_graph(self, path_graphml_broken):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        # simplify path check on Windows
        target_output = (os.path.basename(path_graphml_broken) if
                         sys.platform.startswith('win') else
                         path_graphml_broken)

        with pytest.raises(ValueError, match=target_output):
            xml_io_service.read_graphml(path_graphml_broken)

    # - test logging
    def test_read_directed_graph(self, path_graphml_directed, capsys):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        with inject_pytest_logger():
            graph = xml_io_service.read_graphml(path_graphml_directed)

        captured = capsys.readouterr()[1]

        assert graph.number_of_nodes() == 11
        assert len(graph.edges) == 12

        assert not graph.is_directed()
        warning_line = next(filter(lambda x: x.startswith("WARNING"), str(captured).split("\n")), None)
        assert warning_line is not None
        assert re.search(r"indicated as directed", warning_line)


# test validate_xml
class TestValidateXml(TestMockPaths, TestParameters):

    @pytest.mark.requires_lxml
    def test_correct_xml_1(self, path_xml_components_1):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        assert xml_io_service.validate_xml(path_xml_components_1)

    @pytest.mark.requires_lxml
    def test_correct_xml_2(self, path_xml_components_2):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        assert xml_io_service.validate_xml(path_xml_components_2)

    @pytest.mark.requires_lxml
    def test_broken_xml(self, path_xml_broken):
        from lxml.etree import DocumentInvalid

        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        with pytest.raises(DocumentInvalid, match=r"Missing child element.*\{tram:components\}components"):
            xml_io_service.validate_xml(path_xml_broken)


# test state_context
class TestStateContext(TestMockData, TestParameters):

    def test_base(self, mock_data_state_context_key_1):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        buffer, structure_key = mock_data_state_context_key_1

        context = xml_io_service.state_context(buffer, structure_key)

        with context:
            pass

        buffer.seek(0)
        lines = buffer.readlines()

        assert len(lines) == 2
        assert '"1"' in lines[0]

    def test_wrapping(self, mock_data_state_context_key_minus_1_5):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        buffer, structure_key = mock_data_state_context_key_minus_1_5
        dummy_text = "Test"

        context = xml_io_service.state_context(buffer, structure_key)

        with context:
            buffer.write("%s\n" % dummy_text)

        buffer.seek(0)
        lines = buffer.readlines()

        assert len(lines) == 3
        assert '"-1.5"' in lines[0]
        # test wrapping functionality
        assert lines[1].startswith(dummy_text)


# test write_temp_xml_fragment
class TestWriteTempXmlFragment(TestMockData, TestParameters):

    def test_base(self, mock_data_temp_fragment):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        buffer, halo, sub_components = mock_data_temp_fragment
        # opening component, opening halo => 2
        # closing halo, opening sub_components => 2
        # closing sub_components, closing component => 2
        expected_lines = 2 + len(halo) + 2 + len(sub_components) + 2

        xml_io_service.write_temp_xml_fragment(buffer, halo, sub_components)

        buffer.seek(0)

        lines = buffer.readlines()
        assert len(lines) == expected_lines
        component_lines = lines[2 + len(halo) + 2:-2]
        for a, b in zip(component_lines, sub_components):
            assert a.startswith(XMLConstants.SPECIAL_CHARACTER.value)
            assert a[len(XMLConstants.SPECIAL_CHARACTER.value):].startswith(b)

    def test_no_halo(self, mock_data_temp_fragment_no_halo):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        buffer, halo, sub_components = mock_data_temp_fragment_no_halo
        # opening component, single halo tag, opening sub_components => 3
        # closing sub_components, closing component => 2
        expected_lines = 3 + len(sub_components) + 2

        xml_io_service.write_temp_xml_fragment(buffer, halo, sub_components)

        buffer.seek(0)

        lines = buffer.readlines()
        assert len(lines) == expected_lines

        assert "/>" in lines[1]  # halo should be single tag

        component_lines = lines[2 + 1 + 2:-2]
        for a, b in zip(component_lines, sub_components):
            assert a.startswith(XMLConstants.SPECIAL_CHARACTER.value)
            assert a[len(XMLConstants.SPECIAL_CHARACTER.value):].startswith(b)

    def test_multiple_sub_components(self, mock_data_temp_fragment_multiple_sub_components):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        buffer, halo, sub_components = mock_data_temp_fragment_multiple_sub_components
        # opening component, opening halo => 2
        # closing halo, opening sub_components => 2
        # closing sub_components, closing component => 2
        expected_lines = 2 + len(halo) + 2 + len(sub_components) + 2

        xml_io_service.write_temp_xml_fragment(buffer, halo, sub_components)

        buffer.seek(0)

        lines = buffer.readlines()
        assert len(lines) == expected_lines
        component_lines = lines[2 + len(halo) + 2:-2]
        for a, b in zip(component_lines, sub_components):
            assert a.startswith(XMLConstants.SPECIAL_CHARACTER.value)
            assert a[len(XMLConstants.SPECIAL_CHARACTER.value):].startswith(b)


# test convert_temp_to_xml
class TestConvertTempToXml(TestMockPaths, TestMockData, TestParameters):

    def test_convert_C_0(self, temp_xml_path_C_0, mock_data_base_components_C_0, parameters_xml_io_default):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        temp_path, out_path = temp_xml_path_C_0
        base_components, num_base_components, component_mapping = mock_data_base_components_C_0

        xml_io_service.convert_temp_to_xml(
            temp_path=temp_path,
            xml_path=out_path,
            base_components=base_components,
            num_base_components=num_base_components,
            component_mapping=component_mapping,
            is_trajectory=False,
            parameter_id=parameters_xml_io_default
        )

        # verify output file exists
        assert os.path.exists(out_path)
        with open(out_path) as file:
            xml_lines = file.readlines()

        # verify XML structure
        try:
            components_root = ET.fromstringlist(xml_lines)
        except ParseError:
            pytest.fail(reason="Unable to parse constructed XML file.")

        # verify insertion of base components
        parsed_base_components = components_root.find("{tram:components}components")  # noqa components_root is always set
        assert parsed_base_components is not None
        assert len(parsed_base_components) == 1

        # verify copying of states from temp file
        parsed_states = components_root.find("{tram:components}states")
        assert parsed_states is not None
        assert len(parsed_states) == 1

        # verify replacement of C-0 with component reference
        injected_component = parsed_states.find(".//{tram:components}components/{tram:components}component")
        assert injected_component is not None
        assert injected_component.get("id", None) == "1"


    def test_convert_C_1_C_2(self, temp_xml_path_C_1_C_2, mock_data_base_components_C_1_C_2, parameters_xml_io_default):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        temp_path, out_path = temp_xml_path_C_1_C_2
        base_components, num_base_components, component_mapping = mock_data_base_components_C_1_C_2

        xml_io_service.convert_temp_to_xml(
            temp_path=temp_path,
            xml_path=out_path,
            base_components=base_components,
            num_base_components=num_base_components,
            component_mapping=component_mapping,
            is_trajectory=True,
            parameter_id=parameters_xml_io_default
        )

        # verify output file exists
        assert os.path.exists(out_path)
        with open(out_path) as file:
            xml_lines = file.readlines()

        # verify XML structure
        try:
            components_root = ET.fromstringlist(xml_lines)
        except ParseError:
            pytest.fail(reason="Unable to parse constructed XML file.")

        # verify insertion of base components
        parsed_base_components = components_root.find("{tram:components}components")  # noqa components_root is always set
        assert parsed_base_components is not None
        assert len(parsed_base_components) == 3

        # verify copying of states from temp file
        parsed_states = components_root.find("{tram:components}states")
        assert parsed_states is not None
        assert len(parsed_states) == 2

        # check stride value for trajectory
        assert parsed_states.get("key", None) == "50"

        # check component list of first component in first state
        comp_1 = parsed_states.find(
            "./tram:state[1]/tram:component/tram:components",
            {"tram": "tram:components"}
        )
        assert comp_1 is not None
        # C-1 should expand to 1 component
        # C-2 should expand to 2 component
        assert len(comp_1) == 3

        # check second state
        state_2 = parsed_states.find("./tram:state[2]", {"tram": "tram:components"})
        assert state_2 is not None
        assert state_2.get("key", None) == "2"
        assert len(state_2) == 2

        # check component list of second component in second state
        comp_2 = state_2.find("./tram:component[2]/tram:components", {"tram": "tram:components"})
        assert comp_2 is not None
        # C-2 should expand to 2 component
        assert len(comp_2) == 2

    def test_convert_C_1_C_2_with_discard(self, temp_xml_path_C_1_C_2, mock_data_base_components_C_1_C_2, parameters_xml_io_default):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        temp_path, out_path = temp_xml_path_C_1_C_2
        base_components, num_base_components, component_mapping = mock_data_base_components_C_1_C_2

        xml_io_service.convert_temp_to_xml(
            temp_path=temp_path,
            xml_path=out_path,
            base_components=base_components,
            num_base_components=num_base_components,
            component_mapping=component_mapping,
            is_trajectory=False,
            discarded_keys=["2"],
            parameter_id=parameters_xml_io_default
        )

        # verify output file exists
        assert os.path.exists(out_path)
        with open(out_path) as file:
            xml_lines = file.readlines()

        # verify XML structure
        try:
            components_root = ET.fromstringlist(xml_lines)
        except ParseError:
            pytest.fail(reason="Unable to parse constructed XML file.")

        # verify insertion of base components
        parsed_base_components = components_root.find("{tram:components}components")  # noqa components_root is always set
        assert parsed_base_components is not None
        assert len(parsed_base_components) == 3

        # verify copying of states from temp file
        parsed_states = components_root.find("{tram:components}states")
        assert parsed_states is not None
        assert len(parsed_states) == 1  # discarded state "2"

        # check component list of first component in first state
        comp_1 = parsed_states.find(
            "./tram:state[1]/tram:component/tram:components",
            {"tram": "tram:components"}
        )
        assert comp_1 is not None
        # C-1 should expand to 1 component
        # C-2 should expand to 2 component
        assert len(comp_1) == 3


# test write_pebble_game_results
class TestWritePebbleGameResults(TestMockPaths, TestMockData, TestParameters):

    def test_pebble_game_results_2_3(self, out_path_pebble_game, mock_data_pebble_game_2_3, parameters_xml_io_default):
        xml_io_service = IOServiceRegistry.XML.query_service(self.TESTED_SERVICE)

        category, components = mock_data_pebble_game_2_3

        xml_io_service.write_pebble_game_results(out_path_pebble_game, category, components, parameters_xml_io_default)

        # verify output file exists
        assert os.path.exists(out_path_pebble_game)
        with open(out_path_pebble_game) as file:
            lines = file.readlines()

        # verify XML structure
        try:
            components_root = ET.fromstringlist(lines)
        except ParseError:
            pytest.fail(reason="Unable to parse constructed XML file.")

        # verify graph attributes
        attrib_k = components_root.get("k", None)  # noqa components_root is always set
        assert attrib_k == "2"
        attrib_l = components_root.get("l", None)
        assert attrib_l == "3"

        # check <components>
        assert len(components_root) == 1
        components = components_root.find("./tram:components", {"tram": "tram:pebble"})
        assert components is not None
        assert len(components) == 2
