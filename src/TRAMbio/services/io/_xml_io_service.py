from contextlib import AbstractContextManager
from typing import Tuple, TextIO, List, Dict, Union, Optional

import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ParseError
import networkx as nx

from TRAMbio.services import ParameterRegistry, lock_registry
from TRAMbio.services.parameter import XtcParameter, PebbleGameParameter
from TRAMbio.util.constants.xml import XMLConstants, XML_SCHEME_COMPONENTS, components_namespace
from TRAMbio.services.io import IXmlIOService, IOServiceRegistry
from TRAMbio.util.structure_library.components import StructureRef
from loguru import logger


__all__ = []


class StateContext(AbstractContextManager):

    def __init__(self, xml_file: TextIO, structure_key: str):
        self._xml_file = xml_file
        self._structure_key = structure_key

    def __enter__(self):
        self._xml_file.write(
            f'{XMLConstants.XML_INDENT2.value}<{XMLConstants.STATE_TAG.value}'
            f' {XMLConstants.KEY_ATTRIBUTE_NAME.value}="{self._structure_key}">\n'
        )
        return self

    def __exit__(self, *exc_info):
        self._xml_file.write(f'{XMLConstants.XML_INDENT2.value}</{XMLConstants.STATE_TAG.value}>\n')


class XmlIOService(IXmlIOService):

    def __init__(self):
        super().__init__()
        try:
            from lxml.etree import XMLSchema, fromstring, XML
            self.__xml_schema_components = XMLSchema(fromstring(XML_SCHEME_COMPONENTS.encode("utf8")))
            self.__xml_schema_components_validator = lambda xml_string: \
                self.__xml_schema_components.assertValid(fromstring(xml_string.encode("utf8")))
        except ModuleNotFoundError:
            self.__xml_schema_components_validator = None

    @property
    def name(self):
        return "XmlIOService"

    def read(self, xml_path: str) -> Tuple[ET.Element, ET.Element]:
        # TODO: Check file
        with open(xml_path) as xml_f:
            xml_lines = xml_f.readlines()

        try:
            components_root = ET.fromstringlist(xml_lines)
        except ParseError:
            raise ValueError(f'Unable to parse {xml_path}')

        base_components = components_root.find(*components_namespace(XMLConstants.COMPONENTS_TAG))
        states = components_root.find(*components_namespace(XMLConstants.STATES_TAG))
        if base_components is None or states is None:
            raise KeyError(
                f'Components XML requires both <{XMLConstants.COMPONENTS_TAG.value}>'
                f' and <{XMLConstants.STATES_TAG.value}> elements.'
            )

        return base_components, states

    def read_graphml(self, graphml_path: str) -> Union[nx.Graph, nx.MultiGraph]:
        # TODO: Check file exists
        try:
            graph = nx.read_graphml(graphml_path)
        except nx.NetworkXError as e:
            exc = ValueError(f"Unable to load graph from file {graphml_path}.\nCause: {str(e)}")
            try:
                raise exc
            finally:
                exc.__context__ = None
        if graph.is_directed():
            logger.warning(
                f"Graph {graphml_path} is indicated as directed. Continuing with undirected version undirected graph."
            )
            graph = graph.to_undirected()
        return graph

    def validate_xml(self, xml_path: str) -> bool:
        """
        Validate document in ``xml_path``.

        Parameters
        ----------
        xml_path

        Raises
        ------
        lxml.etree.SyntaxError
            if the data in ``xml_path`` is no XML document.
        lxml.etree.DocumentInvalid
            if the document in ``xml_path`` does not comply with the components XML schema.

        Returns
        -------

        """
        if self.__xml_schema_components_validator is None:
            logger.warning(f"Module lxml not found. Unable to validate {xml_path}.")
            return True

        with open(xml_path) as xml_f:
            xml_string = '\n'.join(xml_f.readlines())
        self.__xml_schema_components_validator(xml_string)
        return True

    def state_context(
            self,
            xml_file: TextIO,
            structure_key: str
    ) -> AbstractContextManager:
        return StateContext(xml_file=xml_file, structure_key=structure_key)

    def write_temp_xml_fragment(
            self,
            xml_file: TextIO,
            halo: Optional[List[str]],
            sub_components: List[str]
    ) -> None:
        # opening tag
        xml_file.write(f'{XMLConstants.XML_INDENT3.value}<{XMLConstants.COMPONENT_TAG.value}>\n')

        # write component halo
        if halo is None or len(halo) == 0:
            xml_file.write(f"{XMLConstants.XML_INDENT4.value}<{XMLConstants.HALO_TAG.value}/>\n")
        else:
            xml_file.write(f"{XMLConstants.XML_INDENT4.value}<{XMLConstants.HALO_TAG.value}>\n")
            for node_id in halo:
                xml_file.write(
                    f"{XMLConstants.XML_INDENT5.value}<{XMLConstants.NODE_TAG.value}>"
                    f"{node_id}</{XMLConstants.NODE_TAG.value}>\n"
                )
            xml_file.write(f"{XMLConstants.XML_INDENT4.value}</{XMLConstants.HALO_TAG.value}>\n")

        # write subcomponents
        xml_file.write(f"{XMLConstants.XML_INDENT4.value}<{XMLConstants.COMPONENTS_TAG.value}>\n")
        for comp_id in sub_components:
            xml_file.write(
                f'{XMLConstants.SPECIAL_CHARACTER.value}{comp_id}\n')
        xml_file.write(f"{XMLConstants.XML_INDENT4.value}</{XMLConstants.COMPONENTS_TAG.value}>\n")

        # closing tag
        xml_file.write(f"{XMLConstants.XML_INDENT3.value}</{XMLConstants.COMPONENT_TAG.value}>\n")

    @lock_registry(kwargs_name='parameter_id')
    def convert_temp_to_xml(
            self,
            xml_path: str,
            temp_path: str,
            base_components: List[str],
            num_base_components: int,
            component_mapping: Dict[str, StructureRef],
            is_trajectory: bool,
            discarded_keys: Optional[List[str]] = None,
            parameter_id: str = ''
    ) -> None:
        if discarded_keys is None:
            discarded_keys = []

        stride = None
        if is_trajectory:
            parameter_registry = ParameterRegistry.get_parameter_set(parameter_id=parameter_id)
            stride = parameter_registry(XtcParameter.STRIDE.value)

        with open(xml_path, 'w') as xml_f:
            xml_f.write(
                f'<?xml version="1.0" encoding="UTF-8"?>\n'
                f'<{XMLConstants.GRAPH_TAG.value} {XMLConstants.XSD_COMPONENTS.value}>\n'
                f'{XMLConstants.XML_INDENT.value}<{XMLConstants.COMPONENTS_TAG.value}'
                f' {XMLConstants.SIZE_ATTRIBUTE_NAME.value}="{num_base_components}">\n'
            )

            # write base components
            for base_comp in base_components:
                xml_f.write(base_comp)

            attribute_string = f' {XMLConstants.KEY_ATTRIBUTE_NAME.value}="{stride}"' if stride is not None else ''
            xml_f.write(
                f'{XMLConstants.XML_INDENT.value}</{XMLConstants.COMPONENTS_TAG.value}>'
                f'\n{XMLConstants.XML_INDENT.value}<{XMLConstants.STATES_TAG.value}{attribute_string}>\n'
            )

            # read, (modify), and copy components from temp path
            within_state = False
            skip_state = False
            key_string = XMLConstants.KEY_ATTRIBUTE_NAME.value + "=\""

            with open(temp_path) as temp_xml_f:
                prefix_len = len(XMLConstants.SPECIAL_CHARACTER.value)
                for line in temp_xml_f:
                    if XMLConstants.STATE_TAG.value in line:
                        within_state = not within_state
                        if within_state:
                            index1 = line.find(key_string)
                            index2 = line.rfind("\"")
                            if index1 == -1 or index2 == -1:
                                raise KeyError
                            structure_key = line[index1 + len(key_string):index2]
                            skip_state = structure_key in discarded_keys

                    if skip_state:
                        continue
                    if len(line) == 0 or line[0] != XMLConstants.SPECIAL_CHARACTER.value:
                        xml_f.write(line)
                        continue

                    # read comp_id
                    comp_id = line[prefix_len:].strip()

                    # write output strings
                    for out_comp in component_mapping[comp_id]['out']:
                        xml_f.write(out_comp)

            xml_f.write(
                f'{XMLConstants.XML_INDENT.value}</{XMLConstants.STATES_TAG.value}>\n</{XMLConstants.GRAPH_TAG.value}>')

    @lock_registry(kwargs_name='parameter_id')
    def write_pebble_game_results(
            self,
            xml_out_path: str,
            category: str,
            components: List[List[str]],
            parameter_id: str = ''
    ) -> None:
        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id=parameter_id)
        param_k = parameter_registry(PebbleGameParameter.K.value)
        param_l = parameter_registry(PebbleGameParameter.L.value)

        with open(xml_out_path, 'w') as xml_f:
            xml_f.write(
                f'<?xml version="1.0" encoding="UTF-8"?>\n<{XMLConstants.GRAPH_TAG.value}'
                f' {XMLConstants.K_ATTRIBUTE_NAME.value}="{param_k}"'
                f' {XMLConstants.L_ATTRIBUTE_NAME.value}="{param_l}"'
                f' {XMLConstants.CATEGORY_ATTRIBUTE_NAME.value}="{category}"'
                f' {XMLConstants.XSD_PEBBLE.value}>\n'
                f'{XMLConstants.XML_INDENT.value}<{XMLConstants.COMPONENTS_TAG.value}'
                f' {XMLConstants.SIZE_ATTRIBUTE_NAME.value}="{len(components)}">\n'
            )
            for component in components:
                xml_f.write(
                    f'{XMLConstants.XML_INDENT2.value}<{XMLConstants.COMPONENT_TAG.value}'
                    f' {XMLConstants.SIZE_ATTRIBUTE_NAME.value}="{len(component)}">\n'
                )
                for node in component:
                    xml_f.write(
                        f'{XMLConstants.XML_INDENT3.value}<{XMLConstants.NODE_TAG.value}>'
                        f'{node}</{XMLConstants.NODE_TAG.value}>\n'
                    )
                xml_f.write(f"{XMLConstants.XML_INDENT2.value}</{XMLConstants.COMPONENT_TAG.value}>\n")
            xml_f.write(
                f"{XMLConstants.XML_INDENT.value}</{XMLConstants.COMPONENTS_TAG.value}>\n"
                f"</{XMLConstants.GRAPH_TAG.value}>"
            )


IOServiceRegistry.XML.register_service(XmlIOService())
