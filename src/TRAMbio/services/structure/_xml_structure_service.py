from typing import Dict, Set, Tuple, List
import xml.etree.ElementTree as ET

from TRAMbio.services.structure import IXmlStructureService, StructureServiceRegistry
from TRAMbio.util.constants.xml import XMLConstants, components_namespace
from TRAMbio.util.structure_library.generator import as_custom_generator, CustomGenerator


__all__ = []


class XmlStructureService(IXmlStructureService):

    @property
    def name(self):
        return "XmlStructureService"

    def resolve_base_mapping(
            self,
            base_components: ET.Element
    ) -> Tuple[Dict[str, Set[str]], Dict[str, int]]:
        base_mapping: Dict[str, Set[str]] = {}
        base_lengths: Dict[str, int] = {}

        for component in base_components:
            comp_id = component.get('id')
            if comp_id is None:
                # sanity check; should never reach
                continue

            nodes = set()

            halo_list = component.find(*components_namespace(XMLConstants.HALO_TAG))
            nodes_list = component.find(*components_namespace(XMLConstants.NODES_TAG))

            if halo_list is not None:
                for node in halo_list:
                    nodes.add(node.text)
            if nodes_list is not None:
                for node in nodes_list:
                    nodes.add(node.text)

            node_len = len(nodes)
            if node_len > 0:
                base_mapping[comp_id] = nodes
                base_lengths[comp_id] = node_len

        return base_mapping, base_lengths

    def create_list_from_state(
            self,
            xml_state_element: ET.Element,
            node_map: Dict[str, Set]
    ) -> Tuple[List[Set[str]], List[Set[str]]]:
        component_list = []

        for component in xml_state_element:
            halo_nodes = component.find(*components_namespace(XMLConstants.HALO_TAG))
            sub_components = component.find(*components_namespace(XMLConstants.COMPONENTS_TAG))
            if sub_components is None:
                # sanity check; should never reach
                continue

            comp = set()
            comp_nodes = set()
            if halo_nodes is not None:
                for node in halo_nodes:
                    comp_nodes.add(node.text)

            for sub_component in sub_components:
                comp_id = sub_component.get('id')
                if comp_id is not None:
                    comp.add(comp_id)
                    # resolve nodes
                    comp_nodes = comp_nodes | node_map[comp_id]

            component_list.append((comp, comp_nodes, len(comp_nodes)))

        # list from largest to smallest component
        component_list, nodes_list = zip(*list(
            (comp, comp_nodes) for comp, comp_nodes, _ in sorted(component_list, reverse=True, key=lambda x: x[2])))
        return list(component_list), list(nodes_list)

    def consistent_color_neighbor_states(
            self,
            state1: List[Set[str]],
            state2: List[Set[str]],
            base_lengths: Dict[str, int],
            color_map: Dict[int, int],
            next_color: int
    ) -> Tuple[Dict[int, int], int]:
        # inputs state1 and state2 are expected to be sorted from largest to smallest component

        len_state1 = len(state1)
        len_state2 = len(state2)

        rank_map1: Dict[int, List] = {i: [] for i in range(len_state1)}
        rank_map2: Dict[int, List] = {j: [] for j in range(len_state2)}

        for i, comp1 in enumerate(state1):
            for j, comp2 in enumerate(state2):
                comp_overlap = comp1.intersection(comp2)
                if len(comp_overlap) > 0:
                    node_overlap = sum(base_lengths[comp_id] for comp_id in comp_overlap)
                    rank_map1[i].append((j, node_overlap))
                    rank_map2[j].append((i, node_overlap))
            # rank descending-order sort (keeping order of larger to smaller component size)
            rank_map1[i] = [x[0] for x in sorted(rank_map1[i], reverse=True, key=lambda x: x[1])]
        for j in range(len_state2):
            # rank descending-order sort
            rank_map2[j] = [x[0] for x in sorted(rank_map2[j], reverse=True, key=lambda x: x[1])]

        new_color_map = {j: -1 for j in range(len_state2)}

        for i in range(len_state1):
            has_mapping = False
            while not has_mapping and len(rank_map1[i]) > 0:
                j = rank_map1[i].pop(0)
                if new_color_map[j] == -1:
                    new_color_map[j] = i
                    has_mapping = True
                elif rank_map2[j].index(i) < rank_map2[j].index(new_color_map[j]):
                    new_color_map[j] = i
                    has_mapping = True

        # assign colors to pairs
        for j in new_color_map.keys():
            if new_color_map[j] == -1:
                # assign entirely new color
                new_color_map[j] = next_color
                next_color += 1
            else:
                # map over color from pairing
                new_color_map[j] = color_map[new_color_map[j]]

        return new_color_map, next_color

    @as_custom_generator(Tuple[str, List[Set[str]], Dict[int, int]], int)
    def consistent_color_components(
            self,
            base_components: ET.Element,
            states: ET.Element
    ) -> CustomGenerator[Tuple[str, List[Set[str]], Dict[int, int]], int]:
        # retrieve base mapping to resolve comp_ids
        base_mapping, base_lengths = self.resolve_base_mapping(base_components)

        initial_state, initial_nodes = self.create_list_from_state(states[0], base_mapping)
        initial_coloring = {i: i + 1 for i in range(len(initial_state))}
        next_color = len(initial_state) + 1

        yield states[0].get(XMLConstants.KEY_ATTRIBUTE_NAME.value, None), initial_nodes, initial_coloring

        for i in range(1, len(states)):
            next_state, next_nodes = self.create_list_from_state(states[i], base_mapping)

            new_color_map, next_color = self.consistent_color_neighbor_states(
                initial_state, next_state, base_lengths, initial_coloring, next_color
            )

            yield states[i].get(XMLConstants.KEY_ATTRIBUTE_NAME.value, None), next_nodes, new_color_map

            initial_state = next_state
            initial_coloring = new_color_map

        return next_color


StructureServiceRegistry.XML.register_service(XmlStructureService())
