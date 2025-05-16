from typing import Dict, List, Tuple, Any
import abc as abc

import os
import pandas as pd
from loguru import logger

from TRAMbio.services import lock_registry, IOServiceRegistry
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.services.workflow import IBaseWorkflowService, IXtcWorkflowService
from TRAMbio.util.constants.xml import XMLConstants
from TRAMbio.util.constants.graph import RIGID_RING_ATOMS
from TRAMbio.util.structure_library.components import StructureRef, PebbleGameResult, IntermediateComponents
from TRAMbio.util.structure_library.generator import CustomGenerator


__all__ = ["BaseWorkflowService"]


class BaseWorkflowService(IBaseWorkflowService, metaclass=abc.ABCMeta):

    def calculate_components_from_generator(
            self,
            structure_generator: CustomGenerator[Tuple[str, List[PebbleGameResult]], Any],
            temp_path
    ) -> Tuple[Dict[str, IntermediateComponents], List[str]]:
        node_dict: Dict[any, str] = {}  # mapping for smallest component a node occurs in
        comp_dict: Dict[str, IntermediateComponents] = {}
        comp_prefix = 'C-'
        next_comp_index = 0

        xml_io_service = IOServiceRegistry.XML.single_service()
        discard_list = []

        old_comp_list = set()

        with open(temp_path, 'w') as xml_f:  # open temp file
            is_first_frame = True
            for structure_key, components in structure_generator:
                # type mark items
                structure_key: str
                components: List[PebbleGameResult]

                no_change_yet = True
                new_comp_list = set()

                with xml_io_service.state_context(xml_f, structure_key=str(structure_key)):

                    # bookkeeping
                    for component in components:
                        """
                        For each component only four types of situations can occur:

                        (1) The component is completely detached from any existing one.
                        (2) The component is identical to an existing one.
                        (3) The component is a strict subset of an existing component.
                        (4) The component partially overlaps with an existing component.
                           (a) The overlapped component is a strict subset.
                           (b) The overlapped component is not a strict subset.
                        """

                        new_component_id = None
                        num_new_nodes = 0
                        overlap = {}

                        for node_id in component['nodes']:
                            if node_id not in node_dict.keys():
                                num_new_nodes += 1

                                if new_component_id is None:
                                    # create new component id
                                    new_component_id = comp_prefix + str(next_comp_index)
                                    next_comp_index += 1
                                    comp_dict[new_component_id] = IntermediateComponents(
                                        size=1, nodes=[node_id], components=None
                                    )
                                else:
                                    # component already has an id
                                    comp = comp_dict[new_component_id]
                                    comp['nodes'].append(node_id)
                                    comp['size'] += 1

                                # register node in new component
                                node_dict[node_id] = new_component_id

                            else:
                                # register overlap
                                comp_id = node_dict[node_id]
                                if comp_id not in overlap.keys():
                                    overlap[comp_id] = [node_id]
                                else:
                                    overlap[comp_id].append(node_id)

                        if not is_first_frame and no_change_yet and num_new_nodes > 0:
                            # Unregistered nodes imply change
                            no_change_yet = False

                        if num_new_nodes == component['size']:
                            """Type 1"""
                            # all nodes were newly registered
                            xml_io_service.write_temp_xml_fragment(
                                xml_file=xml_f,
                                halo=component['halo'],
                                sub_components=[new_component_id]
                            )
                            new_comp_list.add(new_component_id + "".join(sorted(component['halo'])))
                            continue
                        if new_component_id is None and len(overlap) == 1:
                            comp_id = node_dict[component['nodes'][0]]  # get component from any node
                            if comp_dict[comp_id]['size'] == component['size']:
                                """Type 2"""
                                # fully subsumed
                                xml_io_service.write_temp_xml_fragment(
                                    xml_file=xml_f,
                                    halo=component['halo'],
                                    sub_components=[comp_id]
                                )
                                new_comp_list.add(comp_id + "".join(sorted(component['halo'])))
                                continue

                        sub_components = [] if num_new_nodes == 0 else [
                            new_component_id]  # distinction of Types 3 and 4
                        for comp_id, overlap_nodes in overlap.items():
                            overlap_size = len(overlap_nodes)
                            component_size = comp_dict[comp_id]['size']

                            if overlap_size == component_size:
                                """Type 4b"""
                                # overlapped component is fully subsumed
                                sub_components.append(comp_id)
                                continue
                            elif no_change_yet:
                                # Type 3 and 4a imply change in components
                                no_change_yet = False

                            # components only partially overlap
                            # => split into two components

                            sub_component_1 = comp_prefix + str(next_comp_index)
                            sub_component_2 = comp_prefix + str(next_comp_index + 1)
                            next_comp_index += 2

                            comp_dict[sub_component_1] = IntermediateComponents(
                                size=overlap_size, nodes=overlap_nodes, components=None
                            )
                            sub_components.append(sub_component_1)

                            remaining_nodes = []
                            for node_id in comp_dict[comp_id]['nodes']:
                                if node_id in overlap_nodes:
                                    node_dict[node_id] = sub_component_1
                                else:
                                    remaining_nodes.append(node_id)
                                    node_dict[node_id] = sub_component_2

                            comp_dict[sub_component_2] = IntermediateComponents(
                                size=component_size - overlap_size,
                                nodes=remaining_nodes, components=None
                            )

                            comp_dict[comp_id]['nodes'] = None
                            # add components to mark split
                            comp_dict[comp_id]['components'] = [sub_component_1, sub_component_2]

                        xml_io_service.write_temp_xml_fragment(
                            xml_file=xml_f,
                            halo=component['halo'],
                            sub_components=sub_components
                        )
                        new_comp_list.add("".join(sorted(sub_components)) + "".join(sorted(component['halo'])))

                # compare with previous frame
                if not is_first_frame and no_change_yet and old_comp_list == new_comp_list:
                    discard_list.append(structure_key)
                # update reference
                old_comp_list = new_comp_list

                is_first_frame = False


        return comp_dict, discard_list

    @staticmethod
    def _check_structure(nodes, target_atoms):
        if len(nodes) != len(target_atoms):
            # miss-matching length
            return False

        check_query = {atom: False for atom in target_atoms}
        for node_id in nodes:
            atom_name = node_id[10:]
            has_match = False

            # check if atom name is recognized
            for key in check_query.keys():
                if key == atom_name:
                    check_query[key] = True
                    has_match = True
                    break

            if not has_match:
                return False

        # check if all target atoms are present
        return all(check_query.values())

    def _get_meta_info_on_component(self, component_length, node_list):
        if component_length == 3:
            if self._check_structure(node_list, ['N', 'C', 'O']):
                return f' {XMLConstants.STRUCTURE_ATTRIBUTE_NAME.value}="peptide unit"'
        elif component_length in [5, 6, 7, 9]:
            # potential amino acid ring
            present_aas = list({str(node)[6:9] for node in node_list})
            if len(present_aas) == 1 and \
                    present_aas[0] in RIGID_RING_ATOMS.keys() and \
                    self._check_structure(node_list, RIGID_RING_ATOMS[present_aas[0]]):
                # same residue component => mark as ring structure
                return f' {XMLConstants.STRUCTURE_ATTRIBUTE_NAME.value}="ring"'
            if len(present_aas) == 2 and "PRO" in present_aas and \
                    self._check_structure([node for node in node_list if "PRO" in node],
                                     RIGID_RING_ATOMS["PRO"]) and \
                    self._check_structure([node for node in node_list if "PRO" not in node], ['C', 'O']):
                # proline ring with peptide unit
                return f' {XMLConstants.STRUCTURE_ATTRIBUTE_NAME.value}="ring proline"'

        return ''

    def convert_component_archive_to_mapping(
            self,
            comp_dict: Dict[str, IntermediateComponents],
            hydrogen_dict: Dict[Any, List[Any]]
    ) -> Tuple[List[str], int, Dict[str, StructureRef]]:
        component_mapping: Dict[str, StructureRef] = {}
        base_components = []
        num_base_components = 0

        next_base_component_id = 1
        for comp_id, component in comp_dict.items():
            if component['components'] is None:
                # base component
                num_base_components += 1

                component_length = component['size']

                meta_string = self._get_meta_info_on_component(component_length, component['nodes'])

                node_list = []
                halo_list = []

                for node_id in component['nodes']:
                    node_list.append(
                        f'{XMLConstants.XML_INDENT4.value}<{XMLConstants.NODE_TAG.value}>{node_id}</{XMLConstants.NODE_TAG.value}>'
                    )
                    if node_id in hydrogen_dict.keys():
                        # apply hydrogen-halo to base components
                        for h_id in hydrogen_dict[node_id]:
                            halo_list.append(
                                f'{XMLConstants.XML_INDENT4.value}<{XMLConstants.NODE_TAG.value}>{h_id}</{XMLConstants.NODE_TAG.value}>'
                            )

                full_component_length = component_length + len(halo_list)
                attribute_string = f' id="{next_base_component_id}" size="{full_component_length}"{meta_string}'
                node_string = "\n".join(node_list)
                has_halo = len(halo_list) > 0
                halo_string = (
                    f'{XMLConstants.XML_INDENT3.value}<{XMLConstants.HALO_TAG.value}>\n' +
                    "\n".join(halo_list) +
                    f'\n{XMLConstants.XML_INDENT3.value}</{XMLConstants.HALO_TAG.value}>'
                    if has_halo else
                    f'{XMLConstants.XML_INDENT3.value}<{XMLConstants.HALO_TAG.value}/>')

                base_components.append(
                    f'{XMLConstants.XML_INDENT2.value}<{XMLConstants.COMPONENT_TAG.value}{attribute_string}>\n'
                    f'{XMLConstants.XML_INDENT3.value}<{XMLConstants.NODES_TAG.value}>\n'
                    f'{node_string}\n'
                    f'{XMLConstants.XML_INDENT3.value}</{XMLConstants.NODES_TAG.value}>\n'
                    f'{halo_string}\n'
                    f'{XMLConstants.XML_INDENT2.value}</{XMLConstants.COMPONENT_TAG.value}>\n'
                )

                component_mapping[comp_id]: StructureRef = {
                    'out': [
                        f'{XMLConstants.XML_INDENT5.value}<{XMLConstants.COMPONENT_TAG.value} id="{next_base_component_id}"/>\n'],
                    'stack': None
                }
                next_base_component_id += 1

            else:
                # composed component
                out_list = []
                stack_list = iter(component['components']) if len(component['components']) > 0 else None
                component_mapping[comp_id]: StructureRef = {
                    'out': out_list, 'stack': stack_list
                }

        # Unravel reference tree in component_mapping

        def unravel(component_id: str):
            if component_mapping[component_id]['stack'] is None:
                # already unraveled
                return
            for sub_component_id in component_mapping[component_id]['stack']:
                unravel(sub_component_id)
                component_mapping[component_id]['out'] += component_mapping[sub_component_id]['out']
            component_mapping[component_id]['stack'] = None

        # unravel components
        for comp_id in component_mapping.keys():
            unravel(comp_id)

        return base_components, num_base_components, component_mapping

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def run_pipeline_on_generator(
            self,
            generator: CustomGenerator[Tuple[str, List[PebbleGameResult]], pd.DataFrame],
            out_file: str,
            temp_file: str,
            verbose: bool = True,
            parameter_id: str = ''
    ) -> None:
        is_trajectory = isinstance(self, IXtcWorkflowService)
        xml_io_service = IOServiceRegistry.XML.single_service()

        comp_dict, discard_list = self.calculate_components_from_generator(
            structure_generator=generator,
            temp_path=temp_file
        )
        comp_dict: Dict[str, IntermediateComponents]
        discard_list: List[str]

        # exported for halo correction, as most hydrogen atoms are absent within pebble game graph
        hydrogen_mapping: pd.DataFrame = generator.stop()
        # hydrogen_mapping is a DataFrame with columns ['h_id', 'node_id', 'length']
        hydrogen_dict = hydrogen_mapping.copy()\
                            .loc[:, ['h_id', 'node_id']]\
                            .groupby(['node_id'])\
                            .agg(lambda x: x.to_list())["h_id"].to_dict()
        # hydrogen_dict has structure {node_id: [hydrogen_1, hydrogen_2]}

        base_components, num_base_components, component_mapping = self.convert_component_archive_to_mapping(
            comp_dict,
            hydrogen_dict
        )

        if verbose:
            logger.info("Constructing final XML...")

        xml_io_service.convert_temp_to_xml(
            xml_path=out_file,
            temp_path=temp_file,
            base_components=base_components,
            num_base_components=num_base_components,
            component_mapping=component_mapping,
            is_trajectory=is_trajectory,
            discarded_keys=discard_list,
            parameter_id=parameter_id
        )

        # delete temporary file
        os.remove(temp_file)
