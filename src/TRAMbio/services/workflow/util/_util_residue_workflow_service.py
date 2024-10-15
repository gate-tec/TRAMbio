from functools import cmp_to_key
from typing import Dict, List, Tuple, Optional

import xml.etree.ElementTree as ET

import pandas as pd
from loguru import logger

from TRAMbio.util.constants.graph import STANDARD_RESIDUES
from TRAMbio.util.constants.xml import XMLConstants, components_namespace


__all__ = ["export_base_map", "export_residue_states"]


def _custom_compare(item1: Tuple[List[str], int], item2: Tuple[List[str], int]) -> int:
    if item1[1] != item2[1]:
        return item1[1] - item2[1]
    for elem1, elem2 in zip(item1[0], item2[0]):
        if elem1 < elem2:
            return -1
        if elem1 > elem2:
            return 1
    return 0


def _compare_states(state_components1: List[Tuple[List[str], int]], state_components2: List[Tuple[List[str], int]]) -> int:
    len1, len2 = len(state_components1), len(state_components2)

    # Compare number of components
    if len1 != len2:
        return len1 - len2

    # Compare length of corresponding components
    for comp1, comp2 in zip(state_components1, state_components2):
        if comp1[1] != comp2[1]:
            return comp1[1] - comp2[1]

    # Compare content of corresponding components
    for comp1, comp2 in zip(state_components1, state_components2):
        cmp_val = _custom_compare(comp1, comp2)
        if cmp_val != 0:
            return cmp_val

    return 0


def export_base_map(
        base_components: ET.Element,
        resi_counts: Optional[pd.Series],
        threshold: float,
        use_main_chain: bool
) -> Dict[str, Dict[str, float]]:
    base_map = {}
    for component in base_components:
        comp_id = component.attrib['id']
        base_comp: Dict[str, float] = {}

        nodes = component.find(*components_namespace(XMLConstants.NODES_TAG))
        if nodes is None:
            logger.warning(f'No <{XMLConstants.NODES_TAG.value}> found for base component {comp_id}.')
            continue

        for node in nodes:
            atom_id = node.text
            atom_parts = atom_id.strip().split(':')
            resi_name = atom_id[6:9].strip()
            if use_main_chain and resi_name in STANDARD_RESIDUES:
                if atom_parts[1] in ['CA', 'C', 'N']:
                    if atom_parts[0] not in base_comp.keys():
                        base_comp[atom_parts[0]] = 1
                    else:
                        base_comp[atom_parts[0]] += 1
            elif resi_counts is not None:
                chain_id = atom_id[0]
                resi_number = int(atom_id[1:5])

                count = int(resi_counts.loc[(chain_id, resi_number), 'node_id'])

                value = 3 / (count * threshold)
                if atom_parts[0] not in base_comp.keys():
                    base_comp[atom_parts[0]] = value
                else:
                    base_comp[atom_parts[0]] += value

        base_map[comp_id] = base_comp

    return base_map


def export_residue_states(
        states: ET.Element,
        base_map: Dict[str, Dict[str, float]],
        min_key: Optional[str],
        min_key_val: Optional[float],
        is_trajectory: bool
) -> Tuple[List[str], Dict[str, List[List[str]]]]:
    components = {}
    components_counter = 0
    previous_state = None
    key_list = []
    for state in states:
        state_key = state.get(XMLConstants.KEY_ATTRIBUTE_NAME.value)
        if state_key is None:
            raise KeyError("Missing key for state")

        if min_key is not None and state_key != "-INF":
            if min_key == "-INF":
                # include no hydrogen bonds
                break
            if min_key_val is not None and (float(state_key) - min_key_val) * (1 if is_trajectory else -1) < 0:
                break

        key_list.append(state_key)

        state_component_list = []

        for component in state:
            sub_components = component.find(*components_namespace(XMLConstants.COMPONENTS_TAG))
            if sub_components is None:
                logger.warning(f'Empty component in state {state_key}.')
                continue

            residue_buffer: Dict[str, float] = {}

            for sub_component in sub_components:
                comp_id = sub_component.attrib['id']
                # base_comp is mapping of residues to number of present key atoms
                base_comp = base_map[comp_id]

                # only consider base_comp.keys()
                overlapping_keys = residue_buffer.keys() & base_comp.keys()
                for key in overlapping_keys:
                    # increment existing residues
                    residue_buffer[key] += base_comp[key]
                for key in base_comp.keys() - overlapping_keys:
                    # add new residues
                    residue_buffer[key] = base_comp[key]

            if len(residue_buffer) == 0:
                continue

            # filter for residues with all three key atoms present
            residue_list = list(sorted(key for key, value in residue_buffer.items() if value >= 3))

            # discard empty components
            if len(residue_list) > 0:
                state_component_list.append((residue_list, len(residue_list)))

        state_component_list = list(sorted(state_component_list, key=cmp_to_key(_custom_compare)))

        # discard unchanged states (w.r.t. coarse-grained residue focus)
        if previous_state is None or _compare_states(previous_state, state_component_list) != 0:
            components[state_key] = [entry for entry, length in state_component_list]
            components_counter += 1
            previous_state = state_component_list

    return key_list, components
