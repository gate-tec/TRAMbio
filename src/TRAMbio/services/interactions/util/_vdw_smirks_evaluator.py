from typing import List, Optional, Dict, Tuple, Generator, Callable

import networkx as nx
import numpy as np
import pandas as pd

from TRAMbio.util.constants.smirnoff import LENNARD_JONES_12_6
from TRAMbio.util.functions.graph_neighborhood import get_covalent_neighbors, BondType


class PotentialEvaluator:

    def __init__(self, atom_graph: nx.Graph, atom_df: pd.DataFrame, hydrogen_mapping: pd.DataFrame):
        self._atom_graph = atom_graph
        self._atom_df = atom_df
        self._hydrogen_mapping = hydrogen_mapping

    _PREFIXES: Dict[str, str] = {
        '#': '_check_element_type',
        'X': '_check_connectivity',
        '+': '_check_positive_charge',
        '-': '_check_negative_charge',
        'H': '_check_hydrogen_neighbors',
        ':': '_check_index',
    }

    _BONDS: Dict[str, List[BondType]] = {
        '-': [BondType.SINGLE],
        '~': [BondType.SINGLE, BondType.DOUBLE, BondType.RESONANCE]
    }

    _ELEMENTS: Dict[str, str] = {
        '1': 'H',
        '3': 'Li',
        '5': 'B',
        '6': 'C',
        '7': 'N',
        '8': 'O',
        '9': 'F',
        '11': 'Na',
        '12': 'Mg',
        '13': 'Al',
        '14': 'Si',
        '15': 'P',
        '16': 'S',
        '17': 'Cl',
        '19': 'K',
        '35': 'Br',
        '37': 'Rb',
        '53': 'I',
        '55': 'Cs'
    }

    def evaluate_node(self, node_id: str) -> Tuple[float, float]:
        index = None
        for i, potential in enumerate(LENNARD_JONES_12_6):
            if all(self._evaluate([potential['smirks']], node_id, None, None)):
                index = i

        if index is not None:
            potential = LENNARD_JONES_12_6[index]
            return potential['epsilon'], potential['rmin_half']
        raise ValueError(f"No matching for {node_id}")

    def _evaluate(self, smirks: List[str], node_id: str, prev_node_id: Optional[str], bond_type: Optional[BondType]) -> Tuple[bool, ...]:
        results = []
        for smirks_string in smirks:
            if not smirks_string.startswith('['):
                bond_types = self._BONDS[smirks_string[0]]
                smirks_string = smirks_string[1:]
                if bond_type is None:
                    raise KeyError
                if bond_type not in bond_types:
                    results.append(False)
                    continue

            result = True  # default

            for func, args, optional in self._parse_smirks(smirks_string):
                sub_eval = func(args, node_id, prev_node_id)
                if optional and sub_eval:
                    break
                elif not optional and not sub_eval:
                    result = False
                    break

            results.append(result)

        return tuple(results)

    def _parse_smirks(
            self,
            smirks_string: str
    ) -> Generator[Tuple[Callable[[str, str, Optional[str]], bool], str, bool], None, None]:
        if not smirks_string.startswith('['):
            raise ValueError

        prefix_list = list(self._PREFIXES.keys())

        i = 1
        while i < len(smirks_string):
            if smirks_string[i] == ']':
                if i < len(smirks_string) - 1:
                    yield self._check_neighbor, smirks_string[i + 1:], False
                break
            elif smirks_string[i] == '*':
                # wildcard operator
                i += 1
            elif smirks_string[i] in prefix_list:
                j = i + 1
                while j < len(smirks_string) and smirks_string[j] not in prefix_list + [']', ',']:
                    j += 1
                optional = smirks_string[j] == ','
                yield self.__getattribute__(self._PREFIXES[smirks_string[i]]), smirks_string[i + 1:j], optional
                i = j + (1 if optional else 0)
            else:
                raise KeyError(f"Unkown symbol {smirks_string[i]}")

    def _check_element_type(self, args: str, node_id: str, prev_node_id: Optional[str]) -> bool:
        target_element = self._ELEMENTS[args]
        element_symbol = self._atom_df.loc[self._atom_df['node_id'] == node_id, 'element_symbol'].values[0]
        return element_symbol == target_element

    def _check_connectivity(self, args: str, node_id: str, prev_node_id: Optional[str]) -> bool:
        target_connectivity = int(args)
        if not self._atom_graph.has_node(node_id):
            # hydrogen atom
            return 1 == target_connectivity

        counter = len(self._get_neighbors(node_id, None))

        return counter == target_connectivity

    def _check_positive_charge(self, args: str, node_id: str, prev_node_id: Optional[str]) -> bool:
        target_charge = int(args)
        charge = self._atom_df.loc[self._atom_df['node_id'] == node_id, 'charge'].values[0]
        if np.isnan(charge):
            charge = 0.0

        return charge == target_charge

    def _check_negative_charge(self, args: str, node_id: str, prev_node_id: Optional[str]) -> bool:
        target_charge = -int(args)
        charge = self._atom_df.loc[self._atom_df['node_id'] == node_id, 'charge'].values[0]
        if np.isnan(charge):
            charge = 0.0

        return charge == target_charge

    def _check_hydrogen_neighbors(self, args: str, node_id: str, prev_node_id: Optional[str]) -> bool:
        target_hydrogen = int(args)
        num_hydrogen = len(self._hydrogen_mapping.loc[self._hydrogen_mapping['node_id'] == node_id, :])

        return target_hydrogen == num_hydrogen

    def _check_index(self, args: str, node_id: str, prev_node_id: Optional[str]) -> bool:
        # currently unused since all smirks only contain index :1 at first node
        return True

    def _get_smirks_for_neighbor(self, smirks_string: str) -> List[str]:
        # handle sequence of (potentially multiple) neighbors for node
        if smirks_string.startswith('('):
            delimiter = smirks_string.index(')')
            if delimiter < 1:
                raise ValueError

            return [smirks_string[1:delimiter]] + self._get_smirks_for_neighbor(smirks_string[delimiter + 1:])

        # only one neighbor required
        return [smirks_string]

    def _has_permutation(self, grid: List[Tuple], width: int, height: int, row: int, block_columns: List[int]) -> bool:
        if row >= height or len(block_columns) == width:
            return len(block_columns) == width
        for column in range(width):
            if column in block_columns:
                continue
            if grid[row][column] and self._has_permutation(grid, width, height, row + 1, block_columns + [column]):
                return True
        return self._has_permutation(grid, width, height, row + 1, block_columns)

    def _check_neighbor(self, smirks_string: str, node_id: str, prev_node_id: Optional[str]) -> bool:
        neighbor_smirks = self._get_smirks_for_neighbor(smirks_string)
        result = tuple([False] * len(neighbor_smirks))
        node_results = []

        for neighbor, bond_type in self._get_neighbors(node_id, prev_node_id).items():
            sub_result = self._evaluate(neighbor_smirks, neighbor, node_id, bond_type)
            node_results.append(sub_result)
            result = tuple(a or b for a, b in zip(result, sub_result))

        if not all(result):
            return False

        return self._has_permutation(node_results, len(result), len(node_results), 0, [])

    def _get_neighbors(self, node_id: str, prev_node_id: Optional[str]) -> Dict[str, BondType]:
        all_neighbors = get_covalent_neighbors(self._atom_graph, node_id)
        if prev_node_id is not None and prev_node_id in all_neighbors:
            del all_neighbors[prev_node_id]
        return all_neighbors
