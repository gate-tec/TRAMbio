from contextlib import AbstractContextManager
from typing import Dict
from enum import Enum

import networkx as nx

from TRAMbio.util.constants.interaction import InteractionType


class BondType(Enum):
    SINGLE = 'SINGLE'
    DOUBLE = 'DOUBLE'
    RESONANCE = 'RESONANCE'


class NeighborType(Enum):
    """
    Neighborhood characterization.

    For `P_XY` the two nodes are at the `X` and `Y` position in a chain of nodes.
    """
    P_11 = "11"
    """Type for distance between a node and itself."""
    P_12 = "12"
    """Neighbors with a single edge distance (1-2)"""
    P_13 = "13"
    """Neighbors with a two edge distance (1-2-3)"""
    P_14 = "14"
    """Neighbors with a three edge distance (1-2-3-4)"""
    P_15 = "15"
    """Neighbors with a four edge distance (1-2-3-4-5)"""
    UNBOUNDED = "UNBOUNDED"
    """Neighbors are more than four edges apart."""


class AtomNeighborhood(AbstractContextManager):

    def __init__(self, graph: nx.Graph):
        self._graph = graph
        self._path_lengths = None

    def __enter__(self):
        """Return `self` upon entering the runtime context."""
        self._path_lengths = dict(nx.all_pairs_shortest_path_length(self._graph, 4))  # limit to 1-5 pairs
        return self

    def distance_type(self, atom1: str, atom2: str) -> NeighborType:
        if self._path_lengths is None:
            raise ValueError(f"Neighborhood is not calculated.")

        if atom1 == atom2:
            return NeighborType.P_11
        try:
            if atom2 not in self._path_lengths[atom1]:
                return NeighborType.UNBOUNDED
            length = self._path_lengths[atom1][atom2]
        except KeyError:
            raise KeyError(f"Pair ({atom1},{atom2}) does not exist in graph.")

        if length == 1:
            return NeighborType.P_12
        elif length == 2:
            return NeighborType.P_13
        elif length == 3:
            return NeighborType.P_14
        elif length == 4:
            return NeighborType.P_15
        else:
            return NeighborType.UNBOUNDED

    def __exit__(self, *exec_info):
        """Raise any exception triggered within the runtime context."""
        return None


def get_covalent_neighbors(
        atom_graph: nx.Graph,
        node_id: str
) -> Dict[str, BondType]:
    if atom_graph.has_node(node_id):
        bonds = {}
        for neighbor in atom_graph[node_id]:
            bond_kind = atom_graph.get_edge_data(node_id, neighbor)['kind']

            if bond_kind is None or all(
                    x not in bond_kind
                    for x in [InteractionType.COVALENT.value, InteractionType.PEPTIDE_BOND.value]
            ):
                continue

            # currently only supports single bonds
            bonds[neighbor] = BondType.SINGLE

        return bonds
    else:
        raise KeyError(f"Node {node_id} not in atom_graph")
