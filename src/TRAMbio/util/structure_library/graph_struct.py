import sys
from looseversion import LooseVersion
from typing import Tuple, List, Union, Dict

if LooseVersion(sys.version) < LooseVersion("3.11"):
    from typing_extensions import TypedDict, Required
else:
    from typing import TypedDict, Required

import copy
from enum import Enum
import networkx as nx
import pandas as pd
import numpy as np

from loguru import logger

from TRAMbio.util.constants.interaction import InteractionType, map_interaction_ranking


__all__ = [
    "GraphDictionary", "GraphKey", "ProteinGraph", "initialize_graphs_from_dataframe", "add_missing_nodes",
    "copy_graphs_for_dataframe", "process_node_in_template_graph", "add_bonds_from_frame", "export_bond_frame"
]


class GraphDictionary(TypedDict):
    atom: Required[nx.Graph]
    full: Required[nx.Graph]
    pebble: Required[nx.DiGraph]


class GraphKey(Enum):
    STANDARD_EDGES = 'standard_edges'
    COVALENT_EDGES = 'covalent_edges'
    NON_COVALENT_EDGES = 'non_covalent_edges'
    QUANTIFIED_NON_COVALENT_EDGES = 'dilution_edges'


class ProteinGraph:

    def __init__(
            self,
            graphs: GraphDictionary,
            atom_df: pd.DataFrame, others_df: pd.DataFrame,
            heavy_atom_df: pd.DataFrame, hydrogen_df: pd.DataFrame,
            hydrogen_mapping: pd.DataFrame
    ):
        self.__graphs = graphs
        self.__atom_df = atom_df
        self.__others_df = others_df
        self.__heavy_atom_df = heavy_atom_df
        self.__hydrogen_df = hydrogen_df
        self.__hydrogen_mapping = hydrogen_mapping

    @property
    def graphs(self):
        return self.__graphs

    @property
    def atom_df(self):
        return self.__atom_df

    @property
    def others_df(self):
        return self.__others_df

    @property
    def heavy_atom_df(self):
        return self.__heavy_atom_df

    @property
    def hydrogen_df(self):
        return self.__hydrogen_df

    @property
    def hydrogen_mapping(self):
        return self.__hydrogen_mapping


def _fill_atom_graph(
    graph: nx.Graph,
    atom_df: pd.DataFrame,
) -> None:
    atom_nodes = atom_df['node_id']

    graph.add_nodes_from(atom_nodes)
    graph.graph['node_order'] = dict(zip(graph.nodes(), range(1, len(atom_nodes) + 1)))
    nx.set_node_attributes(
        graph, dict(zip(atom_nodes, atom_df["record_name"])), "record_name"
    )
    nx.set_node_attributes(
        graph, dict(zip(atom_nodes, atom_df["residue_name"])), "residue_name"
    )
    nx.set_node_attributes(
        graph, dict(zip(atom_nodes, atom_df["residue_number"])), "residue_number"
    )
    nx.set_node_attributes(graph, dict(zip(atom_nodes, atom_df["atom_name"])), "atom_type")
    nx.set_node_attributes(
        graph, dict(zip(atom_nodes, atom_df["element_symbol"])), "element_symbol"
    )
    nx.set_node_attributes(
        graph, dict(zip(atom_nodes, np.asarray(atom_df[["x_coord", "y_coord", "z_coord"]]))), "coords"
    )


def initialize_graphs_from_dataframe(
        atom_df: pd.DataFrame,
        heavy_atom_df: pd.DataFrame
) -> "GraphDictionary":
    graphs: GraphDictionary = {
        'atom': nx.Graph(),
        'full': nx.Graph(),
        'pebble': nx.DiGraph()
    }

    _fill_atom_graph(graphs['atom'], heavy_atom_df)
    _fill_atom_graph(graphs['full'], atom_df)

    # initialize atom-level graph and pebble graph
    graphs = add_missing_nodes(graphs=graphs, atom_df=heavy_atom_df, clear_key_lists=True)

    return graphs


def add_missing_nodes(
        graphs: GraphDictionary,
        atom_df: pd.DataFrame,
        clear_key_lists: bool = True
) -> "GraphDictionary":
    atom_nodes = atom_df['node_id']

    new_atom_df = atom_df.loc[~atom_df['node_id'].isin(list(graphs['atom'].nodes)), :].copy().reset_index(drop=True)
    if len(new_atom_df) > 0:
        _fill_atom_graph(graphs['atom'], new_atom_df)
        _fill_atom_graph(graphs['full'], new_atom_df)

    # place nodes in pebbled graph
    graphs['pebble'].add_nodes_from(atom_nodes, pebbles=6)
    # insert empty list for each edge type
    if clear_key_lists:
        for edge_type in GraphKey:
            graphs['pebble'].graph[edge_type.value] = []

    return graphs


def copy_graphs_for_dataframe(
        graphs: GraphDictionary,
        atom_df: pd.DataFrame
) -> "GraphDictionary":
    new_graphs: GraphDictionary = {
        'atom': graphs['atom'].copy(),
        'full': graphs['full'].copy(),
        'pebble': graphs['pebble'].copy()
    }

    atom_nodes = atom_df.loc[atom_df['node_id'].isin(list(graphs['atom'].nodes)), :]
    nx.set_node_attributes(
        new_graphs['atom'], dict(zip(atom_nodes.node_id, np.asarray(atom_nodes[["x_coord", "y_coord", "z_coord"]]))),
        "coords"
    )

    atom_nodes = atom_df.loc[atom_df['node_id'].isin(list(graphs['full'].nodes)), :]
    nx.set_node_attributes(
        new_graphs['full'], dict(zip(atom_nodes.node_id, np.asarray(atom_nodes[["x_coord", "y_coord", "z_coord"]]))),
        "coords"
    )

    for edge_type in GraphKey:
        new_graphs['pebble'].graph[edge_type.value] = copy.copy(graphs['pebble'].graph[edge_type.value])

    return new_graphs


def process_node_in_template_graph(
        graphs: GraphDictionary,
        graph_template: List[Union[Tuple[str, int], Tuple[str, int, bool]]],
        node_id: str
):
    # insert intra-residue edge ('atom', 'pebble')
    node_prefix = node_id[:node_id.rfind(':')]

    for i, entry in enumerate(graph_template):
        if len(entry) == 2:
            target, edge_weight = entry
            include = True
        else:
            target, edge_weight, include = entry

        target_node_id = f"{node_prefix}:{target}"
        if not graphs['atom'].has_node(target_node_id):
            # sanity check
             logger.warning(f'Missing target "{target_node_id}" for "{node_id}". Possibly malformed residue.')
        else:
            if include:
                graphs['atom'].add_edge(node_id, target_node_id, kind={InteractionType.COVALENT.value}, base=True)
                graphs['full'].add_edge(node_id, target_node_id, kind={InteractionType.COVALENT.value}, base=True)
            if i == 0:
                # insert multi-edge and corresponding (empty) back-edge
                graphs['pebble'].add_edge(node_id, target_node_id, weight=5)
                graphs['pebble'].add_edge(target_node_id, node_id, weight=0)
                graphs['pebble'].nodes[node_id]['pebbles'] = 1
                # nx.set_node_attributes(graphs['pebble'], {node_id: 1}, name='pebbles')
                if edge_weight == 6:
                    graphs['pebble'].graph[GraphKey.STANDARD_EDGES.value].append((node_id, target_node_id, 1))
            else:
                graphs['pebble'].graph[GraphKey.STANDARD_EDGES.value].append((node_id, target_node_id, edge_weight))


def add_bonds_from_frame(
        graphs: GraphDictionary,
        bond_frame: pd.DataFrame,
        bond_attributes: Dict[str, str] = None,
        pebble_graph_key: str = GraphKey.COVALENT_EDGES.value,
        pebble_graph_weight: int = 5,
        pebble_graph_quantified_keys: List[str] = None
):
    if bond_attributes is None:
        bond_attributes = {'bond_length': 'bond_length'}
    if pebble_graph_quantified_keys is None:
        pebble_graph_quantified_keys = []

    new_edges = 0

    for i, row in bond_frame.iterrows():
        node_1, node_2 = str(row['node_1']), str(row['node_2'])
        bond_type = row['bond_type']

        if not (graphs['atom'].has_node(node_1) and graphs['atom'].has_node(node_2)):
            # sanity check
            continue

        if graphs['atom'].has_edge(node_1, node_2):
            type_set = graphs['atom'].edges[node_1, node_2]["kind"]
            if len(type_set) == 0 or min(map(map_interaction_ranking, type_set)) > map_interaction_ranking(bond_type):
                # overwrite attributes
                for name, key in bond_attributes.items():
                    graphs['atom'].edges[node_1, node_2][name] = row[key]
                    graphs['full'].edges[node_1, node_2][name] = row[key]
            graphs['atom'].edges[node_1, node_2]["kind"].add(bond_type)
            graphs['full'].edges[node_1, node_2]["kind"].add(bond_type)

        else:
            new_edges += 1
            graphs['atom'].add_edge(
                node_1, node_2, kind={bond_type}, **{name: row[key] for name, key in bond_attributes.items()}
            )
            graphs['full'].add_edge(
                node_1, node_2, kind={bond_type}, **{name: row[key] for name, key in bond_attributes.items()}
            )
            graphs['pebble'].graph[pebble_graph_key].append(
                tuple([node_1, node_2, pebble_graph_weight] + [row[key] for key in pebble_graph_quantified_keys])
            )

    return new_edges


def export_bond_frame(graph: nx.Graph, include_base_edges: bool = False):
    bond_list = []
    for node1, node2, data in graph.edges(data=True):
        if 'kind' not in data.keys():
            continue
        if not include_base_edges and 'base' in data.keys():
            continue
        kind = list(data['kind'])
        if len(kind) == 0:
            kind = {InteractionType.UNKNOWN.value}
        main_kind = min(zip(map(map_interaction_ranking, kind), kind))[1]
        energy = data.get('energy', np.nan)
        extra = list(data['extra']) if 'extra' in data.keys() else []
        bond_list.append((
            node1, node2, main_kind, kind, energy, extra
        ))
    if len(bond_list) == 0:
        return pd.DataFrame([], columns=['node1', 'node2', 'type', 'type_set', 'key', 'extra'], dtype=object)
    return pd.DataFrame(bond_list, columns=['node1', 'node2', 'type', 'type_set', 'key', 'extra'])
