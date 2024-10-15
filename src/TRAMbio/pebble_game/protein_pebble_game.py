from typing import Tuple, List, Any, Iterator, Dict

import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

from TRAMbio.pebble_game.base_pebble_game import BasePebbleGame
from TRAMbio.util.structure_library.components import PebbleGameComponent, PebbleGameResult


__all__ = [
    "ProteinPebbleGame"
]

from TRAMbio.util.structure_library.graph_struct import GraphKey


class ProteinPebbleGame(BasePebbleGame):

    def __init__(self, pebble_graph: nx.DiGraph):
        super().__init__(6, 6, pebble_graph, None)

        node_list = list(pebble_graph.nodes)

        self._node_labels: pd.Series = pd.Series(data=[np.nan] * len(node_list), index=node_list, dtype=object)

        # Storage attributes for results
        self._redundant_edges: Dict[any, Dict[any, int]] = {}
        self._pebble_excess: Dict[str, int] = {}

        # Marker for component update
        self._components_changed = False
        self._halo_changed = False

    @classmethod
    def from_simple_graph(cls, simple_graph: nx.Graph):
        pebble_graph = nx.DiGraph()
        pebble_graph.add_nodes_from(simple_graph.nodes, pebbles=6)
        pebble_graph.graph[GraphKey.COVALENT_EDGES.value] = [(x, y, 5) for x, y in simple_graph.edges]

        return cls(pebble_graph=pebble_graph)

    ##################
    # Properties #####
    ##################

    def _handle_redundant_edges(self, vertex_u: any, vertex_v: any, edges: int):
        super()._handle_redundant_edges(vertex_u, vertex_v, edges)
        vertex_1, vertex_2 = vertex_u, vertex_v

        if vertex_1 in self._redundant_edges.keys():
            if vertex_2 in self._redundant_edges[vertex_1].keys():
                self._redundant_edges[vertex_1][vertex_2] += edges
            else:
                self._redundant_edges[vertex_1][vertex_2] = edges
        else:
            self._redundant_edges[vertex_1] = {vertex_2: edges}

        if vertex_2 in self._redundant_edges.keys():
            if vertex_1 in self._redundant_edges[vertex_2].keys():
                self._redundant_edges[vertex_2][vertex_1] += edges
            else:
                self._redundant_edges[vertex_2][vertex_1] = edges
        else:
            self._redundant_edges[vertex_2] = {vertex_1: edges}

    def _is_same_component(
            self,
            vertex_u: any,
            vertex_v: any
    ) -> bool:
        """Checks if two nodes belong to the same component."""
        label_u = self._node_labels[vertex_u]
        return isinstance(label_u, str) and label_u == self._node_labels[vertex_v]

    ###############################################
    # Component Pebble Game utility functions #####
    ###############################################

    def _component_maintenance(self, node_listing: List[any]):
        """

        Should be identical to :py:meth:`TRAMbio.pebble_game.pebble_game.PebbleGame._component_maintenance_l_low`

        Parameters
        ----------
        node_listing

        Returns
        -------

        """
        self._components_changed = True
        component_label = f'{self._next_component_prefix}{self._next_component_label}'
        self._next_component_label += 1

        # collect overwritten components
        overwritten_components = self._node_labels.loc[node_listing].dropna().unique()
        # gather complete new component
        new_node_listing = list(
            set(node_listing).union(set(self._node_labels[self._node_labels.isin(overwritten_components)].index))
        )

        # bookkeeping
        self._node_labels[new_node_listing] = component_label

    ############################
    # Public Run Functions #####
    ############################

    def play_component_pebble_game_dilution(
            self,
            verbose: bool = False,
            desc: str = "Pebble Game",
            edge_key: str = GraphKey.COVALENT_EDGES.value
    ) -> Iterator[Tuple[Any, List[PebbleGameResult]]]:
        if edge_key in self._pebble_excess.keys() or self.l >= 2 * self.k or \
                edge_key not in self._help_graph.graph.keys() or len(self._help_graph.graph[edge_key]) == 0:
            return iter(())

        # counter for excluded edges
        redundant = 0
        last_edge = None

        sum_edge_weights = sum(map(lambda x: x[2], self._help_graph.graph[edge_key]))
        with tqdm(desc=desc, disable=not verbose, total=sum_edge_weights) as progress_bar:
            for edge in self._help_graph.graph[edge_key]:
                last_edge = edge
                newly_redundant = self._play_component_pebble_game_on_edge(edge)
                self._halo_changed = True
                redundant += newly_redundant
                progress_bar.update(edge[2])
                progress_bar.set_postfix({"redundant": redundant})

                if self._components_changed:
                    yield edge, self.get_components()
                    last_edge = None  # if already yielded on final edge, skip next yield

        if last_edge is not None:
            # always yield component state after final edge was inserted
            yield last_edge, self.get_components()

        self._pebble_excess[edge_key] = sum(nx.get_node_attributes(self._help_graph, 'pebbles').values()) - self.l

    def play_component_pebble_game(
            self,
            verbose: bool = False,
            desc: str = "Pebble Game",
            edge_key: str = GraphKey.COVALENT_EDGES.value
    ) -> None:
        if edge_key in self._pebble_excess.keys() or self.l >= 2 * self.k or \
                edge_key not in self._help_graph.graph.keys() or len(self._help_graph.graph[edge_key]) == 0:
            return None

        # counter for excluded edges
        redundant = 0

        sum_edge_weights = sum(map(lambda x: x[2], self._help_graph.graph[edge_key]))
        with tqdm(desc=desc, disable=not verbose, total=sum_edge_weights) as progress_bar:
            for edge in self._help_graph.graph[edge_key]:
                newly_redundant = self._play_component_pebble_game_on_edge(edge)
                self._halo_changed = True
                redundant += newly_redundant
                progress_bar.update(edge[2])
                progress_bar.set_postfix({"redundant": redundant})

        self._pebble_excess[edge_key] = sum(nx.get_node_attributes(self._help_graph, 'pebbles').values()) - self.l

    ##############
    # Getter #####
    ##############

    def get_components(self) -> List[PebbleGameResult]:
        # Reset flag
        components_changed = self._components_changed or self._halo_changed
        self._components_changed = False
        self._halo_changed = False

        if len(self._pebble_excess.keys()) == 0:
            # no games run yet
            return []

        if components_changed or self._components is None:
            component_dict: Dict[str, PebbleGameComponent] = {}

            for node, comp_id in self._node_labels.items():
                node = str(node)
                if not isinstance(comp_id, str):
                    continue
                if comp_id in component_dict.keys():
                    component_dict[comp_id]["nodes"].append(node)
                else:
                    component_dict[comp_id]: PebbleGameComponent = {"nodes": [node], "halo": set()}

                for node_2 in self._help_graph[node]:
                    # if node_2 < node:
                    #     # only count each possible hinge once
                    #     continue
                    comp_id_2 = self._node_labels[node_2]
                    if comp_id_2 == comp_id:
                        continue

                    edge_count = self._help_graph.edges[node, node_2]["weight"] + self._help_graph.edges[node_2, node][
                        "weight"]
                    if edge_count < 5 and node in self._redundant_edges.keys() and node_2 in self._redundant_edges[
                        node].keys():
                        edge_count += self._redundant_edges[node][node_2]
                    if edge_count >= 5:
                        # hinge between two components
                        component_dict[comp_id]["halo"].add(node_2)
                        if not isinstance(comp_id_2, str):
                            continue
                        if comp_id_2 in component_dict.keys():
                            component_dict[comp_id_2]["halo"].add(node)
                        else:
                            component_dict[comp_id_2] = {"nodes": [], "halo": {node}}

                if node in self._redundant_edges.keys():
                    # node has potential halo
                    for node_2, edges in self._redundant_edges[node].items():
                        if self._help_graph.has_edge(node, node_2):
                            continue
                        comp_id_2 = self._node_labels[node_2]
                        if comp_id_2 != comp_id and edges >= 5:
                            # hinge between two components
                            component_dict[comp_id]["halo"].add(node_2)
                            if not isinstance(comp_id_2, str):
                                continue
                            if comp_id_2 in component_dict.keys():
                                component_dict[comp_id_2]["halo"].add(node)
                            else:
                                component_dict[comp_id_2] = {"nodes": [], "halo": {node}}

            self._components: List[PebbleGameResult] = list(sorted([
                {"size": len(comp["nodes"]), "halo_size": len(comp["halo"]),
                 "nodes": list(sorted(comp["nodes"])), "halo": list(sorted(comp["halo"]))}
                for comp in component_dict.values()], key=lambda x: x["size"]))
        return self._components
