"""
This module provides the PebbleGame class and functions for conducting a (Component) Pebble Game on
(`k,l`)-(Multi-)Graphs.
The implementation is based upon the algorithm proposed in
`Lee and Streinu (2008)`_.

 .. _Lee and Streinu (2008):
    https://doi.org/10.1016/j.disc.2007.07.104

"""
from typing import Tuple, Union, List, Optional, Dict

import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

from networkx.utils import not_implemented_for

from TRAMbio.pebble_game.base_pebble_game import BasePebbleGame
from TRAMbio.util.constants.pebble_game import PebbleGameCategory


__all__ = ["PebbleGame", "run_pebble_game", "run_component_pebble_game", "get_evaluation_from_results"]


class PebbleGame(BasePebbleGame):
    """
    Pebble Game Class
    =================

    This class provides methods and functionality for conducting the Component Pebble Game on a given
    (`k,l`)-(Multi-)Graph.

    Parameters:
    -----------
    graph : `networkx.Graph` or `networkx.MultiGraph`
        The input (`k,l`)-(Multi-)Graph on which the Pebble Game will be simulated.
    k : `int`
        The number of degrees of freedom each node adds to the graph.
    l : `int`
        The number of system degrees of freedom.
    component_prefix : `str`, optional (default: 'comp-1A')
        A prefix used for naming identified components.
    use_dfs : `bool`, optional (default: True)
        Flag indicating whether to use Depth-First Search (DFS) if True or Breadth-First Search (BFS) if False.

    Examples:
    ---------
        >>> import networkx as nx
        >>> G = nx.Graph()
        >>> G.add_edges_from([(0, 1), (1, 2), (3, 4)])
        >>> pebble_game = PebbleGame(G, k=2, l=3)
        >>> redundant, excess = pebble_game.play_component_pebble_game(verbose=True)
        >>> components = pebble_game.get_components()
    """

    def __init__(
            self,
            graph: Union[nx.Graph, nx.MultiGraph],
            k: int,
            l: int,
            component_prefix: str = 'comp-1A'
    ):
        node_list = list(graph.nodes)
        super().__init__(k, l, None, node_list, component_prefix)

        if isinstance(graph, nx.MultiGraph):
            self._base_graph = nx.Graph()
            for u, v, data in graph.edges(data=True):
                w = data['weight'] if 'weight' in data else 1.0
                if self._base_graph.has_edge(u, v):
                    self._base_graph[u][v]['weight'] += w
                else:
                    self._base_graph.add_edge(u, v, weight=w)
        else:
            self._base_graph = graph

        self._num_nodes = len(node_list)
        self.time_value = 0.0

        if self.type == 0:
            self._main_component = []
            self._is_same_component = self._is_same_component_l_zero
            self._component_maintenance = self._component_maintenance_l_zero
        elif self.type == 1:
            self._node_labels: pd.Series = pd.Series(data=[np.nan] * self.num_nodes, index=node_list, dtype=object)
            self._is_same_component = self._is_same_component_l_low
            self._component_maintenance = self._component_maintenance_l_low
        else:
            self._edge_labels: Dict[Tuple[any, any], Union[str, None]] = {}
            self._is_same_component = self._is_same_component_l_high
            self._component_maintenance = self._component_maintenance_l_high

        self._pebble_excess: Optional[int] = None

    @property
    def num_nodes(self) -> int:
        return self._num_nodes

    @property
    def pebble_excess(self) -> Optional[int]:
        return self._pebble_excess

    def _is_same_component(self, vertex_u: any, vertex_v: any) -> bool:
        raise NotImplementedError

    def _component_maintenance(self, node_listing: List[any]) -> None:
        raise NotImplementedError

    def _is_same_component_l_zero(self, vertex_u: any, vertex_v: any) -> bool:
        """Checks if two nodes belong to the same component with `l=0`."""
        return all(x in self._main_component for x in [vertex_u, vertex_v])

    def _component_maintenance_l_zero(self, node_listing: List[any]):
        """Performs component maintenance for the Pebble Game with `l=0`."""
        node_listing = list(set(self._main_component + node_listing))
        self._main_component = node_listing

    def _is_same_component_l_low(self, vertex_u: any, vertex_v: any) -> bool:
        """Checks if two nodes belong to the same component with `0<l<=k`."""
        label_u = self._node_labels[vertex_u]
        return isinstance(label_u, str) and label_u == self._node_labels[vertex_v]

    def _component_maintenance_l_low(self, node_listing: List[any]):
        """Performs component maintenance for the Pebble Game with `0<l<=k`."""
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

    def _is_same_component_l_high(self, vertex_u: any, vertex_v: any) -> bool:
        """Checks if two nodes belong to the same component with `l>k`."""
        edge = (min(vertex_u, vertex_v), max(vertex_u, vertex_v))
        return edge in self._edge_labels.keys() and self._edge_labels[edge] is not None

    def _component_maintenance_l_high(self, node_listing: List[any]):
        """Performs component maintenance for the Pebble Game with `l>k`."""
        component_label = f'{self._next_component_prefix}{self._next_component_label}'
        self._next_component_label += 1

        self._edge_labels.update({(min(edge[0], edge[1]), max(edge[0], edge[1])): component_label for edge in
                                  nx.induced_subgraph(self._help_graph, node_listing).edges})

    def play_component_pebble_game(
            self,
            verbose: bool = False,
            desc: str = "Pebble Game") -> Tuple[int, int]:
        """
        Simulates the Component Pebble Game for the edges of the corresponding graph.

        Parameters
        ----------
        verbose : `bool`, optional (default: False)
            If True, display progress information during the simulation.

        desc : `str`, optional (default: 'Pebble Game')
            Description for the progress bar (if verbose=True).

        Returns
        -------
        Tuple[int, int]
            A tuple containing two integers:

            (1) The total count of rejected and redundant edges during the Pebble Game, aka. rho value.

            (2) The pebble excess, which equivalents to the remaining degrees of freedom in the graph
                (minus the number of system degrees of freedom `l`).
        """
        if self._pebble_excess is not None or self.l >= 2 * self.k:
            return -1, -1

        edges = [(u, v, w) for u, v, w in self._base_graph.edges(data="weight", default=1)]
        total_edges = sum(map(lambda x: x[2], edges))

        # counter for excluded edges
        redundant = 0

        with tqdm(total=total_edges, desc=desc, disable=not verbose) as progress:
            for edge in edges:
                newly_redundant = self._play_component_pebble_game_on_edge(edge)

                redundant += newly_redundant
                progress.update(edge[2])
                progress.set_postfix({"redundant": redundant})

        self._pebble_excess = sum(nx.get_node_attributes(self._help_graph, 'pebbles').values()) - self.l

        return self.num_redundant_edges, self._pebble_excess

    def get_components(self) -> List[List[any]]:
        """
        Getter for the calculated rigid components in the corresponding graph.

        Returns
        -------
        `list` [`list`]
            A list of rigid components, where each component is represented as a list of nodes.
        """

        if self._components is not None:
            return self._components

        if self._pebble_excess is None:
            return []

        if self.type == 0:
            self._components = [sorted(self._main_component)] if len(self._main_component) > 0 else []
            return self._components
        elif self.type == 1:
            component_dict = {}

            for node, comp_id in self._node_labels.items():
                node = str(node)
                if not isinstance(comp_id, str):
                    continue
                if comp_id in component_dict.keys():
                    component_dict[comp_id].append(node)
                else:
                    component_dict[comp_id] = [node]

            self._components = [sorted(comp) for comp in component_dict.values()]
            return self._components
        else:
            component_dict = {}
            for edge, comp_id in self._edge_labels.items():
                if comp_id is None:
                    continue
                if comp_id in component_dict.keys():
                    component_dict[comp_id].add(edge[0])
                    component_dict[comp_id].add(edge[1])
                else:
                    component_dict[comp_id] = {edge[0], edge[1]}
            self._components = [list(sorted(comp)) for comp in component_dict.values()]
            return self._components


########################
# Public functions #####
########################

def run_pebble_game(
        graph: Union[nx.Graph, nx.MultiGraph],
        k: int,
        l: int,
        verbose: bool = False) -> Tuple[int, int, str]:
    """
    Calculates the Pebble Game result for the given graph.

    Parameters
    ----------
    graph : `networkx.Graph` or `networkx.MultiGraph`
        The input (`k,l`)-(Multi-)Graph on which the Pebble Game will be simulated.
    k : `int`
        The number of degrees of freedom each node adds to the graph.
    l : `int`
        The number of system degrees of freedom.
    verbose : `bool`, optional (default: False)
       If True, display verbosity output during the Pebble Game.

    Returns
    -------
    Tuple[int, int, str]
        A tuple containing the following elements:

        (1) The total count of rejected and redundant edges during the Pebble Game, aka. rho value.

        (2) The pebble excess, which equivalents to the remaining degrees of freedom in the graph
            (minus the number of system degrees of freedom `l`).

        (3) An evaluation of the Pebble Game result, categorizing the graph as either 'well-constrained',
            'over-constrained', 'under-constrained', or 'other'.

    """
    pebble_game = PebbleGame(graph=graph, k=k, l=l)

    rho, pebble_excess = pebble_game.play_component_pebble_game(verbose=verbose)

    return rho, pebble_excess, get_evaluation_from_results(rho=rho, pebble_excess=pebble_excess)


@not_implemented_for("directed")
def run_component_pebble_game(
        graph: Union[nx.Graph, nx.MultiGraph],
        k: int,
        l: int,
        verbose: bool = False):
    """
    Run the Component Pebble Game on the given graph.

    Parameters
    ----------
    graph : `networkx.Graph` or `networkx.MultiGraph`
        The input (`k,l`)-(Multi-)Graph on which the Component Pebble Game will be simulated.
    k : `int`
        The number of degrees of freedom each node adds to the graph.
    l : `int`
        The number of system degrees of freedom.
    verbose : `bool`, optional (default: False)
        If True, display verbosity output during the Pebble Game.

    Returns
    -------
    Tuple[int, int, List[List[any]]]
        A tuple containing the following elements:

        (1) The total count of rejected and redundant edges during the Pebble Game, aka. rho value.

        (2) The pebble excess, which equivalents to the remaining degrees of freedom in the graph
            (minus the number of system degrees of freedom `l`).

        (3) A list of components, where each component is represented as a list of nodes.

    Raises
    ------
    NetworkXNotImplemented
        If graph is directed.
    """
    pebble_game = PebbleGame(graph=graph, k=k, l=l)

    rho, pebble_excess = pebble_game.play_component_pebble_game(verbose=verbose)

    return rho, pebble_excess, pebble_game.get_components()


def get_evaluation_from_results(rho: int, pebble_excess: int):
    """
    Determines the category descriptor for the Pebble Game result values.

    Parameters
    ----------
    rho : `int`
        The number of rejected or redundant edges during the Pebble Game.
    pebble_excess : `int`
        The pebble excess value calculated in the Pebble Game,
        which equivalents to the remaining degrees of freedom in the graph
        (minus the number of system degrees of freedom `l`).

    Returns
    -------
    str
        A string indicating the evaluation result, which can be one of the following:

        - 'well-constrained' if rho is 0 and pebble_excess is 0.

        - 'over-constrained' if rho is not 0 and pebble_excess is 0.

        - 'under-constrained' if rho is 0 and pebble_excess is not 0.

        - 'other' for any other combination of rho and pebble_excess values.
    """
    if pebble_excess == 0:
        return PebbleGameCategory.WELL_CONSTRAINED.value if rho == 0 else PebbleGameCategory.OVER_CONSTRAINED.value
    else:
        return PebbleGameCategory.UNDER_CONSTRAINED.value if rho == 0 else PebbleGameCategory.OTHER.value
