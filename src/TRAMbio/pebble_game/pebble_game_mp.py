"""
This module provides the MultiPebbleGame class and functions for conducting a faster version of the
(Component) Pebble Game on (`k,l`)-(Multi-)Graphs.

See Also
--------
TRAMbio.pebble_game.pebble_game
"""
import concurrent.futures
import itertools
from typing import Union, List, Tuple

import networkx as nx
import pandas as pd
from networkx.utils import not_implemented_for

from TRAMbio.pebble_game.pebble_game import (
    PebbleGame,
    run_pebble_game,
    run_component_pebble_game,
    get_evaluation_from_results
)
from TRAMbio.pebble_game.utility_pebble_game_mp import (convert_to_sub_graphs_by_seeding)
from loguru import logger
from tqdm import tqdm


__all__ = [
    "MultiPebbleGame",
    "run_pebble_game_mp",
    "run_component_pebble_game_mp"
]


class MultiPebbleGame(PebbleGame):
    """
    Utility class for combining multiple Pebble Games from a partitioned (`k,l`)-(Multi-)Graph.

    Parameters:
    -----------
    graph : `networkx.Graph` or `networkx.MultiGraph`
        The input (`k,l`)-(Multi-)Graph on which the Pebble Game will be simulated.
    k : `int`
        The number of degrees of freedom each node adds to the graph.
    l : `int`
        The number of system degrees of freedom.
    """

    def __init__(
            self,
            graph: Union[nx.Graph, nx.MultiGraph],
            k: int,
            l: int
    ):
        # pass empty graph to superclass
        super().__init__(graph=nx.Graph(), k=k, l=l, component_prefix="multi-A")

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

    def integrate_prev_pebble_game(self, pebble_game: PebbleGame):
        """
        Integrate a pre-run Pebble Game into the combined Pebble Game result.

        Parameters
        ----------
        pebble_game: `PebbleGame`
            Pre-run `PebbleGame` instance to be integrated into Multiple Pebble Game.
        """

        if pebble_game.pebble_excess is None or pebble_game.k != self.k or pebble_game.l != self.l:
            logger.warning("Unable to integrate Pebble Game")
            return

        self._help_graph = nx.union(self._help_graph, pebble_game._help_graph)

        self._num_redundant_edges += pebble_game.num_redundant_edges
        self._num_nodes += pebble_game.num_nodes

        if self.type == 0:
            self._main_component = list(set(self._main_component + pebble_game._main_component))
        elif self.type == 1:
            self._node_labels = pd.concat([self._node_labels, pebble_game._node_labels])
        else:
            self._edge_labels.update(pebble_game._edge_labels)

    def integrate_multiple_prev_pebble_games(self, pebble_games: List[PebbleGame]):
        """
        Integrate a list of multiple pre-run Pebble Game into the combined Pebble Game result.

        Parameters
        ----------
        pebble_games: `list` [`PebbleGame`]
            List of pre-run `PebbleGame` instances to be integrated into Multiple Pebble Game.
        """
        valid_pebble_games = []

        help_graphs = [self._help_graph]

        if self.type == 0:
            main_components = [self._main_component]
        elif self.type == 1:
            node_labels = [self._node_labels]
        else:
            edge_labels = [self._edge_labels]

        for i, pebble_game in enumerate(pebble_games):
            if pebble_game.pebble_excess is None or pebble_game.k != self.k or pebble_game.l != self.l:
                logger.warning(f"Skipping incompatible Pebble Game {i}")
                continue
            valid_pebble_games.append(pebble_game)

            help_graphs.append(pebble_game._help_graph)

            self._num_redundant_edges += pebble_game.num_redundant_edges
            self._num_nodes += pebble_game.num_nodes

            if self.type == 0:
                main_components.append(pebble_game._main_component)
            elif self.type == 1:
                node_labels.append(pebble_game._node_labels)
            else:
                edge_labels.append(pebble_game._edge_labels)

        self._help_graph = nx.union_all(help_graphs)

        if self.type == 0:
            self._main_component = list(itertools.chain.from_iterable(comp for comp in main_components))
        elif self.type == 1:
            self._node_labels = pd.concat(node_labels)
        else:
            self._edge_labels = dict(itertools.chain.from_iterable(x.items() for x in edge_labels))


#############################################
# Utility functions for multiprocessing #####
#############################################

def _pebble_game_function(
        graph: Union[nx.Graph, nx.MultiGraph],
        k: int,
        l: int,
        component_prefix: str
) -> PebbleGame:
    """
    Help function for running a single Pebble Game on the given (`k,l`)-(Multi-)Graph.

    Parameters
    ----------
    graph : `networkx.Graph` or `networkx.MultiGraph`
        The input (`k,l`)-(Multi-)Graph on which the Pebble Game will be simulated.
    k : `int`
        The number of degrees of freedom each node adds to the graph.
    l : `int`
        The number of system degrees of freedom.
    component_prefix : `str`, optional (default: 'comp-1A')
        A prefix used for naming identified components.

    Returns
    -------
    `PebbleGame`
        Pebble Game object in executed state.
    """

    pebble_game = PebbleGame(graph=graph, k=k, l=l, component_prefix=component_prefix)

    _, _ = pebble_game.play_component_pebble_game(verbose=False)

    return pebble_game


def _run_on_multiple_sub_graphs(
        sub_graphs: List[nx.MultiGraph],
        k: int,
        l: int,
        connecting_graph: nx.MultiGraph,
        verbose: bool = False,
        threads: int = 2
) -> Tuple[int, int, MultiPebbleGame]:
    """
    Run the Component Pebble Game on a graph partitioning, using multiprocessing.

    Parameters
    ----------
    sub_graphs: `list` [`networkx.MultiGraph`]
        The graph partitioning. List of node-disjoint sub-graphs.
    k : `int`
        The number of degrees of freedom each node adds to the graph.
    l : `int`
        The number of system degrees of freedom.
    connecting_graph: `networkx.MultiGraph`
        Graph connecting the individual sub-graphs.
    verbose : `bool`, optional (default: False)
        If True, display verbosity output during the Pebble Game.
    threads: `int`, optional (default: 2)
        Number of threads for multiprocessing.

    Returns
    -------
    Tuple[int, int, MultiPebbleGame]
        A tuple containing the following elements:

        (1) The total count of rejected and redundant edges during the Pebble Game, aka. rho value.

        (2) The pebble excess, which equivalents to the remaining degrees of freedom in the graph
            (minus the number of system degrees of freedom `l`).

        (3) The combined `MultiPebbleGame` with the final Pebble Game result.

    Raises
    ------
    ValueError
        If number of threads is less than 1.
    """

    if threads < 1:
        raise ValueError("Number of threads for multiprocessing needs to be greater or equal to 1")

    resulting_pebble_games = []
    processed_edges = 0
    redundant_edges = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor, \
            tqdm(total=len(sub_graphs), desc='Sub graphs', disable=not verbose) as progress:
        work_load = {executor.submit(
            _pebble_game_function,
            sub_graphs[i],
            k,
            l,
            f"comp-{i}A"
        ): i for i in range(len(sub_graphs))}

        for future in work_load:
            future.add_done_callback(lambda p: progress.update())
        for future in concurrent.futures.as_completed(work_load):
            index = work_load[future]
            try:
                sub_graph_pebble_game: PebbleGame = future.result()
            except Exception as exc:
                # catch thrown exceptions and store for printing later on
                progress.write('%d generated an exception: %s' % (index, exc))
            else:
                resulting_pebble_games.append(sub_graph_pebble_game)
                processed_edges += k * sub_graph_pebble_game.num_nodes - l - sub_graph_pebble_game.pebble_excess
                redundant_edges += sub_graph_pebble_game.num_redundant_edges
                progress.set_postfix(
                    {"processed": processed_edges, "redundant": redundant_edges}
                )

    multi_pebble_game = MultiPebbleGame(graph=connecting_graph, k=k, l=l)
    multi_pebble_game.integrate_multiple_prev_pebble_games(resulting_pebble_games)

    rho, pebble_excess = multi_pebble_game.play_component_pebble_game(verbose=verbose, desc='Joining edges')

    return rho, pebble_excess, multi_pebble_game


########################
# Public functions #####
########################


@not_implemented_for("directed")
def run_pebble_game_mp(
        graph: Union[nx.Graph, nx.MultiGraph],
        k: int,
        l: int,
        verbose: bool = False,
        threads: int = 1):
    """
    Calculates the Pebble Game result for the given graph using a faster approach.

    Parameters
    ----------
    graph: `networkx.Graph` or `networkx.MultiGraph`
        The base (`k,l`)-(Multi-)Graph for the Pebble Game.
    k : `int`
        The number of degrees of freedom each node adds to the graph.
    l : `int`
        The number of system degrees of freedom.
    verbose : `bool`, optional (default: False)
        If True, display verbosity output during the Pebble Game.
    threads: `int`, optional (default: 2)
        Number of threads for multiprocessing.

    Returns
    -------
    Tuple[int, int, str]
        A tuple containing the following elements:

        (1) The total count of rejected and redundant edges during the Pebble Game, aka. rho value.

        (2) The pebble excess, which equivalents to the remaining degrees of freedom in the graph
            (minus the number of system degrees of freedom `l`).

        (3) An evaluation of the Pebble Game result, categorizing the graph as either 'well-constrained',
            'over-constrained', 'under-constrained', or 'other'.

    Raises
    ------
    ValueError
        If number of threads is less than 1.
    """

    if threads < 1:
        raise ValueError("Number of threads for multiprocessing needs to be greater or equal to 1")

    if threads == 1 or len(graph.edges) / graph.number_of_nodes() < 250:
        return run_pebble_game(graph=graph, k=k, l=l, verbose=verbose)

    graph_partitioning = convert_to_sub_graphs_by_seeding(G=graph, by='greedy')

    if graph_partitioning is None:
        return run_pebble_game(graph=graph, k=k, l=l, verbose=verbose)

    sub_graphs, node_dict, connecting_graph = graph_partitioning

    rho, pebble_excess, _ = _run_on_multiple_sub_graphs(
        sub_graphs=sub_graphs, k=k, l=l, connecting_graph=connecting_graph,
        threads=threads, verbose=verbose
    )
    return rho, pebble_excess, get_evaluation_from_results(rho=rho, pebble_excess=pebble_excess)


@not_implemented_for("directed")
def run_component_pebble_game_mp(
        graph: Union[nx.Graph, nx.MultiGraph],
        k: int,
        l: int,
        verbose: bool = False,
        threads: int = 1):
    """
    Run a fast version of the Component Pebble Game on the given graph.

    Parameters
    ----------
    graph: `networkx.Graph` or `networkx.MultiGraph`
        The base (`k,l`)-(Multi-)Graph for the Component Pebble Game.
    k : `int`
        The number of degrees of freedom each node adds to the graph.
    l : `int`
        The number of system degrees of freedom.
    verbose : `bool`, optional (default: False)
        If True, display verbosity output during the Pebble Game.
    threads: `int`, optional (default: 2)
        Number of threads for multiprocessing.

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
    ValueError
        If number of threads is less than 1.
    """

    if threads < 1:
        raise ValueError("Number of threads for multiprocessing needs to be greater or equal to 1")

    if threads == 1 or graph.number_of_nodes() < 150:
        return run_component_pebble_game(graph=graph, k=k, l=l, verbose=verbose)

    graph_partitioning = convert_to_sub_graphs_by_seeding(G=graph, by='greedy')

    if graph_partitioning is None:
        return run_component_pebble_game(graph=graph, k=k, l=l, verbose=verbose)

    sub_graphs, node_dict, connecting_graph = graph_partitioning

    rho, pebble_excess, multi_pebble_game = _run_on_multiple_sub_graphs(
        sub_graphs=sub_graphs, k=k, l=l, connecting_graph=connecting_graph,
        threads=threads, verbose=verbose
    )

    return rho, pebble_excess, multi_pebble_game.get_components()
