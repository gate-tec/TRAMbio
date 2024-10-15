"""
Utility functions for heuristic (`k,l`)-(Multi-)Graph partitioning.
"""
import random
from typing import Union, Dict, List, Tuple, Optional, Literal, Iterable

import networkx as nx
from networkx.utils import not_implemented_for
from networkx.algorithms.community.modularity_max import greedy_modularity_communities

__all__ = [
    "convert_to_sub_graphs_by_seeding",
    "convert_to_sub_graphs",
    "select_seeds_randomized",
    "select_seeds_by_components_greedy"
]


######################
# Edge functions #####
######################


def shuffler(iterable: Iterable):
    listing = list(iterable)
    random.shuffle(listing)
    yield from listing


@not_implemented_for("directed")
def edge_random_search(G: Union[nx.Graph, nx.MultiGraph], source: List):
    nodes = list(G.nbunch_iter(source))
    if not nodes:
        return

    # start BFS
    visited_nodes = set(nodes)
    visited_edges = set()
    queue = [(n, iter(shuffler(G.edges(n)))) for n in nodes]

    queue_size = len(queue)
    total = 0
    count = 0
    while queue:
        index = random.randint(0, queue_size - 1)
        parent, children = queue[index]
        total += 1
        if index == 0:
            count += 1
        try:
            edge = next(children)

            child = edge[1]
            if child not in visited_nodes:
                visited_nodes.add(child)
                queue.append((child, iter(shuffler(list(G.edges(child))))))
                queue_size += 1
            edgeid = (frozenset(edge[:2]),) + edge[2:]
            if edgeid not in visited_edges:
                visited_edges.add(edgeid)
                yield edge

        except StopIteration:
            queue.pop(index)
            queue_size -= 1

    print(count / total * 100.0)


#########################
# Seeding functions #####
#########################


@not_implemented_for("directed")
def select_seeds_randomized(
        G: Union[nx.Graph, nx.MultiGraph]
) -> Optional[Tuple[List[any], List[nx.MultiGraph], Dict[any, int]]]:
    """
    Select random seed nodes from graph. Number of seeds is 5% of the number of nodes in the graph.

    Parameters
    ----------
    G: `networkx.Graph` or `networkx.MultiGraph`
        The base (`k,l`)-(Multi-)Graph from which to select seed nodes.

    Returns
    -------
    Optional[Tuple[List, List[networkx.MultiGraph], Dict]]
        `None` if no seeding is possible. Otherwise, a tuple containing the following elements:

        (1) The list of seed nodes.

        (2) The list of empty sub-graphs. One per seed.

        (3) The partitioning dictionary mapping selected nodes to sub-graphs.

    Raises
    ------
    NetworkXNotImplemented
        If graph is directed.
    """

    nodes = list(G.nodes)
    num_choices = int(len(nodes) * 0.05)

    if num_choices < 2:
        return None

    source_nodes = random.choices(nodes, k=num_choices)
    sub_graphs = [nx.MultiGraph() for _ in source_nodes]
    node_dict = {source_nodes[i]: i for i in range(len(source_nodes))}

    return source_nodes, sub_graphs, node_dict


@not_implemented_for("directed")
def select_seeds_by_components_greedy(
        G: Union[nx.Graph, nx.MultiGraph]
) -> Optional[Tuple[List[any], List[nx.MultiGraph], Dict[any, int]]]:
    """
    Select seed nodes by utilizing a partition into disjoint modularity communities.

    Parameters
    ----------
    G: `networkx.Graph` or `networkx.MultiGraph`
        The base (`k,l`)-(Multi-)Graph from which to select seed nodes.

    Returns
    -------
    Tuple[List, List[networkx.MultiGraph], Dict]
        A tuple containing the following elements:

        (1) The list of seed nodes.

        (2) The list of empty sub-graphs. One per seed.

        (3) The partitioning dictionary mapping selected nodes to sub-graphs.

    Raises
    ------
    NetworkXNotImplemented
        If graph is directed.

    See Also
    --------
    networkx.algorithms.community.modularity_max.greedy_modularity_communities
    """

    node_sets: List[frozenset] = greedy_modularity_communities(G, resolution=1, cutoff=8, best_n=25)

    source_nodes = []
    sub_graphs = []
    node_dict = {}
    for i, node_set in enumerate(node_sets):
        representative = next(iter(node_set))
        source_nodes.append(representative)
        sub_graphs.append(nx.MultiGraph())
        node_dict[representative] = i

    return source_nodes, sub_graphs, node_dict


####################################
# Sub-graph creation functions #####
####################################


@not_implemented_for("directed")
def convert_to_sub_graphs_by_seeding(
        G: Union[nx.Graph, nx.MultiGraph],
        by: Literal['greedy', 'random'] = 'greedy'
) -> Optional[Tuple[List[nx.MultiGraph], Dict[any, int], nx.MultiGraph]]:
    """
    Partitions a (`k,l`)-(Multi-)Graph into sub-graphs and a collection of binding edges.

    Parameters
    ----------
    G: `networkx.Graph` or `networkx.MultiGraph`
        The base (`k,l`)-(Multi-)Graph to partition.
    by: {'greedy', 'random'}, optional (default: 'greedy')
        Category by which sub-graphs are determined.

    Returns
    -------
    Optional[Tuple[List[nx.MultiGraph], Dict, nx.MultiGraph]]
        `None` if no partitioning is possible. Otherwise, a tuple containing the following elements:

        (1) The list of filled, node-disjoint sub-graphs.

        (2) The updated partitioning dictionary containing all nodes from the base graph.

        (3) The graph containing the remaining edges which connect the individual sub-graphs.

    Raises
    ------
    NetworkXNotImplemented
        If graph is directed.
    ValueError
        If edge_multiplier is less than 1.
    """

    result = None
    if by == 'greedy':
        result = select_seeds_by_components_greedy(G=G)
    elif by == 'random':
        result = select_seeds_randomized(G=G)

    if result is None:
        return None
    return convert_to_sub_graphs(
        G=G,
        source_nodes=result[0],
        sub_graphs=result[1],
        node_dict=result[2]
    )


@not_implemented_for("directed")
def convert_to_sub_graphs(
        G: Union[nx.Graph, nx.MultiGraph],
        source_nodes: List[any],
        sub_graphs: List[nx.MultiGraph],
        node_dict: Dict[any, int]
) -> Tuple[List[nx.MultiGraph], Dict[any, int], nx.MultiGraph]:
    """
    Partitions a (`k,l`)-(Multi-)Graph into sub-graphs

    Parameters
    ----------
    G: `networkx.Graph` or `networkx.MultiGraph`
        The base (`k,l`)-(Multi-)Graph to partition.
    source_nodes: `list`
        The list of seed nodes.
    sub_graphs: `list` [`networkx.MultiGraph`]
        List of pre-created sub-graphs. One per seed node.
    node_dict: `dict`
        The partitioning dictionary mapping selected nodes to sub-graphs.

    Returns
    -------
    Tuple[List[nx.MultiGraph], Dict, nx.MultiGraph]
        A tuple containing the following elements:

        (1) The list of filled, node-disjoint sub-graphs.

        (2) The updated partitioning dictionary containing all nodes from the base graph.

        (3) The graph containing the remaining edges which connect the individual sub-graphs.

    Raises
    ------
    NetworkXNotImplemented
        If graph is directed.
    ValueError
        If edge_multiplier is less than 1.
    """

    # partition nodes using bfs
    # NOTE: nx.edge_bfs has tendency to not generate ALL edges from graph
    for edge in edge_random_search(G=G, source=source_nodes):
        start_node_graph_id = node_dict[edge[0]]
        if edge[1] not in node_dict.keys():
            # place in same sub graph
            node_dict[edge[1]] = start_node_graph_id

    # create backup subgraph for missing nodes
    remaining_nodes = {}
    num_sub_graphs = len(sub_graphs)
    for node in G:
        if node not in node_dict.keys():
            remaining_nodes[node] = num_sub_graphs
    if len(remaining_nodes) > 0:
        num_sub_graphs += 1
        node_dict.update(remaining_nodes)
        sub_graphs.append(nx.MultiGraph())

    # Graph for edges between two sub graphs
    binding_graph = nx.MultiGraph()
    binding_edges = dict()

    # Fill sub graphs
    for edge in G.edges:
        start_node_graph_id = node_dict[edge[0]]
        end_node_graph_id = node_dict[edge[1]]

        if end_node_graph_id == start_node_graph_id:
            sub_graphs[start_node_graph_id].add_edge(edge[0], edge[1])
        else:
            sub_graphs[start_node_graph_id].add_node(edge[0])
            sub_graphs[end_node_graph_id].add_node(edge[1])
            id_tuple = (min(start_node_graph_id, end_node_graph_id), max(start_node_graph_id, end_node_graph_id))
            if id_tuple not in binding_edges.keys():
                binding_edges[id_tuple] = []
            binding_edges[id_tuple].append((edge[0], edge[1]))

    # Order edges between sub graphs for runtime efficiency
    sub_graph_sizes = {i: len(sub_graphs[i].edges) for i in range(len(sub_graphs))}
    print(sub_graph_sizes)
    while len(binding_edges.keys()) > 0:
        lowest_tuple = None
        lowest_tuple_edge_count = 0
        for tuple_id in binding_edges.keys():
            edge_count = sub_graph_sizes[tuple_id[0]] + sub_graph_sizes[tuple_id[1]] + len(binding_edges[tuple_id])
            if lowest_tuple is None or lowest_tuple_edge_count > edge_count:
                lowest_tuple = tuple_id
                lowest_tuple_edge_count = edge_count

        sub_graph_sizes[lowest_tuple[0]] = lowest_tuple_edge_count
        sub_graph_sizes[lowest_tuple[1]] = lowest_tuple_edge_count
        binding_graph.add_edges_from(binding_edges.pop(lowest_tuple))

    return sub_graphs, node_dict, binding_graph
