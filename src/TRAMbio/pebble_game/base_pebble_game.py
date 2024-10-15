from typing import Tuple, Union, List, Optional, Any, Set
import abc

import networkx as nx


class BasePebbleGame(metaclass=abc.ABCMeta):

    def __init__(
            self,
            k: int,
            l: int,
            help_graph: Optional[nx.DiGraph],
            node_list: Optional[List[Any]],
            component_prefix: str = 'comp-1A'
    ):
        self._k = k
        if k < 1:
            raise ValueError(f"Parameter k for Pebble Game needs to be greater or equal to 1.")
        self._l = l
        if not (0 <= l < 2 * k):
            raise ValueError(f"Parameter l for Pebble Game needs to be in interval [0,2k).")

        self._type = 0 if l == 0 else (1 if l <= k else 2)

        if help_graph is None:
            if node_list is None:
                raise ValueError(f"Either help_graph or node_list must be provided.")

            self._help_graph: nx.DiGraph = nx.DiGraph()
            self._help_graph.add_nodes_from(node_list, pebbles=self._k)
        else:
            self._help_graph: nx.DiGraph = help_graph
            # TODO: validity check

        self._next_component_prefix: str = component_prefix
        self._next_component_label = 1

        self._components: Optional[List] = None
        self._num_redundant_edges = 0

    ##################
    # Properties #####
    ##################

    @property
    def k(self) -> int:
        """`k` property regarding the (`k,l`)-(Multi-)Graph.

        The number of degrees of freedom each node adds to the graph.
        """
        return self._k

    @property
    def l(self) -> int:
        """`l` property regarding the (`k,l`)-(Multi-)Graph.

        The number of system degrees of freedom.
        """
        return self._l

    @property
    def type(self) -> int:
        """An integer indicating the type of the Pebble Game (0, 1, or 2) based on `k` and `l`."""
        return self._type

    @property
    def num_redundant_edges(self) -> int:
        """The number of edges identified as redundant during the Pebble Game run."""
        return self._num_redundant_edges

    #################################
    # General utility functions #####
    #################################

    def _handle_redundant_edges(self, vertex_u: any, vertex_v: any, edges: int):
        self._num_redundant_edges += edges

    @abc.abstractmethod
    def _is_same_component(self,
                           vertex_u: any,
                           vertex_v: any) -> bool:
        """Checks if two nodes belong to the same component."""
        raise NotImplementedError

    def _insert_edge(self, vertex_u: any, vertex_v: any, weight: int):
        if self._help_graph.nodes[vertex_u]['pebbles'] < weight:
            # Safety check
            raise ValueError(
                f'Attempting to insert {weight} edges ({vertex_u},{vertex_v}) with {vertex_u} having not enough pebbles.')
        self._help_graph.nodes[vertex_u]['pebbles'] -= weight
        if self._help_graph.has_edge(vertex_u, vertex_v):
            self._help_graph.edges[vertex_u, vertex_v]['weight'] += weight
        else:
            self._help_graph.add_edge(vertex_u, vertex_v, weight=weight)
            if vertex_u != vertex_v:
                self._help_graph.add_edge(vertex_v, vertex_u, weight=0)

    def _dfs_search(self, source: any, excludes: Optional[Union[List, Set]] = None):
        visited = {source}
        if excludes is not None:
            for excluded_vertex in excludes:
                visited.add(excluded_vertex)

        stack = [(source, iter(self._help_graph.successors(source)))]
        while stack:
            parent, children = stack[-1]
            try:
                child = next(children)
                edge_capacity = self._help_graph.edges[parent, child]['weight']
                if edge_capacity > 0 and child not in visited:
                    yield parent, child
                    visited.add(child)
                    stack.append((child, iter(self._help_graph.successors(child))))
            except StopIteration:
                stack.pop(-1)

    ############################
    # Standard Pebble Game #####
    ############################

    def _find_pebbles(self, vertex_u: any, vertex_v: any) -> bool:
        stack = []
        found_node = None
        for edge in self._dfs_search(vertex_u, [vertex_v]):
            stack.append(edge)
            if self._help_graph.nodes[edge[1]]['pebbles'] > 0:
                found_node = edge[1]
                break
        if found_node is not None:
            # move pebble
            self._help_graph.nodes[found_node]['pebbles'] -= 1
            self._help_graph.nodes[vertex_u]['pebbles'] += 1

            # reverse edges
            for edge in stack[::-1]:
                if edge[1] == found_node:
                    found_node = edge[0]
                    self._help_graph.edges[edge[0], edge[1]]['weight'] -= 1
                    self._help_graph.edges[edge[1], edge[0]]['weight'] += 1
            return True
        else:
            return False

    def _play_pebble_game(self, vertex_u: any, vertex_v: any, edge_weight: int) -> Optional[Tuple[any, any, int]]:
        difference = self._help_graph.nodes[vertex_u]['pebbles'] + \
                     (self._help_graph.nodes[vertex_v]['pebbles'] if vertex_u != vertex_v else 0) - \
                     (self.l + edge_weight)

        searching = True
        while difference < 0 and self._help_graph.nodes[vertex_u]['pebbles'] < self.k and searching:
            # try aggregating pebbles to u
            if self._find_pebbles(vertex_u, vertex_v):
                difference += 1
            else:
                searching = False

        searching = True
        while vertex_u != vertex_v and difference < 0 and self._help_graph.nodes[vertex_v]['pebbles'] < self.k and searching:
            # try aggregating pebbles to v
            if self._find_pebbles(vertex_v, vertex_u):
                difference += 1
            else:
                searching = False

        difference = abs(min(difference, 0))
        if difference < edge_weight:
            weight_to_insert = edge_weight - difference
            weight_for_u = min(self._help_graph.nodes[vertex_u]['pebbles'], weight_to_insert)

            if weight_for_u > 0:
                self._insert_edge(vertex_u, vertex_v, weight_for_u)
            if weight_to_insert - weight_for_u > 0:
                self._insert_edge(vertex_v, vertex_u, weight_to_insert - weight_for_u)
            return vertex_u, vertex_v, weight_to_insert
        return None

    ###############################################
    # Component Pebble Game utility functions #####
    ###############################################

    @abc.abstractmethod
    def _component_maintenance(self, node_listing: List[any]) -> None:
        raise NotImplementedError

    #############################
    # Component Pebble Game #####
    #############################

    def _component_search_1(self, reach: Set):
        visited = set()
        visited.update(reach)
        queue = []

        # Initialize queue
        for node in reach:
            for predecessor in self._help_graph.predecessors(node):
                if predecessor not in visited and self._help_graph.edges[predecessor, node]['weight'] > 0:
                    queue.append(predecessor)

        while queue:
            w = queue.pop(0)
            if w in visited:
                continue
            visited.add(w)
            no_free_pebbles = self._help_graph.nodes[w]['pebbles'] == 0
            if not no_free_pebbles:
                continue
            sub_reach = {w}
            for edge in self._dfs_search(w, reach):
                next_node = edge[1]
                sub_reach.add(next_node)
                if no_free_pebbles and self._help_graph.nodes[next_node]['pebbles'] > 0:
                    no_free_pebbles = False
                    break

            if no_free_pebbles:
                reach.update(sub_reach)
                visited.update(sub_reach)
                for node in sub_reach:
                    for predecessor in self._help_graph.predecessors(node):
                        if predecessor not in visited and self._help_graph.edges[predecessor, node]['weight'] > 0:
                            queue.append(predecessor)

        self._component_maintenance(list(reach))

    def _play_component_pebble_game_on_edge(self, edge) -> int:
        vertex_u, vertex_v, edge_weight = edge[0], edge[1], edge[2]

        newly_rejected = 0

        # reject impossible loops
        if vertex_u == vertex_v and self.k <= self.l - 1:
            self._num_redundant_edges += edge_weight
            return edge_weight

        # same component
        if self._is_same_component(vertex_u, vertex_v):
            # no need to record specific edges, as they are already same component
            self._num_redundant_edges += edge_weight
            return edge_weight

        # attempt Pebble Game
        inserted_edge = self._play_pebble_game(vertex_u, vertex_v, edge_weight)
        if inserted_edge is None:
            # fully redundant interaction
            self._handle_redundant_edges(vertex_u, vertex_v, edge_weight)
            return edge_weight
        elif inserted_edge[2] < edge_weight:
            # not all edges have been added, store rejected amount for later return
            newly_rejected = edge_weight - inserted_edge[2]
            # append overflow to redundant edges
            self._handle_redundant_edges(inserted_edge[0], inserted_edge[1], newly_rejected)

        # loop
        if self._help_graph.nodes[vertex_u]['pebbles'] + \
                (self._help_graph.nodes[vertex_v]['pebbles'] if vertex_u != vertex_v else 0) > self.l:
            return newly_rejected

        # DFS for reach
        still_has_pebble = False
        reach = {vertex_u, vertex_v}
        for dfs_edge in self._dfs_search(vertex_u):
            next_node = dfs_edge[1]
            if next_node != vertex_v and self._help_graph.nodes[next_node]['pebbles'] > 0:
                still_has_pebble = True
                break
            reach.add(next_node)

        if still_has_pebble:
            return newly_rejected

        self._component_search_1(reach)

        return newly_rejected

    @abc.abstractmethod
    def get_components(self) -> List[List[Any]]:
        raise NotImplementedError
