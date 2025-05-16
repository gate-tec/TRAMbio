import random
from typing import Tuple

import pytest

import networkx as nx

from TRAMbio.pebble_game.pebble_game import run_pebble_game, PebbleGame, run_component_pebble_game


#######################################################
# Help class for verifying pebble game invariants #####
#######################################################

class TestInvariants:

    @staticmethod
    def assert_pebble_game_invariants(k, l, graph):
        node_set = [x for x in graph.nodes]
        first_node = node_set[0]
        node_set = node_set[1:]

        permutations = 2 ** len(node_set)

        pebble_game = PebbleGame(graph, k, l)

        for x, y, w in graph.edges(data="weight", default=1):
            # gradually play _play_component_pebble_game_on_edge((x, y, w))
            pebble_game._play_component_pebble_game_on_edge((x, y, w))

            # after each edge, for V' sets of V:
            for i in range(permutations):
                v1, v2 = [first_node], []
                p1, p2 = pebble_game._help_graph.nodes[first_node]["pebbles"], 0

                for j, b in enumerate(map(lambda x: int(x) == 1, bin(i + permutations)[3:])):
                    node = node_set[j]
                    pebbles = pebble_game._help_graph.nodes[node]["pebbles"]
                    if b:
                        v1.append(node)
                        p1 += pebbles
                    else:
                        v2.append(node)
                        p2 += pebbles

                e1, o1, e2, o2 = 0, 0, 0, 0

                for edge in pebble_game._help_graph.edges(data="weight"):
                    weight = edge[2]
                    if weight == 0:
                        continue
                    if edge[0] in v1:
                        if edge[1] in v1:
                            e1 += weight
                        else:
                            o1 += weight
                    else:
                        if edge[1] in v1:
                            o2 += weight
                        else:
                            e2 += weight

                # assert invariants hold
                assert p1 + e1 + o1 == k * len(v1)
                if p1 + e1 >= l:
                    assert e1 <= k * len(v1) - l

                if len(v2) == 0:
                    continue

                assert p2 + e2 + o2 == k * len(v2)
                if p2 + e2 >= l:
                    assert e2 <= k * len(v2) - l


###########################################
# Test generally with complete graphs #####
###########################################

class TestCompleteGraphs(TestInvariants):

    TEST_PARAMETERS = [
        # Type 0 (l==0)
        (1, 0, 3, 0, 0, 'well-constrained'),
        # Type 1 (l<=k)
        (1, 1, 3, 1, 0, 'over-constrained'),
        (2, 1, 3, 0, 2, 'under-constrained'),
        (2, 2, 3, 0, 1, 'under-constrained'),
        # Type 2 (k<l<2k)
        (2, 3, 2, 0, 0, 'well-constrained'),
        (2, 3, 3, 0, 0, 'well-constrained')
    ]

    @pytest.mark.parametrize(
        ("k", "l", "num_nodes", "rho", "excess", "state"),
        TEST_PARAMETERS,
        ids=["%d-%d-%d" % (x[0], x[1], x[2]) for x in TEST_PARAMETERS]
    )
    def test_k_l_pebble_game_on_complete_graph(self, k, l, num_nodes, rho, excess, state):
        graph = nx.complete_graph(num_nodes)
        v_rho, v_excess, v_state = run_pebble_game(graph, k, l, verbose=False)

        assert v_rho == rho
        assert v_excess == excess
        assert v_state == state

    @pytest.mark.parametrize(
        ("k", "l", "num_nodes", "rho", "excess", "state"),
        TEST_PARAMETERS,
        ids=["%d-%d-%d" % (x[0], x[1], x[2]) for x in TEST_PARAMETERS]
    )
    def test_k_l_pebble_game_invariants(self, k, l, num_nodes, rho, excess, state):
        self.assert_pebble_game_invariants(k, l, nx.complete_graph(num_nodes))

    def test_error_on_k_too_low(self):
        graph = nx.complete_graph(3)
        with pytest.raises(ValueError, match="Parameter k"):
            run_pebble_game(graph, 0, 1, verbose=False)

    def test_error_on_l_too_high(self):
        graph = nx.complete_graph(3)
        with pytest.raises(ValueError, match="Parameter l"):
            run_pebble_game(graph, 2, 4, verbose=False)


#################################
# Test with specific graphs #####
#################################

class TestSpecificGraphs(TestInvariants):

    TESTED_GRAPHS = [
        (2, 3, nx.Graph([
            (0, 1), (0, 2), (1, 2), (1, 3), (2, 3)
        ]), 0, 0, 1),
        (2, 3, nx.Graph([
            (0, 1), (0, 2), (0, 3),
            (1, 4), (1, 6),
            (2, 3), (2, 4),
            (3, 4), (3, 8),
            (4, 5),
            (5, 6), (5, 7), (5, 8),
            (6, 7),
            (7, 8)
        ]), 0, 0, 1),
        (2, 3, nx.Graph([
            (0, 1), (0, 2), (1, 2), (1, 6), (1, 8), (2, 3), (3, 4), (3, 5), (4, 5), (4, 7), (5, 8), (6, 7), (6, 8),
            (7, 8)
        ]), 0, 1, 6),
        (3, 3, nx.MultiGraph([
            (0, 1), (0, 2), (0, 2), (0, 3), (0, 4), (0, 4),
            (1, 3), (1, 3),
            (2, 3), (2, 3), (2, 3), (2, 4), (2, 5), (2, 5)
        ]), 0, 1, 1),
        (3, 1, nx.MultiGraph([
            (0, 0), (0, 1), (0, 1), (0, 1), (0, 2),
            (1, 1), (1, 2), (1, 2), (1, 3),
            (2, 3),
            (3, 4),
            (4, 4), (4, 4)
        ]), 0, 1, 2),
        (2, 0, nx.MultiGraph([
            (0, 0), (0, 0), (0, 1),
            (1, 1), (1, 2),
            (2, 3),
            (3, 4),
            (4, 5), (4, 6),
            (5, 6), (5, 7), (5, 7),
            (6, 6), (6, 7),
            (7, 7),
            (8, 8), (8, 8)
        ]), 0, 1, 1)
    ]

    @pytest.mark.parametrize(
        ("k", "l", "graph", "rho", "excess", "comps"),
        TESTED_GRAPHS,
        ids=["%d-%d-%d" % (x[0], x[1], x[2].number_of_nodes()) for x in TESTED_GRAPHS]
    )
    def test_specific_graphs(self, k, l, graph, rho, excess, comps):
        v_rho, v_excess, components = run_component_pebble_game(graph, k=k, l=l)
        assert v_rho == rho
        assert v_excess == excess
        assert len(components) == comps

    @pytest.mark.parametrize(
        ("k", "l", "graph", "rho", "excess", "comps"),
        TESTED_GRAPHS,
        ids=["%d-%d-%d" % (x[0], x[1], x[2].number_of_nodes()) for x in TESTED_GRAPHS]
    )
    def test_invariants_specific_graphs(self, k, l, graph, rho, excess, comps):
        self.assert_pebble_game_invariants(k, l, graph)


###################################################
# Test with random variations of Laman graphs #####
###################################################

class TestRandomGraphs:

    @staticmethod
    def create_random_laman_graph(num_nodes: int = 6) -> nx.Graph:
        assert 3 <= num_nodes

        graph = nx.Graph()
        graph.add_nodes_from([0, 1, 2])
        graph.add_edges_from([(0, 1), (0, 2), (1, 2)])

        for node in range(3, num_nodes):
            node_selection = random.sample(range(node), 3)
            sub_graph = nx.induced_subgraph(graph, node_selection)
            if len(sub_graph.edges) > 0 and random.random() > 0.5:
                # Henneberg type II step
                edge_to_remove = random.choice(list(sub_graph.edges))
                graph.remove_edge(edge_to_remove[0], edge_to_remove[1])
            else:
                # Henneberg type I step
                node_selection = random.sample(node_selection, 2)

            graph.add_node(node)
            graph.add_edges_from([(node, x) for x in node_selection])

        return graph

    @staticmethod
    def add_random_edges(graph: nx.Graph, num_edges: int = 1) -> Tuple[nx.Graph, int]:
        assert 1 <= num_edges

        num_nodes = len(graph.nodes)
        inserted_edges = 0
        for _ in range(num_edges):
            node1 = random.randint(0, num_nodes - 1)
            options = [x for x in range(num_nodes) if x != node1 and x not in graph[node1].keys()]
            if len(options) > 0:
                inserted_edges += 1
                node2 = random.choice(options)
                graph.add_edge(node1, node2)

        return graph, inserted_edges

    @staticmethod
    def remove_random_edges(graph: nx.Graph, num_edges: int = 1) -> nx.Graph:
        assert 1 <= num_edges
        assert num_edges <= len(graph.edges)

        graph.remove_edges_from(random.sample(list(graph.edges), num_edges))

        return graph

    def test_well_constrained(self):
        # construct random well-constrained laman graphs
        for _ in range(20):
            laman_graph = self.create_random_laman_graph(random.randint(6, 20))
            rho, pebble_excess, components = run_component_pebble_game(laman_graph, k=2, l=3)
            assert rho == 0
            assert pebble_excess == 0
            assert len(components) == 1

    def test_two_main_components(self):
        # construct two random laman graphs joined by a single edge
        for _ in range(10):
            laman_graph_1 = self.create_random_laman_graph(10)
            laman_graph_2 = self.create_random_laman_graph(10)
            laman_graph_2 = nx.convert_node_labels_to_integers(
                laman_graph_2, first_label=10, ordering='default', label_attribute=None
            )

            combined_graph: nx.Graph = nx.union(laman_graph_1, laman_graph_2)
            combined_graph.add_edge(random.choice(range(10)), random.choice(range(10)) + 10)

            rho, pebble_excess, components = run_component_pebble_game(combined_graph, k=2, l=3)
            assert rho == 0
            assert pebble_excess == 2
            assert len(components) == 3
            # the three components are the two laman graphs and a single edge
            assert len(list(filter(lambda x: len(x) == 10, components))) == 2
            node1, node2 = tuple(sorted(list(filter(lambda x: len(x) == 2, components))[0]))
            assert node1 < 10
            assert node2 >= 10

    def test_overconstrained(self):
        # construct random overconstrained graphs
        for _ in range(10):
            laman_graph = self.create_random_laman_graph(15)

            over_constrained_graph, inserted_edges = self.add_random_edges(laman_graph, 2)

            rho, pebble_excess, components = run_component_pebble_game(over_constrained_graph, k=2, l=3)
            assert rho == inserted_edges
            assert pebble_excess == 0
            assert len(components) == 1

    def test_underconstrained(self):
        # construct random underconstrained graph
        for _ in range(10):
            laman_graph = self.create_random_laman_graph(15)

            num_edges_to_remove = random.randint(2, 4)
            under_constrained_graph = self.remove_random_edges(laman_graph, num_edges_to_remove)

            rho, pebble_excess, components = run_component_pebble_game(under_constrained_graph, k=2, l=3)
            assert rho == 0
            assert pebble_excess == num_edges_to_remove
            # removing any edge from the graph breaks it into multiple components
            assert len(components) > 1 or (len(components) == 1 and len(components[0]) < 15)
