import pandas as pd
import pytest

import networkx as nx
from TRAMbio.pebble_game.protein_pebble_game import ProteinPebbleGame
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import GraphKey, add_bonds_from_frame, add_missing_nodes
from TRAMbio.util.wrapper.base.list import first
from tests.util.graphs import construct_protein_graph_base, graph_gly_donor, graph_gly_acceptor, graph_arg_donor, \
    graph_glu_acceptor


#######################################################
# Help class for verifying pebble game invariants #####
#######################################################

class TestInvariants:

    @staticmethod
    def assert_pebble_game_invariants(
            pebble_game: ProteinPebbleGame,
            graph: nx.DiGraph,
            edge_key=GraphKey.STANDARD_EDGES.value
    ):
        node_set = [x for x in graph.nodes]
        first_node = node_set[0]
        node_set = node_set[1:]

        permutations = 2 ** len(node_set)

        for e in graph.graph[edge_key]:
            x, y, w = e[0], e[1], e[2]
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
                assert p1 + e1 + o1 == 6 * len(v1)
                if p1 + e1 >= 6:
                    assert e1 <= 6 * len(v1) - 6

                if len(v2) == 0:
                    continue

                assert p2 + e2 + o2 == 6 * len(v2)
                if p2 + e2 >= 6:
                    assert e2 <= 6 * len(v2) - 6

        pebble_game._pebble_excess[edge_key] = sum(nx.get_node_attributes(pebble_game._help_graph, 'pebbles').values()) - 6
        pebble_game._halo_changed = True


###################
# Mock Graphs #####
###################

class TestMockGraphs:

    @pytest.fixture
    def mock_graph_pro_pentagon(self):
        edges = [
            ("A0001-PRO:N", "A0001-PRO:CA", {"weight": 5}),
            ("A0001-PRO:CD", "A0001-PRO:N", {"weight": 5}),
            ("A0001-PRO:CG", "A0001-PRO:CB", {"weight": 5}),
            ("A0001-PRO:CB", "A0001-PRO:CA", {"weight": 5}),
        ]

        pebbled_graph = nx.DiGraph([
            edge for n1, n2, d in edges for edge in ((n1, n2, d), (n2, n1, dict(d, weight=0)))
        ])

        pebbles = {
            "A0001-PRO:N": 1,
            "A0001-PRO:CD": 1,
            "A0001-PRO:CG": 1,
            "A0001-PRO:CB": 1,
            "A0001-PRO:CA": 6
        }

        nx.set_node_attributes(pebbled_graph, pebbles, "pebbles")

        pebbled_graph.graph[GraphKey.STANDARD_EDGES.value] = [
            ("A0001-PRO:CD", "A0001-PRO:CG", 5)
        ]
        pebbled_graph.graph[GraphKey.COVALENT_EDGES.value] = []
        pebbled_graph.graph[GraphKey.NON_COVALENT_EDGES.value] = []
        pebbled_graph.graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value] = []

        yield pebbled_graph

    @pytest.fixture
    def mock_graph_phe_hexagon(self):
        edges = [
            ("A0001-PHE:CZ", "A0001-PHE:CE1", {"weight": 5}),
            ("A0001-PHE:CE1", "A0001-PHE:CD1", {"weight": 5}),
            ("A0001-PHE:CE2", "A0001-PHE:CD2", {"weight": 5}),
            ("A0001-PHE:CD1", "A0001-PHE:CG", {"weight": 5}),
            ("A0001-PHE:CD2", "A0001-PHE:CG", {"weight": 5}),
        ]

        pebbled_graph = nx.DiGraph([
            edge for n1, n2, d in edges for edge in ((n1, n2, d), (n2, n1, dict(d, weight=0)))
        ])

        pebbles = {
            "A0001-PHE:CZ": 1,
            "A0001-PHE:CE2": 1,
            "A0001-PHE:CD2": 1,
            "A0001-PHE:CE1": 1,
            "A0001-PHE:CD1": 1,
            "A0001-PHE:CG": 6
        }

        nx.set_node_attributes(pebbled_graph, pebbles, "pebbles")

        pebbled_graph.graph[GraphKey.STANDARD_EDGES.value] = [
            ("A0001-PHE:CZ", "A0001-PHE:CE2", 5)
        ]
        pebbled_graph.graph[GraphKey.COVALENT_EDGES.value] = []
        pebbled_graph.graph[GraphKey.NON_COVALENT_EDGES.value] = []
        pebbled_graph.graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value] = []

        yield pebbled_graph

    @pytest.fixture
    def mock_graph_gly_gly_peptide(self):
        edges = [
            ("A0001-GLY:H", "A0001-GLY:N", {"weight": 5}),
            ("A0001-GLY:N", "A0001-GLY:CA", {"weight": 5}),
            ("A0001-GLY:HA2", "A0001-GLY:CA", {"weight": 5}),
            ("A0001-GLY:HA3", "A0001-GLY:CA", {"weight": 5}),
            ("A0001-GLY:CA", "A0001-GLY:C", {"weight": 5}),
            ("A0001-GLY:O", "A0001-GLY:C", {"weight": 5}),
            #
            ("A0001-GLY:C", "A0002-GLY:N", {"weight": 5}),
            #
            ("A0002-GLY:H", "A0002-GLY:N", {"weight": 5}),
            ("A0002-GLY:N", "A0002-GLY:CA", {"weight": 5}),
            ("A0002-GLY:HA2", "A0002-GLY:CA", {"weight": 5}),
            ("A0002-GLY:HA3", "A0002-GLY:CA", {"weight": 5}),
            ("A0002-GLY:CA", "A0002-GLY:C", {"weight": 5}),
            ("A0002-GLY:O", "A0002-GLY:C", {"weight": 5}),
        ]

        pebbled_graph = nx.DiGraph([
            edge for n1, n2, d in edges for edge in ((n1, n2, d), (n2, n1, dict(d, weight=0)))
        ])

        pebbles = {
            "A0001-GLY:H": 1,
            "A0001-GLY:N": 1,
            "A0001-GLY:HA2": 1,
            "A0001-GLY:HA3": 1,
            "A0001-GLY:CA": 1,
            "A0001-GLY:O": 1,
            "A0001-GLY:C": 1,
            "A0002-GLY:H": 1,
            "A0002-GLY:N": 1,
            "A0002-GLY:HA2": 1,
            "A0002-GLY:HA3": 1,
            "A0002-GLY:CA": 1,
            "A0002-GLY:O": 1,
            "A0002-GLY:C": 6
        }

        nx.set_node_attributes(pebbled_graph, pebbles, "pebbles")

        pebbled_graph.graph[GraphKey.STANDARD_EDGES.value] = [
            ("A0001-GLY:C", "A0002-GLY:N", 1),
            ("A0001-GLY:O", "A0001-GLY:C", 1),
            ("A0002-GLY:O", "A0002-GLY:C", 1)
        ]
        pebbled_graph.graph[GraphKey.COVALENT_EDGES.value] = []
        pebbled_graph.graph[GraphKey.NON_COVALENT_EDGES.value] = []
        pebbled_graph.graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value] = []

        yield pebbled_graph


class TestMockProteinGraphs:

    @pytest.fixture
    def mock_protein_graph_gly_gly_hbond_5(self):
        """sp2 donor and sp2 acceptor"""
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_gly_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_acceptor(len(coords1) + 2, 2, 1.2, 1.5)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        protein_graph = construct_protein_graph_base(
            coords=coords1 + coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )

        add_missing_nodes(
            graphs=protein_graph.graphs,
            atom_df=protein_graph.hydrogen_df.loc[protein_graph.hydrogen_df.node_id == "A0001-GLY:H", :],
            clear_key_lists=False
        )
        add_bonds_from_frame(
            graphs=protein_graph.graphs,
            bond_frame=pd.DataFrame.from_records([{
                'node_1': "A0001-GLY:H", 'node_2': "A0001-GLY:N",
                'bond_type': InteractionType.COVALENT.value, 'bond_length': 1.07
            }]),
            bond_attributes=None,  # defaults to just length
            pebble_graph_key=GraphKey.COVALENT_EDGES.value,
            pebble_graph_weight=5,
            pebble_graph_quantified_keys=None
        )

        add_bonds_from_frame(
            graphs=protein_graph.graphs,
            bond_frame=pd.DataFrame.from_records([{
                'node_1': "A0001-GLY:H", 'node_2': "A0002-GLY:O",
                'bond_type': InteractionType.H_BOND.value, 'bond_length': 1.92, 'energy': -1.749
            }]),
            bond_attributes={"bond_length": 'bond_length', "energy": 'energy'},
            pebble_graph_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value,
            pebble_graph_weight=5,
            pebble_graph_quantified_keys=['energy']
        )

        # sort quantified edges
        protein_graph.graphs["pebble"].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value] = \
            sorted(protein_graph.graphs["pebble"].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value], key=lambda x: x[3])

        yield protein_graph

    @pytest.fixture
    def mock_protein_graph_arg_glu_salt_bridge_5(self):
        """sp2 donor and sp2 acceptor"""
        coords1, hm1, cov1, pe1, pm1, pb1 = graph_arg_donor(1, 1)
        coords2, hm2, cov2, pe2, pm2, pb2 = graph_glu_acceptor(len(coords1) + 2, 2, 1.6)
        others = [
            ["TER", f"{len(coords1) + 1:5d}", len(coords1) + 1]
        ]

        protein_graph = construct_protein_graph_base(
            coords=coords1 + coords2,
            others=others,
            hm=hm1 + hm2,
            cov=cov1 + cov2,
            pe=pe1 + pe2,
            pm=dict(pm1, **pm2),
            pb=pb1 + pb2
        )

        add_missing_nodes(
            graphs=protein_graph.graphs,
            atom_df=protein_graph.hydrogen_df,
            clear_key_lists=False
        )
        add_bonds_from_frame(
            graphs=protein_graph.graphs,
            bond_frame=pd.DataFrame.from_records([{
                    'node_1': row["node_id"], 'node_2': row["h_id"], 'bond_type': InteractionType.COVALENT.value,
                    'bond_length': row["length"]
                } for _, row in protein_graph.hydrogen_mapping.iterrows()
            ]),
            bond_attributes=None,  # defaults to just length
            pebble_graph_key=GraphKey.COVALENT_EDGES.value,
            pebble_graph_weight=5,
            pebble_graph_quantified_keys=None
        )

        add_bonds_from_frame(
            graphs=protein_graph.graphs,
            bond_frame=pd.DataFrame.from_records([{
                'node_1': "A0001-ARG:NH1", 'node_2': "A0002-GLU:OE1",
                'bond_type': InteractionType.SALT_BRIDGE.value, 'bond_length': 2.67, 'energy': -7.857
            },{
                'node_1': "A0001-ARG:NH2", 'node_2': "A0002-GLU:OE2",
                'bond_type': InteractionType.SALT_BRIDGE.value, 'bond_length': 2.67, 'energy': -7.857
            }]),
            bond_attributes={"bond_length": 'bond_length', "energy": 'energy'},
            pebble_graph_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value,
            pebble_graph_weight=5,
            pebble_graph_quantified_keys=['energy']
        )

        # sort quantified edges
        protein_graph.graphs["pebble"].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value] = \
            sorted(protein_graph.graphs["pebble"].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value], key=lambda x: x[3])

        yield protein_graph

    @pytest.fixture
    def mock_protein_graph_gly_gly_and_arg_glu_5(self):
        coords, hm, cov, pe, pm, pb = graph_gly_donor(1, 1)
        len_coords = len(coords)
        others = [
            ["TER", f"{len_coords + 1:5d}", len_coords + 1]
        ]

        coords2, hm2, cov2, pe2, pm2, pb2 = graph_gly_acceptor(len_coords + 2, 2, 1.2, 1.5)
        coords, hm, cov, pe, pm, pb = coords + coords2, hm + hm2, cov + cov2, pe + pe2, dict(pm, **pm2), pb + pb2
        len_coords += len(coords2) + 2
        others += [
            ["TER", f"{len_coords + 1:5d}", len_coords + 1]
        ]

        coords2, hm2, cov2, pe2, pm2, pb2 = graph_arg_donor(len_coords + 2, 3, 0.0, 6.0)
        coords, hm, cov, pe, pm, pb = coords + coords2, hm + hm2, cov + cov2, pe + pe2, dict(pm, **pm2), pb + pb2
        len_coords += len(coords2) + 2
        others += [
            ["TER", f"{len_coords + 1:5d}", len_coords + 1]
        ]

        coords2, hm2, cov2, pe2, pm2, pb2 = graph_glu_acceptor(len_coords + 2, 4, 1.6, 6.6)
        coords, hm, cov, pe, pm, pb = coords + coords2, hm + hm2, cov + cov2, pe + pe2, dict(pm, **pm2), pb + pb2

        protein_graph = construct_protein_graph_base(
            coords=coords,
            others=others,
            hm=hm,
            cov=cov,
            pe=pe,
            pm=pm,
            pb=pb
        )

        add_missing_nodes(
            graphs=protein_graph.graphs,
            atom_df=protein_graph.hydrogen_df,
            clear_key_lists=False
        )
        add_bonds_from_frame(
            graphs=protein_graph.graphs,
            bond_frame=pd.DataFrame.from_records([{
                'node_1': row["node_id"], 'node_2': row["h_id"], 'bond_type': InteractionType.COVALENT.value,
                'bond_length': row["length"]
            } for _, row in protein_graph.hydrogen_mapping.iterrows()
            ]),
            bond_attributes=None,  # defaults to just length
            pebble_graph_key=GraphKey.COVALENT_EDGES.value,
            pebble_graph_weight=5,
            pebble_graph_quantified_keys=None
        )

        add_bonds_from_frame(
            graphs=protein_graph.graphs,
            bond_frame=pd.DataFrame.from_records([{
                'node_1': "A0001-GLY:H", 'node_2': "A0002-GLY:O",
                'bond_type': InteractionType.H_BOND.value, 'bond_length': 1.92, 'energy': -1.749
            }, {
                'node_1': "A0003-ARG:NH1", 'node_2': "A0004-GLU:OE1",
                'bond_type': InteractionType.SALT_BRIDGE.value, 'bond_length': 2.67, 'energy': -7.857
            }, {
                'node_1': "A0003-ARG:NH2", 'node_2': "A0004-GLU:OE2",
                'bond_type': InteractionType.SALT_BRIDGE.value, 'bond_length': 2.67, 'energy': -7.857
            }]),
            bond_attributes={"bond_length": 'bond_length', "energy": 'energy'},
            pebble_graph_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value,
            pebble_graph_weight=5,
            pebble_graph_quantified_keys=['energy']
        )

        # sort quantified edges
        protein_graph.graphs["pebble"].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value] = \
            sorted(protein_graph.graphs["pebble"].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value], key=lambda x: x[3])

        yield protein_graph


#############
# Tests #####
#############

class TestProteinPebbleGame(TestMockGraphs, TestMockProteinGraphs, TestInvariants):

    # test basic graphs
    def test_game_pro_pentagon(self, mock_graph_pro_pentagon):
        spg = ProteinPebbleGame(pebble_graph=mock_graph_pro_pentagon)

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.STANDARD_EDGES.value)

        components = spg.get_components()

        assert len(components) == 1
        assert components[0]["size"] == 5
        assert components[0]["halo_size"] == 0

    def test_game_phe_hexagon(self, mock_graph_phe_hexagon):
        spg = ProteinPebbleGame(pebble_graph=mock_graph_phe_hexagon)

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.STANDARD_EDGES.value)

        components = spg.get_components()

        assert len(components) == 1
        assert components[0]["size"] == 6
        assert components[0]["halo_size"] == 0

    def test_game_gly_gly_peptide(self, mock_graph_gly_gly_peptide):
        spg = ProteinPebbleGame(pebble_graph=mock_graph_gly_gly_peptide)

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.STANDARD_EDGES.value)

        components = spg.get_components()

        assert len(components) == 2
        double_bond = first(x for x in components if x["size"] == 2)
        peptide_bond = first(x for x in components if x["size"] == 3)
        assert double_bond is not None and peptide_bond is not None
        assert double_bond["halo_size"] == 1  # CA atom
        assert peptide_bond["halo_size"] == 3  # 2 adjacent CA atoms and 1 H atom

    # test basic graphs with assertions for invariants
    def test_game_pro_pentagon_invariant(self, mock_graph_pro_pentagon):
        spg = ProteinPebbleGame(pebble_graph=mock_graph_pro_pentagon)

        self.assert_pebble_game_invariants(spg, mock_graph_pro_pentagon, GraphKey.STANDARD_EDGES.value)

        components = spg.get_components()

        assert len(components) == 1
        assert components[0]["size"] == 5
        assert components[0]["halo_size"] == 0

    def test_game_phe_hexagon_invariant(self, mock_graph_phe_hexagon):
        spg = ProteinPebbleGame(pebble_graph=mock_graph_phe_hexagon)

        self.assert_pebble_game_invariants(spg, mock_graph_phe_hexagon, GraphKey.STANDARD_EDGES.value)

        components = spg.get_components()

        assert len(components) == 1
        assert components[0]["size"] == 6
        assert components[0]["halo_size"] == 0

    def test_game_gly_gly_peptide_invariant(self, mock_graph_gly_gly_peptide):
        spg = ProteinPebbleGame(pebble_graph=mock_graph_gly_gly_peptide)

        self.assert_pebble_game_invariants(spg, mock_graph_gly_gly_peptide, GraphKey.STANDARD_EDGES.value)

        components = spg.get_components()

        assert len(components) == 2
        double_bond = first(x for x in components if x["size"] == 2)
        peptide_bond = first(x for x in components if x["size"] == 3)
        assert double_bond is not None and peptide_bond is not None
        assert double_bond["halo_size"] == 1  # CA atom
        assert peptide_bond["halo_size"] == 3  # 2 adjacent CA atoms and 1 H atom

    # test complex graphs
    def test_game_gly_gly_hbond_5(self, mock_protein_graph_gly_gly_hbond_5):
        protein_graph = mock_protein_graph_gly_gly_hbond_5

        spg = ProteinPebbleGame(pebble_graph=protein_graph.graphs['pebble'])

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.STANDARD_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 1
        assert components[0]["size"] == 2
        assert components[0]["halo_size"] == 1

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 1
        assert components[0]["size"] == 2
        assert components[0]["halo_size"] == 2
        assert "A0001-GLY:H" in components[0]["halo"]

    def test_game_arg_glu_salt_bridge_5(self, mock_protein_graph_arg_glu_salt_bridge_5):
        protein_graph = mock_protein_graph_arg_glu_salt_bridge_5

        spg = ProteinPebbleGame(pebble_graph=protein_graph.graphs['pebble'])

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.STANDARD_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 2
        for comp in components:
            assert comp["halo_size"] == 0

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.COVALENT_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 2
        arg_comp = first(comp for comp in components if comp["halo_size"] == 4)
        glu_comp = first(comp for comp in components if comp["halo_size"] == 0)
        assert arg_comp is not None
        assert glu_comp is not None
        assert arg_comp["size"] == 3
        assert glu_comp["size"] == 3

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 1
        assert components[0]["size"] == 6
        assert components[0]["halo_size"] == 4

    def test_game_gly_gly_and_arg_glu_5(self, mock_protein_graph_gly_gly_and_arg_glu_5):
        protein_graph = mock_protein_graph_gly_gly_and_arg_glu_5

        spg = ProteinPebbleGame(pebble_graph=protein_graph.graphs['pebble'])

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.STANDARD_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 3
        assert len(list(comp for comp in components if comp["halo_size"] == 0)) == 2
        assert len(list(comp for comp in components if comp["halo_size"] == 1)) == 1

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.COVALENT_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 3
        arg_comp = first(comp for comp in components if comp["halo_size"] == 4)
        glu_comp = first(comp for comp in components if comp["halo_size"] == 0)
        gly_comp = first(comp for comp in components if comp["halo_size"] == 1)
        assert arg_comp is not None
        assert glu_comp is not None
        assert gly_comp is not None
        assert arg_comp["size"] == 3
        assert glu_comp["size"] == 3
        assert gly_comp["size"] == 2

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 2
        salt_bridge = first(comp for comp in components if comp["size"] == 6)
        double_bond = first(comp for comp in components if comp["size"] == 2)
        assert salt_bridge is not None
        assert double_bond is not None
        assert salt_bridge["halo_size"] == 4
        assert double_bond["halo_size"] == 2
        assert "A0001-GLY:H" in double_bond["halo"]

    # test complex graphs with assertions for invariants
    def test_game_gly_gly_hbond_5_invariant(self, mock_protein_graph_gly_gly_hbond_5):
        protein_graph = mock_protein_graph_gly_gly_hbond_5

        spg = ProteinPebbleGame(pebble_graph=protein_graph.graphs['pebble'])

        self.assert_pebble_game_invariants(spg, protein_graph.graphs['pebble'], GraphKey.STANDARD_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 1
        assert components[0]["size"] == 2
        assert components[0]["halo_size"] == 1

        self.assert_pebble_game_invariants(spg, protein_graph.graphs['pebble'], GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 1
        assert components[0]["size"] == 2
        assert components[0]["halo_size"] == 2
        assert "A0001-GLY:H" in components[0]["halo"]

    def test_game_arg_glu_salt_bridge_5_invariant(self, mock_protein_graph_arg_glu_salt_bridge_5):
        protein_graph = mock_protein_graph_arg_glu_salt_bridge_5

        spg = ProteinPebbleGame(pebble_graph=protein_graph.graphs['pebble'])

        self.assert_pebble_game_invariants(spg, protein_graph.graphs['pebble'], GraphKey.STANDARD_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 2
        for comp in components:
            assert comp["halo_size"] == 0

        self.assert_pebble_game_invariants(spg, protein_graph.graphs['pebble'], GraphKey.COVALENT_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 2
        arg_comp = first(comp for comp in components if comp["halo_size"] == 4)
        glu_comp = first(comp for comp in components if comp["halo_size"] == 0)
        assert arg_comp is not None
        assert glu_comp is not None
        assert arg_comp["size"] == 3
        assert glu_comp["size"] == 3

        self.assert_pebble_game_invariants(spg, protein_graph.graphs['pebble'], GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 1
        assert components[0]["size"] == 6
        assert components[0]["halo_size"] == 4

    @pytest.mark.slow
    def test_game_gly_gly_and_arg_glu_5_invariant(self, mock_protein_graph_gly_gly_and_arg_glu_5):
        protein_graph = mock_protein_graph_gly_gly_and_arg_glu_5

        spg = ProteinPebbleGame(pebble_graph=protein_graph.graphs['pebble'])

        self.assert_pebble_game_invariants(spg, protein_graph.graphs['pebble'], GraphKey.STANDARD_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 3
        assert len(list(comp for comp in components if comp["halo_size"] == 0)) == 2
        assert len(list(comp for comp in components if comp["halo_size"] == 1)) == 1

        self.assert_pebble_game_invariants(spg, protein_graph.graphs['pebble'], GraphKey.COVALENT_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 3
        arg_comp = first(comp for comp in components if comp["halo_size"] == 4)
        glu_comp = first(comp for comp in components if comp["halo_size"] == 0)
        gly_comp = first(comp for comp in components if comp["halo_size"] == 1)
        assert arg_comp is not None
        assert glu_comp is not None
        assert gly_comp is not None
        assert arg_comp["size"] == 3
        assert glu_comp["size"] == 3
        assert gly_comp["size"] == 2

        self.assert_pebble_game_invariants(spg, protein_graph.graphs['pebble'], GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value)

        # assertions
        components = spg.get_components()

        assert len(components) == 2
        salt_bridge = first(comp for comp in components if comp["size"] == 6)
        double_bond = first(comp for comp in components if comp["size"] == 2)
        assert salt_bridge is not None
        assert double_bond is not None
        assert salt_bridge["halo_size"] == 4
        assert double_bond["halo_size"] == 2
        assert "A0001-GLY:H" in double_bond["halo"]


class TestDilutionPebbleGame(TestMockGraphs, TestMockProteinGraphs):

    def test_dilution_gly_gly_hbond_5(self, mock_protein_graph_gly_gly_hbond_5):
        protein_graph = mock_protein_graph_gly_gly_hbond_5

        spg = ProteinPebbleGame(pebble_graph=protein_graph.graphs['pebble'])

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.STANDARD_EDGES.value)

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.COVALENT_EDGES.value)

        # assertions
        base_components = spg.get_components()
        assert len(base_components) == 1
        assert base_components[0]["size"] == 2
        assert base_components[0]["halo_size"] == 1

        counter = 0
        for edge, components in spg.play_component_pebble_game_dilution(verbose=False,
                                                                        edge_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value):
            counter += 1
            if counter == 1:
                assert len(components) == 1
                assert components[0]["size"] == 2
                assert components[0]["halo_size"] == 2
                assert "A0001-GLY:H" in components[0]["halo"]

            else:
                pytest.fail("Too many dilution steps.")

    def test_dilution_arg_glu_salt_bridge_5(self, mock_protein_graph_arg_glu_salt_bridge_5):
        protein_graph = mock_protein_graph_arg_glu_salt_bridge_5

        spg = ProteinPebbleGame(pebble_graph=protein_graph.graphs['pebble'])

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.STANDARD_EDGES.value)

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.COVALENT_EDGES.value)

        # assertions
        base_components = spg.get_components()
        assert len(base_components) == 2
        arg_comp = first(comp for comp in base_components if comp["halo_size"] == 4)
        glu_comp = first(comp for comp in base_components if comp["halo_size"] == 0)
        assert arg_comp is not None
        assert glu_comp is not None
        assert arg_comp["size"] == 3
        assert glu_comp["size"] == 3

        counter = 0
        for edge, components in spg.play_component_pebble_game_dilution(verbose=False, edge_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value):
            counter += 1
            if counter == 1:
                assert edge[3] == -7.857
                assert len(components) == 1
                assert components[0]["size"] == 6
                assert components[0]["halo_size"] == 4

            else:
                pytest.fail("Too many dilution steps.")

    def test_dilution_gly_gly_and_arg_glu_5(self, mock_protein_graph_gly_gly_and_arg_glu_5):
        protein_graph = mock_protein_graph_gly_gly_and_arg_glu_5

        spg = ProteinPebbleGame(pebble_graph=protein_graph.graphs['pebble'])

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.STANDARD_EDGES.value)

        spg.play_component_pebble_game(verbose=False, edge_key=GraphKey.COVALENT_EDGES.value)

        # assertions
        base_components = spg.get_components()
        assert len(base_components) == 3
        arg_comp = first(comp for comp in base_components if comp["halo_size"] == 4)
        glu_comp = first(comp for comp in base_components if comp["halo_size"] == 0)
        gly_comp = first(comp for comp in base_components if comp["halo_size"] == 1)
        assert arg_comp is not None
        assert glu_comp is not None
        assert gly_comp is not None
        assert arg_comp["size"] == 3
        assert glu_comp["size"] == 3
        assert gly_comp["size"] == 2

        counter = 0
        for edge, components in spg.play_component_pebble_game_dilution(verbose=False, edge_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value):
            counter += 1
            if counter == 1:
                assert edge[3] == -7.857
                assert len(components) == 2
                salt_bridge = first(comp for comp in components if comp["size"] == 6)
                double_bond = first(comp for comp in components if comp["size"] == 2)
                assert salt_bridge is not None
                assert double_bond is not None
                assert salt_bridge["halo_size"] == 4
                assert double_bond["halo_size"] == 1

            elif counter == 2:
                assert edge[3] == -1.749
                assert len(components) == 2
                salt_bridge = first(comp for comp in components if comp["size"] == 6)
                double_bond = first(comp for comp in components if comp["size"] == 2)
                assert salt_bridge is not None
                assert double_bond is not None
                assert salt_bridge["halo_size"] == 4
                assert double_bond["halo_size"] == 2
                assert "A0001-GLY:H" in double_bond["halo"]

            else:
                pytest.fail("Too many dilution steps.")
