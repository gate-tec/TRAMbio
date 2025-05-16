from typing import Tuple, List, Generator

from TRAMbio.util.structure_library.components import PebbleGameResult
from TRAMbio.util.structure_library.graph_struct import ProteinGraph, GraphKey
from TRAMbio.pebble_game.protein_pebble_game import ProteinPebbleGame


__all__ = ["apply_pebble_game"]


def apply_pebble_game(
        protein_graph: ProteinGraph,
        verbose: bool = False
) -> Generator[Tuple[str, List[PebbleGameResult]], None, None]:
    spg = ProteinPebbleGame(pebble_graph=protein_graph.graphs['pebble'])

    spg.play_component_pebble_game(
        verbose=False, desc='Template edges', edge_key=GraphKey.STANDARD_EDGES.value)

    spg.play_component_pebble_game(
        verbose=verbose, desc='Covalent edges', edge_key=GraphKey.COVALENT_EDGES.value)

    spg.play_component_pebble_game(
        verbose=verbose, desc='Non-covalent edges', edge_key=GraphKey.NON_COVALENT_EDGES.value)

    base_components = spg.get_components()
    yield "-INF", base_components

    pg_iterator = spg.play_component_pebble_game_dilution(
        verbose=verbose, desc='Dilution edges', edge_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value)

    for edge, components in pg_iterator:
        # default order of edge is (u, v, weight, key)
        if len(edge) >= 4:
            try:
                key = str(round(edge[3], 3))
            except TypeError:
                key = "({},{})".format(edge[0], edge[1])
        else:
            key = "({},{})".format(edge[0], edge[1])

        yield key, components
