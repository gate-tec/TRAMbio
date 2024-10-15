from typing import Tuple, List

from TRAMbio import set_log_level
from TRAMbio.util.structure_library.components import PebbleGameResult
from TRAMbio.util.structure_library.graph_struct import ProteinGraph, GraphKey
from TRAMbio.pebble_game.protein_pebble_game import ProteinPebbleGame

import TRAMbio.services.parameter.registry as registry_space

from TRAMbio.services import StructureServiceRegistry

__all__ = ["init_mp_pool", "create_graph_from_frame", "apply_pebble_game_mp"]

__base_protein_graphs = None


def init_mp_pool(base_protein_graph, origin_registry, parameter_id, log_level):
    global __base_protein_graphs
    if __base_protein_graphs is None:
        __base_protein_graphs = {}
    __base_protein_graphs[parameter_id] = base_protein_graph

    # Since the new processes are spawned, the parameters need to be re-registered
    registry_space.ParameterRegistry.load_registry(origin_registry)

    set_log_level(log_level)


def create_graph_from_frame(args) -> Tuple[int, ProteinGraph]:
    (frame, raw_df), parameter_id = args
    global __base_protein_graphs

    if not isinstance(__base_protein_graphs, dict):
        raise ValueError("Base protein graphs were not initialized correctly.")
    base_protein_graph = __base_protein_graphs[parameter_id]

    pdb_service = StructureServiceRegistry.PDB.single_service()

    atom_df = pdb_service.export_atom_df(raw_df=raw_df, check_ids=False, parameter_id=parameter_id)
    others_df = pdb_service.export_others_df(raw_df=raw_df, ter_only=False, parameter_id=parameter_id)

    protein_graph = pdb_service.copy_graph_for_frame(atom_df=atom_df, others_df=others_df, protein_graph=base_protein_graph, parameter_id=parameter_id)

    pdb_service.apply_non_covalent_interactions(protein_graph=protein_graph, parameter_id=parameter_id)

    protein_graph.graphs['atom'].graph['frame'] = frame

    return frame, protein_graph


def apply_pebble_game_mp(args) -> Tuple[int, List[PebbleGameResult]]:
    frame, graphs = args
    spg = ProteinPebbleGame(pebble_graph=graphs['pebble'])

    spg.play_component_pebble_game(
        verbose=False, desc=f'', edge_key=GraphKey.STANDARD_EDGES.value)

    spg.play_component_pebble_game(
        verbose=False, desc=f'', edge_key=GraphKey.COVALENT_EDGES.value)

    spg.play_component_pebble_game(
        verbose=False, desc=f'', edge_key=GraphKey.NON_COVALENT_EDGES.value)

    spg.play_component_pebble_game(
        verbose=False, desc=f'', edge_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value)

    components = spg.get_components()

    return frame, components
