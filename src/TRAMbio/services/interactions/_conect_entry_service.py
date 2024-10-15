from typing import List

import pandas as pd

from TRAMbio.services.parameter import ParameterRegistry as ParameterRegistry, lock_registry, PdbEntryInteractionParameter
from TRAMbio.services.interactions.registry import IInteractionService, InteractionServiceRegistry
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.util.constants.interaction import InteractionType as InteractionType
from TRAMbio.util.structure_library.graph_struct import ProteinGraph as ProteinGraph, add_bonds_from_frame, \
    add_missing_nodes

from loguru import logger


__all__ = []


class ConectEntryService(IInteractionService):

    @property
    def name(self):
        return "ConectEntryService"

    @property
    def interaction_types(self) -> List[InteractionType]:
        return [InteractionType.CONECT]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        if not ParameterRegistry.get_parameter_set(parameter_id)(PdbEntryInteractionParameter.CONECT_INCLUDE.value):
            return

        conect_frame = protein_graph.others_df\
            .loc[protein_graph.others_df["record_name"] == 'CONECT', :].reset_index(drop=True)

        conect_edges = []
        missing_nodes = []
        for i, row in conect_frame.iterrows():
            entries = [int(row['entry'][x:x + 5].strip()) for x in range(0, len(row['entry'].rstrip()), 5)]
            atom_ids = [
                protein_graph.atom_df
                .loc[protein_graph.atom_df['atom_number'] == atom_number, 'node_id'].values
                for atom_number in entries
            ]
            atom_ids = [id_array[0] if len(id_array) > 0 else None for id_array in atom_ids]
            node_1 = atom_ids[0]
            if node_1 is None or not protein_graph.graphs['full'].has_node(node_1):
                if verbose:
                    logger.warning(f"Unable to find corresponding node with number {entries[0]} for CONECT entry:\nCONECT{row['entry']}")
                continue
            for j, node_2 in enumerate(atom_ids[1:]):
                if node_2 is None:
                    if verbose:
                        logger.warning(
                            f"Unable to find corresponding node with number {entries[j + 1]} for CONECT entry:\nCONECT{row['entry']}")
                    continue
                if not protein_graph.graphs['atom'].has_node(node_2):
                    missing_nodes.append(node_2)
                conect_edges.append({
                    'node_1': node_1, 'node_2': node_2, 'bond_type': InteractionType.CONECT.value
                })

        if len(missing_nodes) > 0:
            missing_atoms = protein_graph.atom_df.loc[protein_graph.atom_df['node_id'].isin(missing_nodes), :]
            # add relevant atoms as nodes
            add_missing_nodes(graphs=protein_graph.graphs, atom_df=missing_atoms, clear_key_lists=False)

        if len(conect_edges) > 0:
            new_edges = add_bonds_from_frame(
                graphs=protein_graph.graphs, bond_frame=pd.DataFrame.from_records(conect_edges), bond_attributes={}
            )

            if verbose:
                if new_edges > 1:
                    logger.info(f"Added {new_edges} new edges from CONECT entries")
                elif new_edges == 1:
                    logger.info(f"Added 1 new edge from CONECT entries")


InteractionServiceRegistry.COV.register_service(ConectEntryService())
