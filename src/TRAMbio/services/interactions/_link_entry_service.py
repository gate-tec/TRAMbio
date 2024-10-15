from typing import List

import pandas as pd

from TRAMbio.services.parameter import ParameterRegistry as ParameterRegistry, lock_registry, PdbEntryInteractionParameter
from TRAMbio.services.interactions.registry import IInteractionService, InteractionServiceRegistry
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.util.constants.interaction import InteractionType as InteractionType
from TRAMbio.util.constants.pdb import LINK_MAP
from TRAMbio.util.functions.pdb_functions import parse_entry_by_dict
from TRAMbio.util.structure_library.graph_struct import ProteinGraph as ProteinGraph, add_bonds_from_frame, \
    add_missing_nodes

from loguru import logger


__all__ = []


class LinkEntryService(IInteractionService):

    @property
    def name(self):
        return "LinkEntryService"

    @property
    def interaction_types(self) -> List[InteractionType]:
        return [InteractionType.LINK]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        if not ParameterRegistry.get_parameter_set(parameter_id)(PdbEntryInteractionParameter.LINK_INCLUDE.value):
            return

        link_df = protein_graph.others_df\
            .loc[protein_graph.others_df["record_name"] == 'LINK', :].reset_index(drop=True)
        atom_df = protein_graph.atom_df

        link_edges = []
        missing_nodes = []
        for i, row in link_df.iterrows():
            entry_dict = parse_entry_by_dict(row['entry'].rstrip(), LINK_MAP)
            if entry_dict is None:
                if verbose:
                    logger.warning(f"Unable to parse LINK entry:\nLINK  {row['entry']}")
                continue

            node_1_series: pd.Series = atom_df.loc[
                (atom_df['chain_id'] == entry_dict['chain_id_1']) &
                (atom_df['residue_number'] == entry_dict['residue_number_1']) &
                (atom_df['atom_name'] == entry_dict['atom_name_1']),
                'node_id'].reset_index(drop=True)
            node_2_series: pd.Series = atom_df.loc[
                (atom_df['chain_id'] == entry_dict['chain_id_2']) &
                (atom_df['residue_number'] == entry_dict['residue_number_2']) &
                (atom_df['atom_name'] == entry_dict['atom_name_2']),
                'node_id'].reset_index(drop=True)
            if len(node_1_series) != 1 or len(node_2_series) != 1:
                if verbose:
                    logger.warning(f"Unable to find corresponding nodes for LINK entry:\nLINK  {row['entry']}")
                continue
            node_1 = str(node_1_series.values[0])
            node_2 = str(node_2_series.values[0])

            if not protein_graph.graphs['atom'].has_node(node_1):
                missing_nodes.append(node_1)
            if not protein_graph.graphs['atom'].has_node(node_2):
                missing_nodes.append(node_2)

            link_edges.append({
                'node_1': node_1, 'node_2': node_2, 'bond_type': InteractionType.LINK.value
            })

        if len(missing_nodes) > 0:
            missing_atoms = protein_graph.atom_df.loc[protein_graph.atom_df['node_id'].isin(missing_nodes), :]
            # add relevant atoms as nodes
            add_missing_nodes(graphs=protein_graph.graphs, atom_df=missing_atoms, clear_key_lists=False)

        if len(link_edges) > 0:
            new_edges = add_bonds_from_frame(
                graphs=protein_graph.graphs, bond_frame=pd.DataFrame.from_records(link_edges), bond_attributes={}
            )
            if verbose:
                if new_edges > 1:
                    logger.info(f"Added {new_edges} new edges from LINK entries")
                elif new_edges == 1:
                    logger.info(f"Added 1 new edge from LINK entries")


InteractionServiceRegistry.COV.register_service(LinkEntryService())
