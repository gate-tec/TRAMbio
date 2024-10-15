from typing import List

import pandas as pd

from TRAMbio.services.parameter import ParameterRegistry as ParameterRegistry, lock_registry, PdbEntryInteractionParameter
from TRAMbio.services.interactions.registry import IInteractionService, InteractionServiceRegistry
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.util.constants.interaction import InteractionType as InteractionType
from TRAMbio.util.constants.pdb import SS_BOND_MAP
from TRAMbio.util.functions.pdb_functions import parse_entry_by_dict
from TRAMbio.util.structure_library.graph_struct import ProteinGraph as ProteinGraph, add_bonds_from_frame

from loguru import logger


__all__ = []


class SSBondEntryService(IInteractionService):

    @property
    def name(self):
        return "SSBondEntryService"

    @property
    def interaction_types(self) -> List[InteractionType]:
        return [InteractionType.SS_BOND]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        if not ParameterRegistry.get_parameter_set(parameter_id)(PdbEntryInteractionParameter.SSBOND_INCLUDE.value):
            return

        ss_bond_df = protein_graph.others_df\
            .loc[protein_graph.others_df["record_name"] == 'SSBOND', :].reset_index(drop=True)
        atom_df = protein_graph.atom_df

        ss_edges = []
        sg_index = atom_df['atom_name'] == 'SG'
        for i, row in ss_bond_df.iterrows():
            entry_dict = parse_entry_by_dict(row['entry'].rstrip(), SS_BOND_MAP)
            if entry_dict is None:
                if verbose:
                    logger.warning(f"Unable to parse SSBOND entry:\nSSBOND{row['entry']}")
                continue

            node_1_series: pd.Series = atom_df.loc[
                (atom_df['chain_id'] == entry_dict['chain_id_1']) &
                (atom_df['residue_number'] == entry_dict['residue_number_1']) &
                sg_index,
                'node_id'].reset_index(drop=True)
            node_2_series: pd.Series = atom_df.loc[
                (atom_df['chain_id'] == entry_dict['chain_id_2']) &
                (atom_df['residue_number'] == entry_dict['residue_number_2']) &
                sg_index,
                'node_id'].reset_index(drop=True)
            if len(node_1_series) != 1 or len(node_2_series) != 1:
                if verbose:
                    logger.warning(f"Unable to find corresponding nodes for SSBOND entry {entry_dict['ss_entry_number']}")
                continue
            node_1 = str(node_1_series.values[0])
            node_2 = str(node_2_series.values[0])

            ss_edges.append({
                'node_1': node_1, 'node_2': node_2, 'bond_type': InteractionType.SS_BOND.value
            })

        if len(ss_edges) > 0:
            add_bonds_from_frame(
                graphs=protein_graph.graphs, bond_frame=pd.DataFrame.from_records(ss_edges), bond_attributes={}
            )
            if verbose:
                logger.info(f"Added {len(ss_edges)} disulphide interactions from SSBOND entries")


InteractionServiceRegistry.COV.register_service(SSBondEntryService())
