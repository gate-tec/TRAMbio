from typing import List

import pandas as pd
import numpy as np

from TRAMbio.services.parameter import ParameterRegistry as ParameterRegistry, lock_registry, DisulphideBridgeParameter
from TRAMbio.services.interactions.registry import IInteractionService, InteractionServiceRegistry
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.util.constants.graph import DISULPHIDE_RESIS, DISULPHIDE_ATOMS
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.graphein.functions import compute_distmat
from TRAMbio.util.structure_library.graph_struct import ProteinGraph as ProteinGraph, GraphKey, add_bonds_from_frame

from loguru import logger as logger


__all__ = []


class DisulphideBridgeService(IInteractionService):

    _EPSILON = 1e-5
    """Allowed distance variance in Angstroms"""

    @property
    def name(self):
        return "DisulphideBridgeService"

    @property
    def interaction_types(self) -> List[InteractionType]:
        return [InteractionType.SS_BOND]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        if not ParameterRegistry.get_parameter_set(parameter_id)(DisulphideBridgeParameter.INCLUDE.value):
            return
        cutoff_distance = ParameterRegistry.get_parameter_set(parameter_id)(DisulphideBridgeParameter.CUTOFF_DISTANCE.value)
        heavy_atom_df = protein_graph.heavy_atom_df.copy()

        # expected distance: 2.04 Angstroms
        valid_sulfur = heavy_atom_df.loc[
                       (heavy_atom_df['residue_name'].isin(DISULPHIDE_RESIS)) &
                       (heavy_atom_df['atom_name'].isin(DISULPHIDE_ATOMS)) &
                       (~heavy_atom_df['node_id'].isin(protein_graph.hydrogen_mapping['node_id'])),
                       :].reset_index(drop=True)

        if len(valid_sulfur) < 2:
            if verbose:
                logger.debug("Cannot add disulphide interaction with fewer than two suitable Cystein sulfurs.")
            return

        distmat = compute_distmat(valid_sulfur)
        distmat = distmat[distmat > 0.4]
        t_distmat = distmat[distmat <= cutoff_distance]

        disulphide_records = []

        inds = list(i for i in zip(*np.where(~np.isnan(t_distmat))) if i[0] < i[1])
        for i in inds:
            length = t_distmat[i[0]][i[1]]
            node_1 = valid_sulfur["node_id"][i[0]]
            node_2 = valid_sulfur["node_id"][i[1]]
            if not (2.04 - self._EPSILON <= length <= 2.04 + self._EPSILON) and verbose:
                logger.debug(
                    f"Potential SS-bond between {node_1} and {node_2} doesn't match ({length:.3f} instead of 2.04).")

            disulphide_records.append({
                'node_1': node_1, 'node_2': node_2, 'bond_length': length, 'bond_type': InteractionType.SS_BOND.value
            })

        disulphide_frame = pd.DataFrame.from_records(disulphide_records)
        probable_bonds = len(disulphide_frame)

        if probable_bonds < 1:
            if verbose:
                logger.debug(f"No valid disulphide interaction found within {cutoff_distance} Angstroms.")
            return

        cond_1 = disulphide_frame\
                     .groupby(by=['node_1'])['bond_length']\
                     .transform("min") == disulphide_frame['bond_length']
        cond_2 = disulphide_frame\
                     .groupby(by=['node_2'])['bond_length']\
                     .transform("min") == disulphide_frame['bond_length']

        disulphide_frame = disulphide_frame.loc[cond_1 & cond_2, :].reset_index(drop=True)

        new_edges = add_bonds_from_frame(
            graphs=protein_graph.graphs, bond_frame=disulphide_frame,
            bond_attributes=None,  # defaults to just length
            pebble_graph_key=GraphKey.COVALENT_EDGES.value,
            pebble_graph_weight=5,
            pebble_graph_quantified_keys=None
        )
        if verbose:
            if new_edges > 1:
                logger.info(f"Found {new_edges} disulphide interactions")
            elif new_edges == 1:
                logger.info(f"Found 1 disulphide interaction")


InteractionServiceRegistry.COV.register_service(DisulphideBridgeService())
