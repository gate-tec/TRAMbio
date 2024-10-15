import math
from typing import List

import pandas as pd
import numpy as np

from TRAMbio.services.parameter import ParameterRegistry, lock_registry, AromaticInteractionParameter
from TRAMbio.services.interactions.registry import IInteractionService, InteractionServiceRegistry
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.functions.selection_functions import get_aromatic_rings
from TRAMbio.util.graphein.functions import compute_distmat
from TRAMbio.util.structure_library.graph_struct import ProteinGraph as ProteinGraph, GraphKey, add_bonds_from_frame
from TRAMbio.util.functions.numpy import angle_between

from loguru import logger as logger


__all__ = []


class AromaticInteractionService(IInteractionService):

    @property
    def name(self):
        return "AromaticInteractionService"

    @property
    def interaction_types(self) -> List[InteractionType]:
        return [InteractionType.PI_STACKING, InteractionType.T_STACKING]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        if not ParameterRegistry.get_parameter_set(parameter_id)(AromaticInteractionParameter.INCLUDE.value):
            return
        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id)
        angle_variance = parameter_registry(AromaticInteractionParameter.ANGLE_VARIANCE.value)
        cutoff_distance_pi = parameter_registry(AromaticInteractionParameter.CUTOFF_DISTANCE_PI.value)
        cutoff_distance_t = parameter_registry(AromaticInteractionParameter.CUTOFF_DISTANCE_T.value)
        bar_count = parameter_registry(AromaticInteractionParameter.BAR_COUNT.value)

        pdb_df = protein_graph.heavy_atom_df.copy()

        ring_centroids, ring_normals = get_aromatic_rings(pdb_df=pdb_df)

        if len(ring_centroids) < 2:
            if verbose:
                logger.debug("Found 0 aromatic interactions")
            return

        distmat = compute_distmat(ring_centroids)
        distmat = distmat[distmat > 0.4]
        t_distmat = distmat[distmat <= max(cutoff_distance_pi, cutoff_distance_t)]

        aromatic_records = []

        inds = list(i for i in zip(*np.where(~np.isnan(t_distmat))) if i[0] < i[1])
        for i in inds:
            length = t_distmat[i[0]][i[1]]
            ref_nodes_1 = [ring_normals[f'ref_node_{x}'][i[0]] for x in range(3)]
            ref_nodes_2 = [ring_normals[f'ref_node_{x}'][i[1]] for x in range(3)]
            bond_type = InteractionType.PI_STACKING.value

            ring_1_normal = ring_normals.loc[i[0], ['x_normal', 'y_normal', 'z_normal']].to_numpy()
            ring_1_center = ring_centroids.loc[i[0], ['x_coord', 'y_coord', 'z_coord']].to_numpy()
            ring_2_normal = ring_normals.loc[i[1], ['x_normal', 'y_normal', 'z_normal']].to_numpy()
            ring_2_center = ring_centroids.loc[i[1], ['x_coord', 'y_coord', 'z_coord']].to_numpy()

            dist_ring_centroids = (ring_2_center - ring_1_center)

            interaction_angle = angle_between(ring_1_normal, ring_2_normal)
            # use supplement of the angle if > 90 degrees
            interaction_angle = min(math.fabs(interaction_angle), math.fabs(interaction_angle - 180))

            if 60 - angle_variance < interaction_angle < 90 + angle_variance and length < cutoff_distance_t:
                # correct angle for t-stacking
                bond_type = InteractionType.T_STACKING.value
            elif interaction_angle >= 30 + angle_variance or length > cutoff_distance_pi:
                # angle too wide or distance too large for pi-stacking
                continue

            angle_1 = angle_between(dist_ring_centroids, ring_1_normal)
            angle_2 = angle_between(dist_ring_centroids, ring_2_normal)
            # use smaller supplement angle if > 90 degrees
            angle_1 = min(math.fabs(angle_1), math.fabs(angle_1 - 180))
            angle_2 = min(math.fabs(angle_2), math.fabs(angle_2 - 180))

            psi_angle = min(angle_1, angle_2)

            if psi_angle >= 45 + angle_variance:
                # angle too wide for interaction
                continue

            aromatic_records.append({
                'node_1': ref_nodes_1[0], 'node_2': ref_nodes_2[0],
                'bond_length': length, 'bond_type': bond_type,
                'extra': [interaction_angle, angle_1, angle_2]
            })

        aromatic_frame = pd.DataFrame.from_records(aromatic_records)
        if len(aromatic_frame) < 1:
            if verbose:
                logger.debug(f"Found 0 aromatic interactions")
            return

        cond_1 = aromatic_frame.groupby(by=['node_1'])['bond_length'].transform("min") == aromatic_frame['bond_length']
        cond_2 = aromatic_frame.groupby(by=['node_2'])['bond_length'].transform("min") == aromatic_frame['bond_length']

        aromatic_frame: pd.DataFrame = aromatic_frame.loc[cond_1 & cond_2, :].reset_index(drop=True)

        new_edges = add_bonds_from_frame(
            graphs=protein_graph.graphs, bond_frame=aromatic_frame,
            bond_attributes={"bond_length": 'bond_length', "extra": 'extra'},
            pebble_graph_key=GraphKey.NON_COVALENT_EDGES.value,
            pebble_graph_weight=bar_count,
            pebble_graph_quantified_keys=None
        )

        if verbose:
            if new_edges > 1:
                logger.info(f"Found {new_edges} aromatic interactions")
            elif new_edges == 1:
                logger.info(f"Found 1 aromatic interaction")


InteractionServiceRegistry.NON_COV.register_service(AromaticInteractionService())
