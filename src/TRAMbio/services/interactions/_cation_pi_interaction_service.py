import math
from typing import List

import pandas as pd
import numpy as np

from TRAMbio.services.parameter import ParameterRegistry as ParameterRegistry, lock_registry, CationPiInteractionParameter
from TRAMbio.services.interactions.registry import IInteractionService, InteractionServiceRegistry
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.functions.selection_functions import get_valid_cation_atoms, get_aromatic_rings
from TRAMbio.util.graphein.functions import compute_distmat
from TRAMbio.util.structure_library.graph_struct import ProteinGraph as ProteinGraph, GraphKey, add_bonds_from_frame
from TRAMbio.util.functions.numpy import angle_between

from loguru import logger as logger


__all__ = []


class CationPiInteractionService(IInteractionService):

    @property
    def name(self):
        return "CationPiInteractionService"

    @property
    def interaction_types(self) -> List[InteractionType]:
        return [InteractionType.CATION_PI]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        if not ParameterRegistry.get_parameter_set(parameter_id)(CationPiInteractionParameter.INCLUDE.value):
            return
        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id)
        cutoff_distance = parameter_registry(CationPiInteractionParameter.CUTOFF_DISTANCE.value)
        bar_count = parameter_registry(CationPiInteractionParameter.BAR_COUNT.value)

        pdb_df = protein_graph.heavy_atom_df.copy()

        # get suitable cation atoms
        cation_atoms = get_valid_cation_atoms(
            pdb_df=pdb_df,
            h_frame=protein_graph.hydrogen_mapping,
            columns=['node_id', 'x_coord', 'y_coord', 'z_coord']
        )
        cation_atoms = cation_atoms.loc[
                       cation_atoms['node_id'].apply(lambda node_id: len(list(protein_graph.graphs['atom'][node_id])) <= 1),
                       :
                       ].reset_index(drop=True)

        ring_centroids, ring_normals = get_aromatic_rings(pdb_df=pdb_df)

        coord_frame = pd.concat([ring_centroids, cation_atoms], axis=0, ignore_index=True)

        if len(coord_frame) < 2:
            if verbose:
                logger.debug(f"Found 0 cation-pi interactions")
            return

        distmat = compute_distmat(coord_frame)
        distmat = distmat[distmat > 0.4]
        t_distmat = distmat[distmat <= cutoff_distance]

        cation_pi_records = []

        inds = list(i for i in zip(*np.where(~np.isnan(t_distmat))) if i[0] < i[1])
        for i in inds:
            length = t_distmat[i[0]][i[1]]
            if isinstance(coord_frame['node_id'][i[0]], str) and isinstance(coord_frame['residue_id'][i[1]], str):
                index_ring = i[1]
                index_cation = i[0]
            elif isinstance(coord_frame['node_id'][i[1]], str) and isinstance(coord_frame['residue_id'][i[0]], str):
                index_ring = i[0]
                index_cation = i[1]
            else:
                # two rings or two cations
                continue

            ring_normal = ring_normals.loc[index_ring, ['x_normal', 'y_normal', 'z_normal']].to_numpy()
            ring_center = coord_frame.loc[index_ring, ['x_coord', 'y_coord', 'z_coord']].to_numpy()
            node_cation = coord_frame['node_id'][index_cation]
            cation_pos = coord_frame.loc[index_cation, ['x_coord', 'y_coord', 'z_coord']].to_numpy()
            dist_cation_ring = (ring_center - cation_pos)

            interaction_angle = angle_between(dist_cation_ring, ring_normal)
            # use supplement of the angle if > 90 degrees
            interaction_angle = min(math.fabs(interaction_angle), math.fabs(interaction_angle - 180))

            if interaction_angle >= 60.0:
                continue

            # insert bond between cation and first reference atom of ring (e.g., CG for PHE)
            cation_pi_records.append({
                'node_1': ring_normals['ref_node_0'][index_ring], 'node_2': node_cation,
                'bond_length': length, 'bond_type': InteractionType.CATION_PI.value, 'extra': [interaction_angle]
            })

        cation_pi_frame = pd.DataFrame.from_records(cation_pi_records)
        if len(cation_pi_frame) < 1:
            if verbose:
                logger.debug(f"Found 0 cation-pi interactions")
            return

        cond_1 = cation_pi_frame.groupby(by=['node_1'])['bond_length'].transform("min") == cation_pi_frame['bond_length']
        cond_2 = cation_pi_frame.groupby(by=['node_2'])['bond_length'].transform("min") == cation_pi_frame['bond_length']

        cation_pi_frame: pd.DataFrame = cation_pi_frame.loc[cond_1 & cond_2, :].reset_index(drop=True)

        new_edges = add_bonds_from_frame(
            graphs=protein_graph.graphs, bond_frame=cation_pi_frame,
            bond_attributes={"bond_length": 'bond_length', "extra": 'extra'},
            pebble_graph_key=GraphKey.NON_COVALENT_EDGES.value,
            pebble_graph_weight=bar_count,
            pebble_graph_quantified_keys=None
        )

        if verbose:
            if new_edges > 1:
                logger.info(f"Found {new_edges} cation-pi interactions")
            if new_edges == 1:
                logger.info(f"Found 1 cation-pi interaction")


InteractionServiceRegistry.NON_COV.register_service(CationPiInteractionService())
