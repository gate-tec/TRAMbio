from loguru import logger
import pandas as pd
import numpy as np

from typing import Optional

from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.functions.graph_neighborhood import NeighborType, AtomNeighborhood
from TRAMbio.util.graphein.constants import VDW_RADII
from TRAMbio.util.graphein.functions import compute_distmat
from TRAMbio.util.structure_library.graph_struct import GraphDictionary


__all__ = ["calculate_hydrophobic_interactions"]


def calculate_hydrophobic_interactions(
        graphs: GraphDictionary,
        heavy_atom_df: pd.DataFrame,
        cutoff_distance: float = 0.25,
        minimal_length: bool = True,
        verbose: bool = False
) -> "Optional[pd.DataFrame]":
    # a hydrophobic interaction is
    # placed for every carbon–carbon, carbon–sulfur, or
    # sulfur–sulfur pair, when their van der Waals surfaces are
    # within a cutoff distance of 0.25 Å.

    index_1 = heavy_atom_df['element_symbol'].isin(['C', 'S'])
    index_2 = heavy_atom_df['residue_name'].isin(
        ['ALA', 'CYS', 'CYH', 'CSS', 'PHE', 'GLY', 'ILE', 'LEU', 'MET', 'PRO', 'VAL', 'TRP']
    )
    # Graphein listing
    # index_2 = heavy_atom_df['residue_name'].isin(
    #     ['ALA', 'PHE', 'ILE', 'LEU', 'MET', 'PRO', 'VAL', 'TRP', 'TYR']
    # )

    pdb_df = heavy_atom_df.copy().loc[index_1 & index_2].reset_index(drop=True)
    pdb_df['vdw_radius'] = pd.Series(pdb_df["element_symbol"].map(VDW_RADII))

    if len(pdb_df) < 2:
        if verbose:
            logger.debug(f"Found 0 hydrophobic interactions")
        return None

    # compute eucl distances between valid atoms
    dist_mat = compute_distmat(pdb_df)

    # Create a van der Waals 'distance' matrix by adding the radius arrays with its transpose
    vdw_radius_distance_matrix = np.add(
        np.array(pdb_df["vdw_radius"]).reshape(-1, 1),
        np.array(pdb_df["vdw_radius"]).reshape(1, -1),
    )

    # Add the tolerance
    vdw_radius_distance_matrix = (
            vdw_radius_distance_matrix + cutoff_distance
    )

    # Threshold Distance Matrix to entries where the eucl distance is less than or equal to
    # the van der Waals radius plus tolerance and larger than or equal to 2.2 (in order to exclude covalent bonds)
    dist_mat = dist_mat[dist_mat >= 2.2]
    t_distmat = dist_mat[dist_mat <= vdw_radius_distance_matrix]

    hydrophobic_records = []

    inds = list(i for i in zip(*np.where(~np.isnan(t_distmat))) if i[0] < i[1])
    with AtomNeighborhood(graphs['full']) as neighborhood:
        for i in inds:
            length = t_distmat[i[0]][i[1]]
            node_1 = pdb_df["node_id"][i[0]]
            node_2 = pdb_df["node_id"][i[1]]

            if neighborhood.distance_type(node_1, node_2) not in [NeighborType.P_15, NeighborType.UNBOUNDED]:
                continue

            hydrophobic_records.append({
                'node_1': node_1, 'node_2': node_2, 'bond_length': length, 'bond_type': InteractionType.HYDROPHOBIC.value
            })

    hydrophobic_frame = pd.DataFrame.from_records(hydrophobic_records)

    if len(hydrophobic_frame) < 1:
        if verbose:
            logger.debug(f"Found 0 hydrophobic interactions within {cutoff_distance} Angstroms.")
        return None

    if minimal_length:
        hydrophobic_frame = hydrophobic_frame.sort_values(by=['bond_length'], ignore_index=True).reset_index(drop=True)

        minimal_indices = []
        node_ids = []
        for i, row in hydrophobic_frame.iterrows():
            node_1, node_2 = row["node_1"], row["node_2"]
            if node_1 not in node_ids and node_2 not in node_ids:
                node_ids.append(node_1)
                node_ids.append(node_2)
                minimal_indices.append(i)

        hydrophobic_frame = hydrophobic_frame.loc[minimal_indices, :].reset_index(drop=True)

    return hydrophobic_frame
