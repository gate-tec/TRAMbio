from typing import List, Tuple

import pandas as pd
import numpy as np

from TRAMbio.util.constants.graph import RING_NORMAL_ATOMS


###########################
# Selection Functions #####
###########################

def get_valid_cation_atoms(
        h_frame: pd.DataFrame,
        pdb_df: pd.DataFrame,
        columns: List[str] = None
) -> "pd.DataFrame":
    if columns is None:
        columns = ['chain_id', 'residue_number', 'node_id']

    h_bond_count: pd.DataFrame = h_frame[['node_id', 'h_id']].groupby(by=['node_id']).count().reset_index()  # don't reset index here
    charged_atoms: pd.DataFrame = pdb_df[['node_id', 'charge']].loc[pdb_df['charge'].str.contains('+', regex=False), :].reset_index(drop=True)

    cation_nitrogens = h_bond_count.loc[
                       (h_bond_count['node_id'].str.contains(':N', regex=False)) &
                       ((h_bond_count['h_id'] == 3) | h_bond_count['node_id'].isin(charged_atoms['node_id'])), :]
    # if cation in ["LYS/NZ", "ARG/NH1 + NH2", "HIP/ND1 + NE2", "/NT1 + NT2 + NT"]
    cation_atoms: pd.DataFrame = pdb_df.loc[
        (pdb_df['node_id'].isin(cation_nitrogens['node_id'])) |
        ((pdb_df['residue_name'] == 'HIP') & (pdb_df['atom_name'].isin(['ND1', 'NE2']))) |
        (pdb_df['atom_name'].isin(['NH1', 'NH2', 'NT', 'NT1', 'NT2'])),
        columns]

    return cation_atoms


def get_aromatic_rings(
        pdb_df: pd.DataFrame
) -> "Tuple[pd.DataFrame, pd.DataFrame]":
    ring_atoms = (
        pd.DataFrame(RING_NORMAL_ATOMS)
        .unstack()
        .rename_axis(("residue_name", "atom_name"))
        .rename("ring_order")
    )

    # filter to ring atoms
    pdb_df = pdb_df.join(ring_atoms, on=["residue_name", "atom_name"]).dropna(subset=['ring_order'])
    pdb_df['residue_id'] = pdb_df.loc[:, 'node_id'].apply(lambda x: str(x)[:9])  # get residue_id from node_id

    # calculate ring centroids
    ring_centroids = pd.pivot_table(
        pdb_df, index=['residue_id'],
        values=['x_coord', 'y_coord', 'z_coord'],
        aggfunc={'x_coord': "mean", 'y_coord': "mean", 'z_coord': "mean"}
    ).reset_index()

    normal_records = []

    # calculate ring normals
    for residue_id in pdb_df['residue_id'].drop_duplicates().values:
        index_c0 = (pdb_df['residue_id'] == residue_id) & (pdb_df['ring_order'] == 0)
        index_c1 = (pdb_df['residue_id'] == residue_id) & (pdb_df['ring_order'] == 1)
        index_c2 = (pdb_df['residue_id'] == residue_id) & (pdb_df['ring_order'] == 2)
        c0_pos = pdb_df.loc[
            index_c0,
            ['x_coord', 'y_coord', 'z_coord']
        ].to_numpy()[0]
        c1_pos = pdb_df.loc[
            (pdb_df['residue_id'] == residue_id) & (pdb_df['ring_order'] == 1),
            ['x_coord', 'y_coord', 'z_coord']
        ].to_numpy()[0]
        c2_pos = pdb_df.loc[
            (pdb_df['residue_id'] == residue_id) & (pdb_df['ring_order'] == 2),
            ['x_coord', 'y_coord', 'z_coord']
        ].to_numpy()[0]

        dist_01 = c1_pos - c0_pos
        dist_02 = c2_pos - c0_pos

        ring_normal = np.cross(dist_01, dist_02)

        normal_records.append({
            'residue_id': residue_id, 'x_normal': ring_normal[0],
            'y_normal': ring_normal[1], 'z_normal': ring_normal[2],
            'ref_node_0': pdb_df.loc[index_c0, 'node_id'].values[0],
            'ref_node_1': pdb_df.loc[index_c1, 'node_id'].values[0],
            'ref_node_2': pdb_df.loc[index_c2, 'node_id'].values[0]
        })

    ring_normals = pd.DataFrame.from_records(normal_records)

    return ring_centroids, ring_normals
