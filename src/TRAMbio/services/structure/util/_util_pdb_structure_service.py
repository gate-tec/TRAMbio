from typing import List

import pandas as pd
import numpy as np

from TRAMbio.util.constants.graph import STANDARD_RESIDUES, RESI_GRAPH_TEMPLATES, HETATM_GRAPH_TEMPLATES, \
    MAX_NUMBER_OF_BONDS, AA_RESIDUES, NA_RESIDUES
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.graphein.constants import COVALENT_RADII
from TRAMbio.util.graphein.functions import compute_distmat, assign_bond_states_to_dataframe, \
    assign_covalent_radii_to_dataframe
from TRAMbio.util.structure_library.graph_struct import GraphDictionary, initialize_graphs_from_dataframe, \
    process_node_in_template_graph, GraphKey

from loguru import logger


_EPSILON = 0.0001
"""Allowed distance variance in Angstroms"""


def create_base_graphs(heavy_atom_df: pd.DataFrame, all_atom_df: pd.DataFrame, ter_flags: List[str], verbose: bool = False):
    dist_mat = compute_distmat(all_atom_df)

    graphs: GraphDictionary = initialize_graphs_from_dataframe(atom_df=all_atom_df, heavy_atom_df=heavy_atom_df)

    peptide_bond_df: pd.DataFrame = heavy_atom_df.loc[heavy_atom_df['atom_name'].isin(['N', 'C']), :]
    phosphodiester_bond_df: pd.DataFrame = heavy_atom_df.loc[heavy_atom_df['atom_name'].isin(["O3'", 'P']), :]

    # include atomic edges from RESI_GRAPH_TEMPLATES
    queue = []
    prev_residue = None
    prev_code = None
    for i, row in all_atom_df.iterrows():
        node_id: str = row['node_id']
        if node_id.rfind(':') != 9:
            # sanity check on the formatting
            raise ValueError(f'Ill formatted node_id: {node_id}')
        current_residue = node_id[:9]
        current_residue_name = row['residue_name']

        if current_residue != prev_residue:
            if queue:
                # recheck previous residue
                if prev_code is None:
                    raise ValueError(f'Unable to determine reference for residue name "{prev_residue}".')
                # update heavy atom dataframe
                heavy_atom_df_corrected = False
                while queue:
                    j, prev_row = queue.pop(0)
                    if not heavy_atom_df_corrected:
                        heavy_atom_df.loc[
                            (heavy_atom_df['chain_id'] == prev_row['chain_id']) & (heavy_atom_df['residue_number'] == prev_row['residue_number']),
                            'residue_name'
                        ] = prev_code[1]
                    graphs['atom'].nodes[prev_row['node_id']]['residue_name'] = prev_code[1]
                    graphs['full'].nodes[prev_row['node_id']]['residue_name'] = prev_code[1]

            if prev_residue is not None and prev_residue not in ter_flags \
                    and prev_residue.startswith(row['chain_id']):
                # Check for valid peptide bond:
                # 1. no TER flag present
                # 2. same chain
                if prev_residue[-3:] in AA_RESIDUES and current_residue[-3:] in AA_RESIDUES:
                    # 3. both are standard residues allowing a regular peptide bond

                    # insert inter-residue edge ('atom', 'pebble')
                    _insert_peptide_bond(
                        graphs,
                        peptide_bond_df=peptide_bond_df,
                        residue_a=prev_residue,
                        residue_b=current_residue
                    )

                if prev_residue[-3:] in NA_RESIDUES and current_residue[-3:] in NA_RESIDUES:
                    # 3. both are standard residues allowing a regular phosphodiester bond

                    # insert inter-residue edge ('atom', 'pebble')
                    _insert_phosphodiester_bond(
                        graphs,
                        phosphodiester_bond_df=phosphodiester_bond_df,
                        residue_a=prev_residue,
                        residue_b=current_residue
                    )

            prev_residue = current_residue
            prev_code = None

        if 'HIS' in current_residue:
            # delay processing of residue until its code is determined
            if row['element_symbol'] != 'H':
                # skip H atoms
                queue.append((i, row))
            elif 'HD1' in node_id:
                # either HID or HIP
                if prev_code is None:
                    prev_code = (current_residue, 'HID')
                else:
                    prev_code = (current_residue, 'HIP')
            elif 'HE2' in node_id:
                # either HIE or HIP
                if prev_code is None:
                    prev_code = (current_residue, 'HIE')
                else:
                    prev_code = (current_residue, 'HIP')
        if 'CYS' in current_residue:
            # delay processing of residue until its code is determined
            prev_code = (current_residue, 'CSS')
            if row['element_symbol'] != 'H':
                # skip H atoms
                queue.append((i, row))
            elif 'HG' in node_id:
                prev_code = (current_residue, 'CYH')

        # skip H atoms
        if row['element_symbol'] == 'H':
            continue

        if current_residue_name in RESI_GRAPH_TEMPLATES.keys():
            residue_template = RESI_GRAPH_TEMPLATES[current_residue_name]
            if row['atom_name'] in residue_template.keys():
                # insert intra-residue edge ('atom', 'pebble')
                process_node_in_template_graph(
                    graphs,
                    residue_template[row['atom_name']],
                    node_id
                )
        elif current_residue_name in HETATM_GRAPH_TEMPLATES.keys():
            residue_template = HETATM_GRAPH_TEMPLATES[current_residue_name]
            if row['atom_name'] in residue_template.keys():
                # insert intra-residue edge ('atom', 'pebble')
                process_node_in_template_graph(
                    graphs,
                    residue_template[row['atom_name']],
                    node_id
                )
        else:
            if verbose:
                logger.debug(f"Unknown residue name {current_residue_name}")
            pass

    # double-check atomic edges with graphein-based calculation
    hydrogen_mapping: pd.DataFrame = compute_missing_atomic_edges(
        graphs=graphs, atom_df=all_atom_df, dist_mat=dist_mat, verbose=verbose
    )

    return graphs, hydrogen_mapping


def compute_missing_atomic_edges(graphs: GraphDictionary, atom_df: pd.DataFrame, dist_mat,
                                 tolerance: float = 0.56, verbose: bool = False):
    """
    Function for computing (probable) covalent edges and adding the ones not already present in the base graph.
    Implementation based on `graphein.protein.edges.atomic.add_atomic_edges`_.

     .. _graphein.protein.edges.atomic.add_atomic_edges:
        https://graphein.ai/modules/graphein.protein.html#graphein.protein.edges.atomic.add_atomic_edges
    """
    # We assign bond states to the dataframe, and then map these to covalent
    # radii
    pdb_df = atom_df.copy()

    # naively assign single bond state to all remaining unknown atoms
    extra_naive_bond_states = pd.Series(pdb_df["element_symbol"].map(
        {'C': 'Csb', 'O': 'Osb', 'N': 'Nsb', 'H': 'Hsb', 'S': 'Ssb', 'P': 'Psb'}
    ))

    pdb_df = assign_bond_states_to_dataframe(pdb_df)
    pdb_df = pdb_df.fillna(value={"atom_bond_state": extra_naive_bond_states})
    pdb_df = assign_covalent_radii_to_dataframe(pdb_df)

    # Create a covalent 'distance' matrix by adding the radius arrays with its
    # transpose
    covalent_radius_distance_matrix = np.add(
        np.array(pdb_df["covalent_radius"]).reshape(-1, 1),
        np.array(pdb_df["covalent_radius"]).reshape(1, -1),
    )

    # Add the tolerance
    covalent_radius_distance_matrix = (
            covalent_radius_distance_matrix + tolerance
    )

    # Threshold Distance Matrix to entries where the eucl distance is less than
    # the covalent radius plus tolerance and larger than 0.4
    dist_mat = dist_mat[dist_mat > 0.4]
    t_distmat = dist_mat[dist_mat < covalent_radius_distance_matrix]

    # recorded bonds to H-atoms
    h_records = {}
    potential_double_bonds = []

    # Get node IDs from non NaN entries in the thresholded distance matrix and
    # add the edge to the graph
    inds = list(i for i in zip(*np.where(~np.isnan(t_distmat))) if i[0] < i[1])
    for i in inds:
        length = t_distmat[i[0]][i[1]]
        node_1 = pdb_df["node_id"][i[0]]
        node_2 = pdb_df["node_id"][i[1]]
        chain_1 = pdb_df["chain_id"][i[0]]
        chain_2 = pdb_df["chain_id"][i[1]]

        # Check atoms are in the same chain
        if chain_1 != chain_2:
            continue

        # record bonded H-atoms
        if pdb_df["element_symbol"][i[0]] == 'H' or pdb_df["element_symbol"][i[1]] == 'H':
            index_h, index_n = i[0], i[1]
            node_h, node_n = node_1, node_2
            if pdb_df["element_symbol"][i[1]] == 'H':
                index_h, index_n = index_n, index_h
                node_h, node_n = node_n, node_h

            # if two atoms could be covalently bonded to this hydrogen
            # select the one with the lower tolerance usage
            if node_h in h_records.keys() and h_records[node_h][0] != node_n:
                covalent_radius_h = pdb_df['covalent_radius'][index_h]
                covalent_radius_new = pdb_df['covalent_radius'][index_n]
                covalent_radius_old = float(pdb_df.loc[pdb_df['node_id'] == h_records[node_h][0], 'covalent_radius'].to_numpy()[0])
                variance_new = length - covalent_radius_new - covalent_radius_h
                variance_old = h_records[node_h][1] - covalent_radius_old - covalent_radius_h
                if variance_old < variance_new:
                    continue
            h_records[node_h] = (node_n, length)
            continue

        # Check nodes are in graph (sanity check)
        if not (graphs['atom'].has_node(node_1) and graphs['atom'].has_node(node_2)):
            continue

        if graphs['atom'].has_edge(node_1, node_2):

            graphs['atom'].edges[node_1, node_2]["kind"].add(InteractionType.COVALENT.value)
            graphs['atom'].edges[node_1, node_2]["bond_length"] = length
            graphs['full'].edges[node_1, node_2]["kind"].add(InteractionType.COVALENT.value)
            graphs['full'].edges[node_1, node_2]["bond_length"] = length
        else:
            if verbose:
                logger.info(f'Missing covalent edge between "{node_1}" and "{node_2}"')
            # Check bond length for potential double bond
            element_1 = pdb_df['element_symbol'][i[0]]
            element_2 = pdb_df['element_symbol'][i[1]]
            bond_type = 'single'
            if ('ATOM' in pdb_df.record_name[i[0]] or 'ATOM' in pdb_df.record_name[i[1]]) and \
                    element_1 in ['C', 'N', 'O'] and element_2 in ['C', 'N', 'O']:
                expected_double = COVALENT_RADII[element_1 + 'db'] + COVALENT_RADII[element_2 + 'db']
                if length <= expected_double:
                    bond_type = 'double'

            if bond_type == 'double':
                potential_double_bonds.append((node_1, pdb_df['element_symbol'][i[0]], node_2, pdb_df['element_symbol'][i[1]], length))
            else:
                graphs['atom'].add_edge(node_1, node_2, kind={InteractionType.COVALENT.value}, bond_length=length)
                graphs['full'].add_edge(node_1, node_2, kind={InteractionType.COVALENT.value}, bond_length=length)
                graphs['pebble'].graph[GraphKey.STANDARD_EDGES.value].append((node_1, node_2, 5))

    # print(f"Detected {len(h_records)} bonded H-atoms")
    hydrogen_mapping = pd.DataFrame.from_records(
        list((h_id, *entry) for h_id, entry in h_records.items()),
        columns=['h_id', 'node_id', 'length']
    )
    for _, row in hydrogen_mapping.iterrows():
        graphs['full'].add_edge(row['node_id'], row['h_id'], kind={InteractionType.COVALENT.value}, bond_length=row['length'])

    # evaluate candidates for double bonds by comparing against expected number of (covalent) bonds
    hydrogen_counts = hydrogen_mapping.value_counts()
    for node_1, element_1, node_2, element_2, length in potential_double_bonds:
        num_bonded_h_1 = hydrogen_counts.get(node_1, 0)
        num_bonded_h_2 = hydrogen_counts.get(node_2, 0)
        if element_1 in MAX_NUMBER_OF_BONDS.keys() and element_2 in MAX_NUMBER_OF_BONDS.keys() and \
                MAX_NUMBER_OF_BONDS[element_1] - graphs['atom'].degree(node_1) - num_bonded_h_1 >= 2 and \
                MAX_NUMBER_OF_BONDS[element_2] - graphs['atom'].degree(node_2) - num_bonded_h_2 >= 2:
            graphs['pebble'].graph[GraphKey.STANDARD_EDGES.value].append((node_1, node_2, 6))
            if verbose:
                logger.debug(f"Detected double bond between {node_1} and {node_2}.")
        else:
            graphs['pebble'].graph[GraphKey.STANDARD_EDGES.value].append((node_1, node_2, 5))
        graphs['atom'].add_edge(node_1, node_2, kind={InteractionType.COVALENT.value}, bond_length=length)
        graphs['full'].add_edge(node_1, node_2, kind={InteractionType.COVALENT.value}, bond_length=length)

    return hydrogen_mapping


def _insert_peptide_bond(
        graphs: GraphDictionary, peptide_bond_df: pd.DataFrame,
        residue_a: str, residue_b: str,
        tolerance: float = 0.56
):
    target_node_id_a = residue_a + ":C"
    target_node_id_b = residue_b + ":N"

    target_frame: pd.DataFrame = peptide_bond_df[peptide_bond_df['node_id'].isin([target_node_id_a, target_node_id_b])]
    if len(target_frame) != 2:
        raise KeyError(f'Missing nodes for peptide bond between "{target_node_id_a}" and "{target_node_id_b}"')

    target_frame = target_frame.loc[:, ['residue_number', 'node_id', 'x_coord', 'y_coord', 'z_coord']].reset_index(
        drop=True)

    sequence_distance = int(target_frame['residue_number'][1]) - int(target_frame['residue_number'][0])

    pos_node_a = target_frame.loc[0, ['x_coord', 'y_coord', 'z_coord']].to_numpy()
    pos_node_b = target_frame.loc[1, ['x_coord', 'y_coord', 'z_coord']].to_numpy()

    dist_peptide_bond = pos_node_b - pos_node_a
    len_peptide_bond = np.linalg.norm(dist_peptide_bond)

    if sequence_distance == 1 or len_peptide_bond < 0.67 + 0.70 + tolerance:
        graphs['atom'].add_edge(target_node_id_a, target_node_id_b, kind={InteractionType.PEPTIDE_BOND.value}, base=True, bond_length=len_peptide_bond)
        graphs['full'].add_edge(target_node_id_a, target_node_id_b, kind={InteractionType.PEPTIDE_BOND.value}, base=True, bond_length=len_peptide_bond)

        graphs['pebble'].add_edge(target_node_id_a, target_node_id_b, weight=5)
        graphs['pebble'].add_edge(target_node_id_b, target_node_id_a, weight=0)
        graphs['pebble'].nodes[target_node_id_a]['pebbles'] = 1
        # Add 6th edge for rigid resonance/ peptide bond
        graphs['pebble'].graph[GraphKey.STANDARD_EDGES.value].append((target_node_id_a, target_node_id_b, 1))


def _insert_phosphodiester_bond(
        graphs: GraphDictionary, phosphodiester_bond_df: pd.DataFrame,
        residue_a: str, residue_b: str,
        tolerance: float = 0.56
):
    target_node_id_a = residue_a + ":O3'"
    target_node_id_b = residue_b + ":P"

    target_frame: pd.DataFrame = phosphodiester_bond_df[phosphodiester_bond_df['node_id'].isin([target_node_id_a, target_node_id_b])]
    if len(target_frame) != 2:
        raise KeyError(f'Missing nodes for phosphodiester bond between "{target_node_id_a}" and "{target_node_id_b}"')

    target_frame = target_frame.loc[:, ['residue_number', 'node_id', 'x_coord', 'y_coord', 'z_coord']].reset_index(
        drop=True)

    sequence_distance = int(target_frame['residue_number'][1]) - int(target_frame['residue_number'][0])

    pos_node_a = target_frame.loc[0, ['x_coord', 'y_coord', 'z_coord']].to_numpy()
    pos_node_b = target_frame.loc[1, ['x_coord', 'y_coord', 'z_coord']].to_numpy()

    dist_phosphodiester_bond = pos_node_b - pos_node_a
    len_phosphodiester_bond = np.linalg.norm(dist_phosphodiester_bond)

    if sequence_distance == 1 or len_phosphodiester_bond < 0.67 + 1.06 + tolerance:
        graphs['atom'].add_edge(target_node_id_a, target_node_id_b, kind={InteractionType.PHOSPHODIESTER_BOND.value},
                                base=True, bond_length=len_phosphodiester_bond)
        graphs['full'].add_edge(target_node_id_a, target_node_id_b, kind={InteractionType.PHOSPHODIESTER_BOND.value},
                                base=True, bond_length=len_phosphodiester_bond)

        graphs['pebble'].add_edge(target_node_id_a, target_node_id_b, weight=5)
        graphs['pebble'].add_edge(target_node_id_b, target_node_id_a, weight=0)
        graphs['pebble'].nodes[target_node_id_a]['pebbles'] = 1


def export_ter_flagged_residues(atom_df: pd.DataFrame, others_df: pd.DataFrame):
    ter_df = others_df.loc[others_df["record_name"] == 'TER', :].reset_index(drop=True)

    ter_flags = []
    for i, row in ter_df.iterrows():
        line_idx = row['line_idx']
        flagged_residue = atom_df.loc[
            atom_df['line_idx'] == line_idx - 1,
            ["chain_id", "residue_number", "insertion", "residue_name"]
        ].reset_index(drop=True)
        flag_id = "".join([
            str(flagged_residue.loc[0, "chain_id"]),
            str(flagged_residue.loc[0, "residue_number"]).rjust(4, '0'),
            str(flagged_residue.loc[0, "insertion"]) if len(str(flagged_residue.loc[0, "insertion"])) > 0 else '-',
            str(flagged_residue.loc[0, "residue_name"])
        ])
        ter_flags.append(flag_id)

    return ter_flags
