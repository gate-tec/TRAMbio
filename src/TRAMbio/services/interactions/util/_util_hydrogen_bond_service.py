import sys
from looseversion import LooseVersion
from typing import Optional, Tuple, List

if LooseVersion(sys.version) < LooseVersion("3.11"):
    from typing_extensions import TypedDict, Required, NotRequired
else:
    from typing import TypedDict, Required, NotRequired

import pandas as pd
import numpy as np

from loguru import logger
from tqdm import tqdm

from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.constants.graph import RESIDUE_HBOND_DONORS, RESIDUE_HBOND_ACCEPTORS, \
    DEFAULT_HYBRIDIZATION, RESIDUE_ATOM_HYBRIDIZATION
from TRAMbio.util.functions.selection_functions import get_valid_cation_atoms
from TRAMbio.util.graphein.functions import compute_distmat
from TRAMbio.util.structure_library.graph_struct import GraphDictionary

from TRAMbio.util.functions.numpy import angle_between
_tqdm_logger = logger.bind(task="tqdm")


__all__ = ["calculate_hydrogen_and_salt_bridge_bonds"]


##########################
# Labeling Functions #####
##########################

def _label_h_bond_donors_acceptors(pdb_df: pd.DataFrame, h_frame: pd.DataFrame):
    # Main chain NHs are default donors
    naive_donors = pd.Series(pdb_df['atom_name'].map(
        lambda x: 1 if x == 'N' else (3 if x in ['NT', 'NT1', 'NT2'] else 0))
    )

    h_bond_count: pd.DataFrame = h_frame[['node_id', 'h_id']].groupby(by=['node_id']).count().reset_index()
    suitable_nitrogens = h_bond_count.loc[
                       (h_bond_count['node_id'].str.contains(':N', regex=False)) & (h_bond_count['h_id'] >= 1), :
                         ].copy().set_index('node_id').to_dict('index')
    suitable_nitrogens_map = {key: entry['h_id'] for key, entry in suitable_nitrogens.items()}

    # Create series of known donor states
    ss = (
        pd.DataFrame(RESIDUE_HBOND_DONORS)
        .unstack()
        .rename_axis(("residue_name", "atom_name"))
        .rename("hbond_donor")
    )

    # Map known states to the DataFrame based on the residue and atom name
    pdb_df = pdb_df.join(ss, on=["residue_name", "atom_name"])

    # Fill the NaNs for suitable nitrogens known hydrogen atoms
    pdb_df = pdb_df.fillna(value={"hbond_donor": suitable_nitrogens_map})

    # Fill the NaNs with the naive states
    pdb_df = pdb_df.fillna(value={"hbond_donor": naive_donors})

    # Main chain O=Cs are default acceptors
    naive_acceptors = pd.Series(pdb_df['atom_name'].map(lambda x: 2 if x in ['O', 'OXT'] else 0))

    # Create series of known acceptor states
    ss = (
        pd.DataFrame(RESIDUE_HBOND_ACCEPTORS)
        .unstack()
        .rename_axis(("residue_name", "atom_name"))
        .rename("hbond_acceptor")
    )

    # Map known states to the DataFrame based on the residue and atom name
    pdb_df = pdb_df.join(ss, on=["residue_name", "atom_name"])

    # Fill the NaNs with the naive states
    pdb_df = pdb_df.fillna(value={"hbond_acceptor": naive_acceptors})

    return pdb_df


def _label_hybridization(pdb_df: pd.DataFrame):
    # Map atoms to naive states
    hetatm_hybridization = pd.Series(pdb_df["record_name"].map({'HETATM': 3}))
    naive_hybridization = pd.Series(pdb_df["atom_name"].map(DEFAULT_HYBRIDIZATION)).fillna(value=3)

    # Create series of known hybridization states
    ss = (
        pd.DataFrame(RESIDUE_ATOM_HYBRIDIZATION)
        .unstack()
        .rename_axis(("residue_name", "atom_name"))
        .rename("hybridization")
    )

    # Map known states to the DataFrame based on the residue and atom name
    pdb_df = pdb_df.join(ss, on=["residue_name", "atom_name"])

    # Fill the NaNs with the naive states
    pdb_df = pdb_df.fillna(value={"hybridization": hetatm_hybridization})
    pdb_df = pdb_df.fillna(value={"hybridization": naive_hybridization})

    return pdb_df


#####################################
# Dataframe selection functions #####
#####################################


def _get_valid_donor_acceptor_atoms(
        pdb_df: pd.DataFrame,
        h_frame: pd.DataFrame
) -> "tuple[pd.DataFrame, pd.DataFrame]":
    valid_donors = pdb_df.loc[(pdb_df['hbond_donor'] > 0) & (pdb_df['node_id'].isin(h_frame['node_id'])), :] \
        .reset_index(drop=True)
    valid_acceptors = pdb_df.loc[pdb_df['hbond_acceptor'] > 0, :].reset_index(drop=True)

    return valid_donors, valid_acceptors


def _get_valid_cation_anion_atoms(
        h_frame: pd.DataFrame,
        valid_donors: pd.DataFrame, valid_acceptors: pd.DataFrame
) -> "tuple[pd.DataFrame, pd.DataFrame]":
    # if cation in ["LYS/NZ", "ARG/NH1 + NH2", "HIP/ND1 + NE2", "/NT1 + NT2 + NT"]
    cation_atoms: pd.DataFrame = get_valid_cation_atoms(h_frame=h_frame, pdb_df=valid_donors)

    # if anion in ["ASP/OD1 + OD2", "GLU/OE1 + OE2", "/O + OXT"]
    oxt_frame = list(
        valid_acceptors.loc[valid_acceptors['atom_name'] == 'OXT', ['chain_id', 'residue_number']].apply(tuple, 1))
    anion_atoms: pd.DataFrame = valid_acceptors.loc[
        ((valid_acceptors['residue_name'].isin(['ASP', 'GLU'])) & (
            valid_acceptors['atom_name'].isin(['OD1', 'OD2', 'OE1', 'OE2']))) |
        ((valid_acceptors[['chain_id', 'residue_number']].apply(tuple, 1).isin(oxt_frame)) & (
            valid_acceptors['atom_name'].isin(['O', 'OXT']))) |
        (valid_acceptors['charge'].str.contains('-', regex=False)),
        ['chain_id', 'residue_number', 'node_id']]

    return cation_atoms, anion_atoms


############################
# Mayo-Energy Function #####
############################

def _calculate_gamma_angle(
        graphs: GraphDictionary,
        node_donor: str, node_accep: str,
        donor_pos: np.ndarray, acceptor_pos: np.ndarray,
        dh_dist: np.ndarray,
        hybridization_df: pd.DataFrame
):
    # for donor use plane containing selected H-atom and base-atom
    try:
        donor_base_atom = next(iter(graphs['atom'][node_donor]))
    except StopIteration:
        raise ValueError(f"Insufficient number of atoms to calculate plane of sp2-center at {node_donor}")
    # since no H-atoms are present in graphs['atom'], the base atom is always a distinct third atom
    donor_base_pos = graphs['atom'].nodes[donor_base_atom]["coords"]
    donor_base_dist = donor_base_pos - donor_pos

    donor_plane_normal = np.cross(dh_dist, donor_base_dist)

    # for acceptor use plane containing base atom (and possibly dist-2-atoms)
    acceptor_base_atoms = list(iter(graphs['atom'][node_accep]))
    acceptor_plane_normal = None

    if len(acceptor_base_atoms) == 0:
        raise ValueError(f"No neighbor found for sp2-acceptor {node_accep}")
    if len(acceptor_base_atoms) == 1 and hybridization_df.loc[hybridization_df['node_id'] == acceptor_base_atoms[0], ['hybridization']].to_numpy()[0] == 3:
        # If the sole neighbor is sp3, then in general it is not possible to calculate the
        # plane of the sp2-center.
        #
        # As a special case:
        # For OP1 and OP2 in DNA or RNA backbone, the plane of their sp2-centers is assumed
        # to be at maximum distance from each other, thus, perpendicular to the plane formed
        # by OP1-P-OP2.
        if (node_accep.endswith(":OP1")
                and acceptor_base_atoms[0].endswith(":P")
                and graphs['atom'].has_node(node_accep[:-3] + "OP2")):
            phosphor_pos = graphs['atom'].nodes[acceptor_base_atoms[0]]["coords"]
            acceptor_base_1_dist = phosphor_pos - acceptor_pos

            op2_pos = graphs['atom'].nodes[node_accep[:-3] + "OP2"]["coords"]
            op2_dist = op2_pos - phosphor_pos

            acceptor_base_2_dist = np.cross(acceptor_base_1_dist, op2_dist)

            acceptor_plane_normal = np.cross(acceptor_base_1_dist, acceptor_base_2_dist)
        elif (node_accep.endswith(":OP2")
                and acceptor_base_atoms[0].endswith(":P")
                and graphs['atom'].has_node(node_accep[:-3] + "OP1")):
            phosphor_pos = graphs['atom'].nodes[acceptor_base_atoms[0]]["coords"]
            acceptor_base_1_dist = phosphor_pos - acceptor_pos

            op1_pos = graphs['atom'].nodes[node_accep[:-3] + "OP1"]["coords"]
            op1_dist = op1_pos - phosphor_pos

            acceptor_base_2_dist = np.cross(acceptor_base_1_dist, op1_dist)

            acceptor_plane_normal = np.cross(acceptor_base_1_dist, acceptor_base_2_dist)
        else:
            raise ValueError(f"Unable to compute plane of sp2-center at {node_accep} since its sole neighbor is sp3.")

    if len(acceptor_base_atoms) < 2:
        acceptor_base_atoms += list(n for n in graphs['atom'][acceptor_base_atoms[0]] if n != node_accep)
    if len(acceptor_base_atoms) < 2:
        # If still missing third atom in plane, gamma cannot be calculated
        raise ValueError(f"Insufficient number of atoms to calculate plane of sp2-center at {node_accep}")

    if acceptor_plane_normal is None:
        acceptor_base_1_pos = graphs['atom'].nodes[acceptor_base_atoms[0]]["coords"]
        acceptor_base_1_dist = acceptor_base_1_pos - acceptor_pos

        acceptor_base_2_pos = graphs['atom'].nodes[acceptor_base_atoms[1]]["coords"]
        acceptor_base_2_dist = acceptor_base_2_pos - acceptor_pos

        acceptor_plane_normal = np.cross(acceptor_base_1_dist, acceptor_base_2_dist)

    gamma_angle = angle_between(donor_plane_normal, acceptor_plane_normal)
    if gamma_angle < 90.0:
        gamma_angle = 180 - gamma_angle  # use supplement of the angle

    return gamma_angle


def _mayo_energy_hydrogen_bond_angular_term_sp3_sp3(theta, phi):
    f = (np.cos(theta) ** 2) * np.exp(-(np.pi - theta) ** 6) * (np.cos(phi - np.radians(109.5)) ** 2)
    return f


def _mayo_energy_hydrogen_bond_angular_term_sp3_sp2(theta, phi):
    f = (np.cos(theta) ** 2) * np.exp(-(np.pi - theta) ** 6) * (np.cos(phi) ** 2)
    return f


def _mayo_energy_hydrogen_bond_angular_term_sp2_sp3(theta):
    f = (np.cos(theta) ** 4) * np.exp(-2 * (np.pi - theta) ** 6)
    return f


def _mayo_energy_hydrogen_bond_angular_term_sp2_sp2(theta, phi, gamma):
    f = (np.cos(theta) ** 2) * np.exp(-(np.pi - theta) ** 6) * (np.cos(max(phi, gamma)) ** 2)
    return f


def _mayo_energy_hydrogen_bond(d: float, angular_term: float):
    """ Angular dependent hydrogen bond energy function (kcal/mol) """
    r0_div_da = 2.8 / d  # R_0 / d
    rhs = r0_div_da ** 10  # power of 10
    lhs = rhs * r0_div_da * r0_div_da  # power of 12
    return 8.0 * (5 * lhs - 6 * rhs) * angular_term


def _mayo_energy_salt_bridge(d: float):
    """ Non-angular dependent salt bridge energy function (kcal/mol) """
    rs_div_da = 3.20 / (d + 0.375)  # R_s / (d + a)
    rhs = rs_div_da ** 10  # power of 10
    lhs = rhs * rs_div_da * rs_div_da  # power of 12
    return 10.0 * (5 * lhs - 6 * rhs)


################################
# Multiprocessing Function #####
################################

class HBondDictionary(TypedDict):
    donor: Required[str]
    acceptor: Required[str]
    hydrogen: NotRequired[str]
    energy: Required[float]
    bond_length: Required[float]
    bond_type: Required[str]
    extra: Required[List[float]]
    multi_base: NotRequired[bool]

def _handle_potential_hydrogen_bond(
        i: tuple,
        length: float,
        energy_threshold: float,
        graphs: GraphDictionary,
        pdb_df: pd.DataFrame, hydrogen_mapping: pd.DataFrame,
        hydrogen_df: pd.DataFrame,
        cation_atoms: pd.DataFrame, anion_atoms: pd.DataFrame,
        hybridizations: pd.DataFrame
) -> Optional[List[HBondDictionary]]:
    node_1 = pdb_df["node_id"][i[0]]
    node_2 = pdb_df["node_id"][i[1]]
    chain_1 = pdb_df["chain_id"][i[0]]
    chain_2 = pdb_df["chain_id"][i[1]]
    resi_1 = pdb_df["residue_number"][i[0]]
    resi_2 = pdb_df["residue_number"][i[1]]

    if chain_1 == chain_2 and resi_1 == resi_2:
        return None

    bond_type = InteractionType.H_BOND.value
    donor_id, accep_id = i[0], i[1]
    node_donor, chain_donor, resi_donor = node_1, chain_1, resi_1
    node_accep, chain_accep, resi_accep = node_2, chain_2, resi_2

    if cation_atoms['node_id'].eq(node_1).any() and anion_atoms['node_id'].eq(node_2).any():
        bond_type = InteractionType.SALT_BRIDGE.value
        r_threshold = 3.6
    elif (':S' in node_1 or ':S' in node_2) and length <= 4.2:
        # distance constraints for sulfur case hbond
        r_threshold = 3.0  # max hydrogen acceptor distance
    elif length <= 3.6:
        # distance constraints for standard case hbond
        r_threshold = 2.6  # max hydrogen acceptor distance
    else:
        return None

    donor_pos = pdb_df.loc[pdb_df['node_id'] == node_donor, ["x_coord", "y_coord", "z_coord"]].to_numpy()[0]
    acceptor_pos = pdb_df.loc[pdb_df['node_id'] == node_accep, ["x_coord", "y_coord", "z_coord"]].to_numpy()[0]

    # storage anchor for gamma angle between two sp2-planes
    gamma_angle = None

    if not hydrogen_mapping['node_id'].eq(node_donor).any():
        _tqdm_logger.info(f"No valid hydrogen for potential donor {node_donor}.")

    hbond_records = []

    # iterate over all H atoms bonded to donor
    for j, row in hydrogen_mapping.loc[hydrogen_mapping['node_id'] == node_donor, :].iterrows():
        # h_pos = np.array((row['x_coord'], row['y_coord'], row['z_coord']))
        selected_h_atom = row['h_id']
        h_pos = hydrogen_df.loc[hydrogen_df['node_id'] == selected_h_atom, ["x_coord", "y_coord", "z_coord"]].to_numpy()[0]

        # calculate hydrogen - acceptor vector
        ha_dist = h_pos - acceptor_pos
        ha_len = np.linalg.norm(ha_dist)

        # check hydrogen - acceptor distance
        if ha_len > r_threshold:
            continue

        # get donor - hydrogen vector
        dh_dist = h_pos - donor_pos
        # donor - hydrogen - acceptor angle
        angle_theta = angle_between(dh_dist, ha_dist)
        if angle_theta < 80:
            continue

        if bond_type == InteractionType.SALT_BRIDGE.value:
            bond_energy = _mayo_energy_salt_bridge(length)

            # sanity check
            if bond_energy > energy_threshold:
                # don't need to consider other H-atoms as energy function is non-angular dependent
                return None

            # don't need to consider other H-atoms as energy function is non-angular dependent
            return [HBondDictionary(
                donor=node_donor, acceptor=node_accep, energy=bond_energy, bond_length=length,
                bond_type=bond_type, extra=[angle_theta]
            )]

        donor_sp, accep_sp = pdb_df["hybridization"][donor_id], pdb_df["hybridization"][accep_id]

        #############################
        # 4 Hybridization cases #####
        #############################

        # bool for h-bond variant (multi-base-acceptor)
        multi_base_acc = len(graphs['atom'][node_accep]) > 1

        if donor_sp == 2 and accep_sp == 3:
            angular_term = _mayo_energy_hydrogen_bond_angular_term_sp2_sp3(np.radians(angle_theta))
            bond_energy = _mayo_energy_hydrogen_bond(length, angular_term)

            if bond_energy > energy_threshold:
                continue

            # store info for later filtering
            hbond_records.append(HBondDictionary(
                donor=node_donor, acceptor=node_accep, hydrogen=selected_h_atom,
                energy=bond_energy, bond_length=ha_len, bond_type=bond_type, extra=[angle_theta],
                multi_base=multi_base_acc
            ))
            continue

        if donor_sp == 2 and accep_sp == 2 and gamma_angle is None:
            # pre-calculate gamma angle
            gamma_angle = _calculate_gamma_angle(
                graphs=graphs,
                node_donor=node_donor, node_accep=node_accep,
                donor_pos=donor_pos, acceptor_pos=acceptor_pos,
                dh_dist=dh_dist,
                hybridization_df=hybridizations
            )

        # angle phi required
        hbond_record = None
        for base_atom in graphs['atom'][node_accep]:
            # select base atom for the highest
            base_atom_pos = graphs['atom'].nodes[base_atom]["coords"]
            base_atom_dist = base_atom_pos - acceptor_pos

            phi_angle = angle_between(ha_dist, base_atom_dist)

            if donor_sp == 3 and accep_sp == 3:
                angular_term = _mayo_energy_hydrogen_bond_angular_term_sp3_sp3(
                    np.radians(angle_theta), np.radians(phi_angle)
                )
            elif donor_sp == 3 and accep_sp == 2:
                angular_term = _mayo_energy_hydrogen_bond_angular_term_sp3_sp2(
                    np.radians(angle_theta), np.radians(phi_angle)
                )
            elif donor_sp == 2 and accep_sp == 2 and gamma_angle is not None:
                angular_term = _mayo_energy_hydrogen_bond_angular_term_sp2_sp2(
                    np.radians(angle_theta), np.radians(phi_angle), np.radians(gamma_angle)
                )
            else:
                raise ValueError('Should not happen')
            bond_energy = _mayo_energy_hydrogen_bond(length, angular_term)

            if bond_energy > energy_threshold:
                continue

            if hbond_record is None or hbond_record['energy'] > bond_energy:
                extra_list = [angle_theta, phi_angle, gamma_angle] if gamma_angle is not None else [angle_theta, phi_angle]
                hbond_record = HBondDictionary(
                    donor=node_donor, acceptor=node_accep, hydrogen=selected_h_atom,
                    energy=bond_energy, bond_length=ha_len, bond_type=bond_type, extra=extra_list,
                    multi_base=multi_base_acc
                )

        if hbond_record is not None:
            hbond_records.append(hbond_record)

    return None if len(hbond_records) == 0 else hbond_records


#####################
# Main Function #####
#####################

def calculate_hydrogen_and_salt_bridge_bonds(
        graphs: GraphDictionary,
        heavy_atom_df: pd.DataFrame,
        hydrogen_mapping: pd.DataFrame,
        hydrogen_df: pd.DataFrame,
        minimum_distance: float = 2.6,
        energy_threshold: float = -0.1,
        verbose: bool = False
) -> "Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]":
    pdb_df = heavy_atom_df.copy()

    # label dataframe
    pdb_df = _label_h_bond_donors_acceptors(pdb_df=pdb_df, h_frame=hydrogen_mapping)
    # filter and label hybridization
    pdb_df: pd.DataFrame = _label_hybridization(pdb_df=pdb_df)
    hybridizations: pd.DataFrame = pdb_df.loc[:, ['node_id', 'hybridization']]
    pdb_df = pdb_df.loc[(pdb_df['hbond_donor'] > 0) | (pdb_df['hbond_acceptor'] > 0)].reset_index(drop=True)

    # get valid donor and acceptor atoms
    valid_donors, valid_acceptors = _get_valid_donor_acceptor_atoms(
        pdb_df=pdb_df, h_frame=hydrogen_mapping
    )

    if len(valid_donors) < 1 or len(valid_acceptors) < 1:
        return None, None

    # get valid cation and anion atoms for salt bridges
    cation_atoms, anion_atoms = _get_valid_cation_anion_atoms(
        h_frame=hydrogen_mapping, valid_donors=valid_donors, valid_acceptors=valid_acceptors
    )

    num_donors = len(valid_donors)
    pdb_df = pd.concat([valid_donors, valid_acceptors], ignore_index=True)
    # pdb_df may contain duplicates, DON'T remove

    # calculate distances with cut-off
    distmat = compute_distmat(pdb_df=pdb_df)
    distmat = distmat[distmat >= minimum_distance]
    t_distmat = distmat[distmat <= 4.6]  # max distance for salt bridges

    salt_bridge_records = []
    hbond_records = []

    inds = list(i for i in zip(*np.where(~np.isnan(t_distmat))) if i[0] < num_donors <= i[1])
    with tqdm(iterable=enumerate(inds), total=len(inds), desc="Test for H-Bonds", disable=not verbose) as progress_bar:
        for i, ind in progress_bar:
            result = _handle_potential_hydrogen_bond(
                ind,
                t_distmat[ind[0]][ind[1]],
                energy_threshold,
                graphs,
                pdb_df, hydrogen_mapping,
                hydrogen_df,
                cation_atoms, anion_atoms,
                hybridizations
            )

            if result is None:
                continue

            if result[0]['bond_type'] == InteractionType.H_BOND.value:
                hbond_records += result
            else:
                salt_bridge_records += result
            if verbose:
                progress_bar.set_postfix({'HB': len(hbond_records), 'SB': len(salt_bridge_records)})

    # combine and filter records
    salt_bridge_frame = pd.DataFrame.from_records(salt_bridge_records)
    hbond_frame = pd.DataFrame.from_records(hbond_records)

    if len(salt_bridge_frame) == 0 and len(hbond_frame) == 0:
        return None, None

    if len(salt_bridge_frame) > 0:
        salt_bridge_frame.rename(columns={'donor': 'node_1', 'acceptor': 'node_2'}, inplace=True)

    if len(hbond_frame) > 0:
        # sort with strongest bonds first
        hbond_frame = hbond_frame.sort_values(by=['energy'], ascending=True).reset_index(drop=True)

        # export covalent bonds for required hydrogen atoms
        hydrogen_list = []
        covalent_records = []

        for i, row in hbond_frame.iterrows():
            node_1, node_2 = row['donor'], row['hydrogen']
            if node_2 not in hydrogen_list:
                covalent_records.append({
                    'node_1': node_1, 'node_2': node_2, 'bond_type': InteractionType.COVALENT.value,
                    'bond_length': hydrogen_mapping.loc[hydrogen_mapping.h_id == node_2, 'length'].to_numpy()[0]
                })
                hydrogen_list.append(node_2)

        covalent_frame = pd.DataFrame.from_records(covalent_records)

        # format hydrogen bonds between hydrogen atom and acceptor atom
        hbond_frame: pd.DataFrame = hbond_frame.loc[
            :, [col for col in hbond_frame.columns if col not in ['donor', 'multi_base']]].reset_index(drop=True)
        hbond_frame.rename(columns={'hydrogen': 'node_1', 'acceptor': 'node_2'}, inplace=True)
    else:
        covalent_frame = None

    bond_frame = pd.concat([hbond_frame, salt_bridge_frame], ignore_index=True)

    return bond_frame, covalent_frame
