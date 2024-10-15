import pandas as pd

from scipy.spatial.distance import pdist, squareform

from TRAMbio.util.graphein.constants import COVALENT_RADII, RESIDUE_ATOM_BOND_STATE, DEFAULT_BOND_STATE


#################
# Functions #####
#################


def compute_distmat(pdb_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute pairwise Euclidean distances between every atom.

    Design choice: passed in a ``pd.DataFrame`` to enable easier testing on
    dummy data.

    :param pdb_df: Dataframe containing protein structure. Must contain columns
        ``["x_coord", "y_coord", "z_coord"]``.
    :type pdb_df: pd.DataFrames
    :raises: ValueError if ``pdb_df`` does not contain the required columns.
    :return: pd.Dataframe of Euclidean distance matrix.
    :rtype: pd.DataFrame
    """
    if (
        not pd.Series(["x_coord", "y_coord", "z_coord"])
        .isin(pdb_df.columns)
        .all()
    ):
        raise ValueError(
            "Dataframe must contain columns ['x_coord', 'y_coord', 'z_coord']"
        )
    eucl_dists = pdist(
        pdb_df[["x_coord", "y_coord", "z_coord"]], metric="euclidean"
    )
    eucl_dists = pd.DataFrame(squareform(eucl_dists))
    eucl_dists.index = pdb_df.index
    eucl_dists.columns = pdb_df.index

    return eucl_dists


def assign_bond_states_to_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Takes a ``PandasPDB`` atom DataFrame and assigns bond states to each atom
    based on:

        *Atomic Structures of all the Twenty Essential Amino Acids and a*
        *Tripeptide, with Bond Lengths as Sums of Atomic Covalent Radii*
        Heyrovska, 2008

    First, maps atoms to their standard bond states
    (:const:`~tram.io_interface.graphein_copy.DEFAULT_BOND_STATE`). Second, maps
    non-standard bonds states
    (:const:`~tram.io_interface.graphein_copy.RESIDUE_ATOM_BOND_STATE`). Fills
    ``NaNs`` with standard bond states.

    :param df: Pandas PDB DataFrame.
    :type df: pd.DataFrame
    :return: DataFrame with added ``atom_bond_state`` column.
    :rtype: pd.DataFrame
    """

    # Map atoms to their standard bond states
    naive_bond_states = pd.Series(df["atom_name"].map(DEFAULT_BOND_STATE))

    # Create series of bond states for the non-standard states
    ss = (
        pd.DataFrame(RESIDUE_ATOM_BOND_STATE)
        .unstack()
        .rename_axis(("residue_name", "atom_name"))
        .rename("atom_bond_state")
    )

    # Map non-standard states to the DataFrame based on the residue and atom
    # name
    df = df.join(ss, on=["residue_name", "atom_name"])

    # Fill the NaNs with the standard states
    df = df.fillna(value={"atom_bond_state": naive_bond_states})

    return df


def assign_covalent_radii_to_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Assigns covalent radius
    (:const:`tram.io_interface.graphein_copy.COVALENT_RADII`) to each atom based
    on its bond state. Adds a ``covalent_radius`` column. Using values from:

        *Atomic Structures of all the Twenty Essential Amino Acids and a*
        *Tripeptide, with Bond Lengths as Sums of Atomic Covalent Radii*
        Heyrovska, 2008

    :param df: Pandas PDB DataFrame with a ``bond_states_column``.
    :type df: pd.DataFrame
    :return: Pandas PDB DataFrame with added ``covalent_radius`` column.
    :rtype: pd.DataFrame
    """
    # Assign covalent radius to each atom
    df["covalent_radius"] = df["atom_bond_state"].map(COVALENT_RADII)

    return df
