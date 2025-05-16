from typing import List, Union
import pandas as pd

__all__ = [
    "convert_to_pdb_lines",
    "load_as_atom_df", "split_into_heavy_atoms_and_hydrogens", "load_as_others_df", "load_as_hydrogen_mapping"
]


def convert_to_pdb_lines(coords: List[List[Union[str, int, float]]], others: List[List[Union[str, int]]]):
    """
    coords in the order of:

    ``"record_name"``, ``"atom_number"``,
    ``"chain_id"``, ``"residue_number"``, ``"insertion"``, ``"residue_name"``, ``"atom_name"``,
    ``"node_id"``,
    ``"element_symbol"``, ``"charge"``,
    ``"x_coord"``, ``"y_coord"``, ``"z_coord"``,
    ``"line_idx"``

    others in the order of:

    ``"record_name"``, ``"entry"``, ``"line_idx"``

    Parameters
    ----------
    coords
    others

    Returns
    -------

    """
    return [
        line + "\n" for _, line in sorted(
            [
                (
                    entry[-1],
                    f"{entry[0]: <6}{entry[1]:5d} {(' 'if len(entry[6]) < 4 else '') + entry[6]: <4} {entry[5]: >3} {entry[2]}{entry[3]:4d}    {entry[10]:8.3f}{entry[11]:8.3f}{entry[12]:8.3f}  0.00  0.00{' ' * 10}{entry[8]: >2}{entry[9]: >2}"
                )
                for entry in coords
            ] + [
                (entry[-1], f"{entry[0]: <6}{entry[1]: <74}")
                for entry in others
            ],
            key=lambda x: x[0]
        )
    ]

def load_as_atom_df(data: List[List[Union[str, int, float]]]):
    """
    in the order of:

    ``"record_name"``, ``"atom_number"``,
    ``"chain_id"``, ``"residue_number"``, ``"insertion"``, ``"residue_name"``, ``"atom_name"``,
    ``"node_id"``,
    ``"element_symbol"``, ``"charge"``,
    ``"x_coord"``, ``"y_coord"``, ``"z_coord"``,
    ``"line_idx"``

    Parameters
    ----------
    data

    Returns
    -------

    """
    return pd.DataFrame(data, columns=[
        "record_name", "atom_number",
        "chain_id", "residue_number", "insertion", "residue_name", "atom_name",
        "node_id",
        "element_symbol", "charge",
        "x_coord", "y_coord", "z_coord",
        "line_idx"
    ])

def split_into_heavy_atoms_and_hydrogens(atom_df: pd.DataFrame):
    return (atom_df.copy().loc[atom_df.element_symbol != "H", :].reset_index(drop=True),
            atom_df.copy().loc[atom_df.element_symbol == "H", :].reset_index(drop=True))

def load_as_others_df(data: List[List[Union[str, int]]]):
    """

    ``"record_name"``, ``"entry"``, ``"line_idx"``

    Parameters
    ----------
    data

    Returns
    -------

    """
    return pd.DataFrame(data, columns=["record_name", "entry", "line_idx"])

def load_as_hydrogen_mapping(data: List[List[Union[str, float]]]):
    """

    ``"h_id"``, ``"node_id"``, ``"length"``

    Parameters
    ----------
    data

    Returns
    -------

    """
    return pd.DataFrame(data, columns=["h_id", "node_id", "length"])
