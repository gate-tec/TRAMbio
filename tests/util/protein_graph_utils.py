from typing import List, Literal, Union, Optional, Tuple
import pandas as pd
import networkx as nx
import numpy as np

from TRAMbio.util.structure_library.graph_struct import GraphDictionary, GraphKey, ProteinGraph, \
    initialize_graphs_from_dataframe

__all__ = [
    "construct_protein_graph_base", "modify_coords", "format_coords", "apply_rotations",
    "AXIS_OPTIONS",
    "convert_to_pdb_lines",
    "load_as_atom_df", "split_into_heavy_atoms_and_hydrogens", "load_as_others_df", "load_as_hydrogen_mapping"
]


def construct_protein_graph_base(coords, others, hm, cov, pe, pm, pb):
    # Load dataframes
    atom_df = load_as_atom_df(coords)
    heavy_atom_df, hydrogen_df = split_into_heavy_atoms_and_hydrogens(atom_df=atom_df)
    others_df = load_as_others_df(others)
    hydrogen_mapping = load_as_hydrogen_mapping(hm)

    # Construct base graph
    graphs: GraphDictionary = initialize_graphs_from_dataframe(atom_df=atom_df, heavy_atom_df=heavy_atom_df)

    # Place covalent edges
    graphs["full"].add_edges_from(cov)
    graphs["atom"].add_edges_from([entry for entry in cov if entry[0][10] != "H" and entry[1][10] != "H"])

    # Construct pebbled graph
    graphs["pebble"].add_edges_from(pe)
    nx.set_node_attributes(graphs["pebble"], pm, "pebbles")
    graphs["pebble"].graph[GraphKey.STANDARD_EDGES.value] += pb  # noqa: PyChram doesn't recognize the graph dictionary

    return ProteinGraph(
        graphs,
        atom_df=atom_df,
        others_df=others_df,
        heavy_atom_df=heavy_atom_df,
        hydrogen_df=hydrogen_df,
        hydrogen_mapping=hydrogen_mapping
    )


def format_coords(
        coords, hydrogen_mapping, covalent_edges, pebble_edges, pebble_map, pebble_bars, exclude
):
    if exclude is not None:
        coords = [entry for entry in coords if entry[6] not in exclude]
        hydrogen_mapping = [entry for entry in hydrogen_mapping if entry[0][10:] not in exclude and entry[1][10:] not in exclude]
        covalent_edges = [entry for entry in covalent_edges if entry[0][10:] not in exclude and entry[1][10:] not in exclude]
        pebble_edges = [entry for entry in pebble_edges if entry[0][10:] not in exclude and entry[1][10:] not in exclude]
        pebble_map = {k: v for k, v in pebble_map.items() if k[10:] not in exclude}
        pebble_bars = [entry for entry in pebble_bars if entry[0][10:] not in exclude and entry[1][10:] not in exclude]

    coords = [
        entry[:10] + [float(f"{entry[10]:8.3f}"), float(f"{entry[11]:8.3f}"), float(f"{entry[12]:8.3f}")] + entry[13:]
        for entry in coords
    ]

    return coords, hydrogen_mapping, covalent_edges, pebble_edges, pebble_map, pebble_bars


def modify_coords(coords, x_offset, y_offset, z_offset):
    if x_offset == y_offset == z_offset == 0.0:
        return coords

    return [
        entry[:10] +
        [
            a + b
            for a, b in zip([x_offset, y_offset, z_offset], entry[10:13])] +
        entry[13:]
        for entry in coords
    ]


def _apply_z_rotation(coords, z_rotation):
    if z_rotation == 0.0:
        return coords

    theta = np.radians(z_rotation)
    new_coords = []
    for entry in coords:
        x = entry[10]
        y = entry[11]
        z = entry[12]

        new_x = round(x * np.cos(theta) - y * np.sin(theta), ndigits=4)
        new_y = round(x * np.sin(theta) + y * np.cos(theta), ndigits=4)
        new_z = z

        new_coords.append(entry[:10] + [new_x, new_y, new_z] + entry[13:])

    return new_coords


def _apply_y_rotation(coords, y_rotation):
    if y_rotation == 0.0:
        return coords

    theta = np.radians(y_rotation)
    new_coords = []
    for entry in coords:
        x = entry[10]
        y = entry[11]
        z = entry[12]

        new_x = round(x * np.cos(theta) + z * np.sin(theta), ndigits=4)
        new_y = y
        new_z = round(-x * np.sin(theta) + z * np.cos(theta), ndigits=4)

        new_coords.append(entry[:10] + [new_x, new_y, new_z] + entry[13:])

    return new_coords


def _apply_x_rotation(coords, x_rotation):
    if x_rotation == 0.0:
        return coords

    theta = np.radians(x_rotation)
    new_coords = []
    for entry in coords:
        x = entry[10]
        y = entry[11]
        z = entry[12]

        new_x = x
        new_y = round(y * np.cos(theta) - z * np.sin(theta), ndigits=4)
        new_z = round(y * np.sin(theta) + z * np.cos(theta), ndigits=4)

        new_coords.append(entry[:10] + [new_x, new_y, new_z] + entry[13:])

    return new_coords

AXIS_OPTIONS = Literal['X', 'Y', 'Z']
_AXIS_MAPPING = {
    'X': _apply_x_rotation,
    'Y': _apply_y_rotation,
    'Z': _apply_z_rotation
}


def apply_rotations(coords, rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None):
    if rotations is not None:
        for angle, axis in rotations:
            coords = _AXIS_MAPPING[axis](coords, angle)
    return coords


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
