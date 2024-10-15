import math
from typing import Optional, List, Literal, Tuple

import networkx as nx
import numpy as np

from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import GraphDictionary, GraphKey, ProteinGraph, \
    initialize_graphs_from_dataframe
from tests.util.protein_graph_utils import *


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


def _format(
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
        entry[:10] + [float("{:8.3f}".format(entry[10])), float("{:8.3f}".format(entry[11])), float("{:8.3f}".format(entry[12]))] + entry[13:]
        for entry in coords
    ]

    return coords, hydrogen_mapping, covalent_edges, pebble_edges, pebble_map, pebble_bars


def _modify_coords(coords, x_offset, y_offset, z_offset):
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

_AXIS_OPTIONS = Literal['X', 'Y', 'Z']
_AXIS_MAPPING = {
    'X': _apply_x_rotation,
    'Y': _apply_y_rotation,
    'Z': _apply_z_rotation
}

def _apply_rotations(coords, rotations: Optional[List[Tuple[float, _AXIS_OPTIONS]]] = None):
    if rotations is not None:
        for angle, axis in rotations:
            coords = _AXIS_MAPPING[axis](coords, angle)
    return coords


def graph_gly_donor(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        rotations: Optional[List[Tuple[float, _AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    coords = _modify_coords(_apply_rotations([
        ["ATOM", start_idx, "A", res, "", "GLY", "CA", f"A{res:04d}-GLY:CA", "C", "",
         -1.07 -math.cos(math.radians(180 - 109.5)) * 1.47, +math.sin(math.radians(180 - 109.5)) * 1.47, 0, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", "GLY", "N", f"A{res:04d}-GLY:N", "N", "",
         -1.07, 0, 0, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", "GLY", "H", f"A{res:04d}-GLY:H", "H", "",
         0, 0, 0, start_idx + 2]
    ], rotations), x_offset, y_offset, z_offset)
    hydrogen_mapping = [
        (f"A{res:04d}-GLY:H", f"A{res:04d}-GLY:N", 1.07)
    ]
    covalent_edges = [
        (f"A{res:04d}-GLY:N", f"A{res:04d}-GLY:H",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-GLY:CA", f"A{res:04d}-GLY:N",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.47, "base": True})
    ]
    pebble_edges = [
        (f"A{res:04d}-GLY:CA", f"A{res:04d}-GLY:N", {"weight": 0}),
        (f"A{res:04d}-GLY:N", f"A{res:04d}-GLY:CA", {"weight": 5})
    ]
    pebble_map = {
        f"A{res:04d}-GLY:N": 1,
        f"A{res:04d}-GLY:CA": 6
    }
    pebble_bars = []

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )


def graph_gly_acceptor(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        rotations: Optional[List[Tuple[float, _AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    coords = _modify_coords(_apply_rotations([
        ["ATOM", start_idx, "A", res, "", "GLY", "CA", f"A{res:04d}-GLY:CA", "C", "",
         1.27 + (0.5 * 1.54), +math.sin(math.radians(60)) * 1.54, 0, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", "GLY", "C", f"A{res:04d}-GLY:C", "C", "",
         1.27, 0, 0, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", "GLY", "O", f"A{res:04d}-GLY:O", "O", "",
         0, 0, 0, start_idx + 2]
    ], rotations), x_offset, y_offset, z_offset)
    hydrogen_mapping = []
    covalent_edges = [
        (f"A{res:04d}-GLY:CA", f"A{res:04d}-GLY:C",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.54, "base": True}),
        (f"A{res:04d}-GLY:C", f"A{res:04d}-GLY:O",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.27, "base": True})
    ]
    pebble_edges = [
        (f"A{res:04d}-GLY:C", f"A{res:04d}-GLY:CA", {"weight": 0}),
        (f"A{res:04d}-GLY:CA", f"A{res:04d}-GLY:C", {"weight": 5}),
        (f"A{res:04d}-GLY:C", f"A{res:04d}-GLY:O", {"weight": 0}),
        (f"A{res:04d}-GLY:O", f"A{res:04d}-GLY:C", {"weight": 5})
    ]
    pebble_map = {
        f"A{res:04d}-GLY:CA": 1,
        f"A{res:04d}-GLY:O": 1,
        f"A{res:04d}-GLY:C": 6
    }
    pebble_bars = [
        (f"A{res:04d}-GLY:O", f"A{res:04d}-GLY:C", 1)
    ]

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )


def graph_asn_donor(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        exclude: Optional[List[str]] = None
):
    coords = [
        ["ATOM", start_idx, "A", res, "", "ASN", "CG", f"A{res:04d}-ASN:CG", "C", "",
         -(1.07 * 0.5) - 1.38 + x_offset, y_offset, z_offset, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", "ASN", "OD1", f"A{res:04d}-ASN:OD1", "O", "",
         -(1.07 * 0.5) - 1.38 + x_offset, 1.27 + y_offset, z_offset, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", "ASN", "ND2", f"A{res:04d}-ASN:ND2", "N", "",
         -(1.07 * 0.5) + x_offset, y_offset, z_offset, start_idx + 2],
        ["ATOM", start_idx + 3, "A", res, "", "ASN", "HD21", f"A{res:04d}-ASN:HD21", "H", "",
         x_offset, +math.sin(math.radians(60)) * 1.07 + y_offset, z_offset, start_idx + 3],
        ["ATOM", start_idx + 4, "A", res, "", "ASN", "HD22", f"A{res:04d}-ASN:HD22", "H", "",
         x_offset, -math.sin(math.radians(60)) * 1.07 + y_offset, z_offset, start_idx + 4]
    ]
    hydrogen_mapping = [
        (f"A{res:04d}-ASN:HD21", f"A{res:04d}-ASN:ND2", 1.07),
        (f"A{res:04d}-ASN:HD22", f"A{res:04d}-ASN:ND2", 1.07)
    ]
    covalent_edges = [
        (f"A{res:04d}-ASN:ND2", f"A{res:04d}-ASN:HD21",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-ASN:ND2", f"A{res:04d}-ASN:HD22",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-ASN:CG", f"A{res:04d}-ASN:ND2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.38, "base": True}),
        (f"A{res:04d}-ASN:CG", f"A{res:04d}-ASN:OD1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.27, "base": True})
    ]
    pebble_edges = [
        (f"A{res:04d}-ASN:ND2", f"A{res:04d}-ASN:HD21", {"weight": 0}),
        (f"A{res:04d}-ASN:HD21", f"A{res:04d}-ASN:ND2", {"weight": 5}),
        (f"A{res:04d}-ASN:ND2", f"A{res:04d}-ASN:HD22", {"weight": 0}),
        (f"A{res:04d}-ASN:HD22", f"A{res:04d}-ASN:ND2", {"weight": 5}),
        (f"A{res:04d}-ASN:CG", f"A{res:04d}-ASN:ND2", {"weight": 0}),
        (f"A{res:04d}-ASN:ND2", f"A{res:04d}-ASN:CG", {"weight": 5}),
        (f"A{res:04d}-ASN:CG", f"A{res:04d}-ASN:OD1", {"weight": 0}),
        (f"A{res:04d}-ASN:OD1", f"A{res:04d}-ASN:CG", {"weight": 5})
    ]
    pebble_map = {
        f"A{res:04d}-ASN:HD21": 1,
        f"A{res:04d}-ASN:HD22": 1,
        f"A{res:04d}-ASN:ND2": 1,
        f"A{res:04d}-ASN:OD1": 1,
        f"A{res:04d}-ASN:CG": 6
    }
    pebble_bars = [
        (f"A{res:04d}-ASN:OD1", f"A{res:04d}-ASN:CG", 1)
    ]

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )


def graph_arg_donor(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        exclude: Optional[List[str]] = None
):
    coords = [
        ["ATOM", start_idx, "A", res, "", "ARG", "NH1", f"A{res:04d}-ARG:NH1", "N", "1+",
         -1.07 + x_offset, +math.sin(math.radians(60)) * 1.38 + y_offset, z_offset, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", "ARG", "HH11", f"A{res:04d}-ARG:HH11", "H", "",
         -1.07 - (0.5 * 1.07) + x_offset, +math.sin(math.radians(60)) * (1.38 + 1.07) + y_offset, z_offset, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", "ARG", "HH12", f"A{res:04d}-ARG:HH12", "H", "",
         x_offset, +math.sin(math.radians(60)) * 1.38 + y_offset, z_offset, start_idx + 2],
        ["ATOM", start_idx + 3, "A", res, "", "ARG", "NH2", f"A{res:04d}-ARG:NH2", "N", "",
         -1.07 + x_offset, -math.sin(math.radians(60)) * 1.38 + y_offset, z_offset, start_idx + 3],
        ["ATOM", start_idx + 4, "A", res, "", "ARG", "HH21", f"A{res:04d}-ARG:HH21", "H", "",
         x_offset, -math.sin(math.radians(60)) * 1.38 + y_offset, z_offset, start_idx + 4],
        ["ATOM", start_idx + 5, "A", res, "", "ARG", "HH22", f"A{res:04d}-ARG:HH22", "H", "",
         -1.07 - (0.5 * 1.07) + x_offset, -math.sin(math.radians(60)) * (1.38 + 1.07) + y_offset, z_offset, start_idx + 5],
        ["ATOM", start_idx + 6, "A", res, "", "ARG", "CZ", f"A{res:04d}-ARG:CZ", "C", "",
         -1.07 - (0.5 * 1.38) + x_offset, y_offset, z_offset, start_idx + 6]
    ]
    hydrogen_mapping = [
        [f"A{res:04d}-ARG:HH11", f"A{res:04d}-ARG:NH1", 1.07],
        [f"A{res:04d}-ARG:HH12", f"A{res:04d}-ARG:NH1", 1.07],
        [f"A{res:04d}-ARG:HH21", f"A{res:04d}-ARG:NH2", 1.07],
        [f"A{res:04d}-ARG:HH22", f"A{res:04d}-ARG:NH2", 1.07]
    ]
    covalent_edges = [
        (f"A{res:04d}-ARG:NH1", f"A{res:04d}-ARG:HH11",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-ARG:NH1", f"A{res:04d}-ARG:HH12",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-ARG:NH2", f"A{res:04d}-ARG:HH21",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-ARG:NH2", f"A{res:04d}-ARG:HH22",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-ARG:CZ", f"A{res:04d}-ARG:NH1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.38, "base": True}),
        (f"A{res:04d}-ARG:CZ", f"A{res:04d}-ARG:NH2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.38, "base": True})
    ]
    pebble_edges = [
        (f"A{res:04d}-ARG:CZ", f"A{res:04d}-ARG:NH1", {"weight": 0}),
        (f"A{res:04d}-ARG:NH1", f"A{res:04d}-ARG:CZ", {"weight": 5}),
        (f"A{res:04d}-ARG:CZ", f"A{res:04d}-ARG:NH2", {"weight": 0}),
        (f"A{res:04d}-ARG:NH2", f"A{res:04d}-ARG:CZ", {"weight": 5})
    ]
    pebble_map = {
        f"A{res:04d}-ARG:NH1": 1,
        f"A{res:04d}-ARG:NH2": 1,
        f"A{res:04d}-ARG:CZ": 6,
    }
    pebble_bars = [
        (f"A{res:04d}-ARG:NH2", f"A{res:04d}-ARG:NH1", 2)
    ]

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )


def graph_glu_acceptor(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        rotations: Optional[List[Tuple[float, _AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    coords = _modify_coords(_apply_rotations([
        ["ATOM", start_idx, "A", res, "", "GLU", "OE1", f"A{res:04d}-GLU:OE1", "O", "",
         0, +math.sin(math.radians(60)) * 1.355, 0, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", "GLU", "OE2", f"A{res:04d}-GLU:OE2", "O", "1-",
         0, -math.sin(math.radians(60)) * 1.355, 0, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", "GLU", "CD", f"A{res:04d}-GLU:CD", "C", "",
         (0.5 * 1.355), 0, 0, start_idx + 2]
    ], rotations), x_offset, y_offset, z_offset)
    hydrogen_mapping = []
    covalent_edges = [
        (f"A{res:04d}-GLU:CD", f"A{res:04d}-GLU:OE1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.35, "base": True}),
        (f"A{res:04d}-GLU:CD", f"A{res:04d}-GLU:OE2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.35, "base": True})
    ]
    pebble_edges = [
        (f"A{res:04d}-GLU:CD", f"A{res:04d}-GLU:OE1", {"weight": 0}),
        (f"A{res:04d}-GLU:OE1", f"A{res:04d}-GLU:CD", {"weight": 5}),
        (f"A{res:04d}-GLU:CD", f"A{res:04d}-GLU:OE2", {"weight": 0}),
        (f"A{res:04d}-GLU:OE2", f"A{res:04d}-GLU:CD", {"weight": 5})
    ]
    pebble_map = {
        f"A{res:04d}-GLU:OE1": 1,
        f"A{res:04d}-GLU:OE2": 1,
        f"A{res:04d}-GLU:CD": 6
    }
    pebble_bars = [
        (f"A{res:04d}-GLU:OE2", f"A{res:04d}-GLU:OE1", 2)
    ]

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )


def graph_ser_acceptor(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        exclude: Optional[List[str]] = None
):
    coords = [
        ["ATOM", start_idx, "A", res, "", "SER", "HG", f"A{res:04d}-SER:HG", "H", "",
         x_offset, 1.04 + y_offset, z_offset, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", "SER", "OG", f"A{res:04d}-SER:OG", "O", "",
         x_offset, y_offset, z_offset, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", "SER", "CB", f"A{res:04d}-SER:CB", "C", "",
         +math.cos(math.radians(104.5 - 90)) * 1.44 + x_offset, -math.sin(math.radians(104.5 - 90)) * 1.44 + y_offset, z_offset, start_idx + 2]
    ]
    hydrogen_mapping = [
        [f"A{res:04d}-SER:HG", f"A{res:04d}-SER:OG", 1.04]
    ]
    covalent_edges = [
        (f"A{res:04d}-SER:OG", f"A{res:04d}-SER:HG",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.35, "base": True}),
        (f"A{res:04d}-SER:CB", f"A{res:04d}-SER:OG",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.35, "base": True})
    ]
    pebble_edges = [
        (f"A{res:04d}-SER:CB", f"A{res:04d}-SER:OG", {"weight": 0}),
        (f"A{res:04d}-SER:OG", f"A{res:04d}-SER:CB", {"weight": 5})
    ]
    pebble_map = {
        f"A{res:04d}-SER:OG": 1,
        f"A{res:04d}-SER:CB": 6
    }
    pebble_bars = []

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )


def graph_lys_donor(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        exclude: Optional[List[str]] = None
):
    unit_dist = math.sin(math.radians(180 - 109.5)) * 1.07
    coords = [
        ["ATOM", start_idx, "A", res, "", "LYS", "NZ", f"A{res:04d}-LYS:NZ", "N", "1+",
         -math.cos(math.radians(180 - 109.5)) * 1.07 + x_offset, y_offset, z_offset, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", "LYS", "HZ1", f"A{res:04d}-LYS:HZ1", "H", "",
         x_offset, +unit_dist + y_offset, z_offset, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", "LYS", "HZ2", f"A{res:04d}-LYS:HZ2", "H", "",
         x_offset, -math.sin(math.radians(30)) * unit_dist + y_offset, +math.cos(math.radians(30)) * unit_dist + z_offset, start_idx + 2],
        ["ATOM", start_idx + 3, "A", res, "", "LYS", "HZ3", f"A{res:04d}-LYS:HZ3", "H", "",
         x_offset, -math.sin(math.radians(30)) * unit_dist + y_offset, -math.cos(math.radians(30)) * unit_dist + z_offset, start_idx + 3],
        ["ATOM", start_idx + 4, "A", res, "", "LYS", "CE", f"A{res:04d}-LYS:CE", "C", "",
         -math.cos(math.radians(180 - 109.5)) * 1.07 - 1.47 + x_offset, y_offset, z_offset, start_idx + 4]
    ]
    hydrogen_mapping = [
        [f"A{res:04d}-LYS:HZ1", f"A{res:04d}-LYS:NZ", 1.07],
        [f"A{res:04d}-LYS:HZ2", f"A{res:04d}-LYS:NZ", 1.07],
        [f"A{res:04d}-LYS:HZ3", f"A{res:04d}-LYS:NZ", 1.07]
    ]
    covalent_edges = [
        (f"A{res:04d}-LYS:NZ", f"A{res:04d}-LYS:HZ1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-LYS:NZ", f"A{res:04d}-LYS:HZ2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-LYS:NZ", f"A{res:04d}-LYS:HZ3",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-LYS:CE", f"A{res:04d}-LYS:NZ",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.47, "base": True})
    ]
    pebble_edges = [
        (f"A{res:04d}-LYS:CE", f"A{res:04d}-LYS:NZ", {"weight": 0}),
        (f"A{res:04d}-LYS:NZ", f"A{res:04d}-LYS:CE", {"weight": 5})
    ]
    pebble_map = {
        f"A{res:04d}-LYS:NZ": 1,
        f"A{res:04d}-LYS:CE": 6,
    }
    pebble_bars = []

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )


def graph_cyh_acceptor(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        rotations: Optional[List[Tuple[float, _AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    coords = _modify_coords(_apply_rotations([
        ["ATOM", start_idx, "A", res, "", "CYH", "SG", f"A{res:04d}-CYH:SG", "S", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", "CYH", "HG", f"A{res:04d}-CYH:HG", "H", "",
         0, 1.41, 0, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", "CYH", "CB", f"A{res:04d}-CYH:CB", "C", "",
         +math.cos(math.radians(96 - 90)) * 1.81, -math.sin(math.radians(96 - 90)) * 1.81, 0, start_idx + 2]
    ], rotations), x_offset, y_offset, z_offset)
    # Cysteine bond angle of 96 degrees at sulphur atom, from:
    # Heyrovska, Raji (2011) : 'Precise Molecular Structures of Cysteine, Cystine, Hydrogen-Bonded Dicysteine, Cysteine Dipeptide, Glutathione and Acetyl Cysteine Based on Additivity of Atomic Radii'
    hydrogen_mapping = [
        (f"A{res:04d}-CYH:HG", f"A{res:04d}-CYH:SG", 1.41)
    ]
    covalent_edges = [
        (f"A{res:04d}-CYH:SG", f"A{res:04d}-CYH:HG",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.41, "base": True}),
        (f"A{res:04d}-CYH:CB", f"A{res:04d}-CYH:SG",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.81, "base": True})
    ]
    pebble_edges = [
        (f"A{res:04d}-CYH:CB", f"A{res:04d}-CYH:SG", {"weight": 0}),
        (f"A{res:04d}-CYH:SG", f"A{res:04d}-CYH:CB", {"weight": 5})
    ]
    pebble_map = {
        f"A{res:04d}-CYH:SG": 1,
        f"A{res:04d}-CYH:CB": 6
    }
    pebble_bars = []

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )


def graph_css(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        rotations: Optional[List[Tuple[float, _AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    coords = _modify_coords(_apply_rotations([
        ["ATOM", start_idx, "A", res, "", "CSS", "SG", f"A{res:04d}-CSS:SG", "S", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", "CSS", "CB", f"A{res:04d}-CSS:CB", "C", "",
         +math.cos(math.radians(96 - 90)) * 1.81, -math.sin(math.radians(96 - 90)) * 1.81, 0, start_idx + 1]
    ], rotations), x_offset, y_offset, z_offset)
    # Cysteine bond angle of 96 degrees at sulphur atom, from:
    # Heyrovska, Raji (2011) : 'Precise Molecular Structures of Cysteine, Cystine, Hydrogen-Bonded Dicysteine, Cysteine Dipeptide, Glutathione and Acetyl Cysteine Based on Additivity of Atomic Radii'
    hydrogen_mapping = []
    covalent_edges = [
        (f"A{res:04d}-CSS:CB", f"A{res:04d}-CSS:SG",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.81, "base": True})
    ]
    pebble_edges = [
        (f"A{res:04d}-CSS:CB", f"A{res:04d}-CSS:SG", {"weight": 0}),
        (f"A{res:04d}-CSS:SG", f"A{res:04d}-CSS:CB", {"weight": 5})
    ]
    pebble_map = {
        f"A{res:04d}-CSS:SG": 1,
        f"A{res:04d}-CSS:CB": 6
    }
    pebble_bars = []

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )


def graph_phe_aromatic_xy_plane(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        rotations: Optional[List[Tuple[float, _AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    c_dist = 1.39
    h_dist = c_dist + 1.14
    coords = _modify_coords(_apply_rotations([
        ["ATOM", start_idx, "A", res, "", "PHE", "CG", f"A{res:04d}-PHE:CG", "C", "",
         c_dist, 0, 0, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", "PHE", "CD1", f"A{res:04d}-PHE:CD1", "C", "",
         (c_dist * 0.5), -math.sin(math.radians(60)) * c_dist , 0, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", "PHE", "CD2", f"A{res:04d}-PHE:CD2", "C", "",
         -(c_dist * 0.5), -math.sin(math.radians(60)) * c_dist , 0, start_idx + 2],
        ["ATOM", start_idx + 3, "A", res, "", "PHE", "CE1", f"A{res:04d}-PHE:CE1", "C", "",
         (c_dist * 0.5), +math.sin(math.radians(60)) * c_dist , 0, start_idx + 3],
        ["ATOM", start_idx + 4, "A", res, "", "PHE", "CE2", f"A{res:04d}-PHE:CE2", "C", "",
         -(c_dist * 0.5), +math.sin(math.radians(60)) * c_dist , 0, start_idx + 4],
        ["ATOM", start_idx + 5, "A", res, "", "PHE", "CZ", f"A{res:04d}-PHE:CZ", "C", "",
         -c_dist, 0, 0, start_idx + 5],
        ["ATOM", start_idx + 6, "A", res, "", "PHE", "HD1", f"A{res:04d}-PHE:HD1", "H", "",
         (h_dist * 0.5), -math.sin(math.radians(60)) * h_dist , 0, start_idx + 6],
        ["ATOM", start_idx + 7, "A", res, "", "PHE", "HD2", f"A{res:04d}-PHE:HD2", "H", "",
         -(h_dist * 0.5), -math.sin(math.radians(60)) * h_dist , 0, start_idx + 7],
        ["ATOM", start_idx + 8, "A", res, "", "PHE", "HE1", f"A{res:04d}-PHE:HE1", "H", "",
         (h_dist * 0.5), +math.sin(math.radians(60)) * h_dist , 0, start_idx + 8],
        ["ATOM", start_idx + 9, "A", res, "", "PHE", "HE2", f"A{res:04d}-PHE:HE2", "H", "",
         -(h_dist * 0.5), +math.sin(math.radians(60)) * h_dist , 0, start_idx + 9],
        ["ATOM", start_idx + 10, "A", res, "", "PHE", "HZ", f"A{res:04d}-PHE:HZ", "H", "",
         -h_dist, 0, 0, start_idx + 10]
    ], rotations), x_offset, y_offset, z_offset)
    hydrogen_mapping = []
    covalent_edges = [
        (f"A{res:04d}-PHE:CG", f"A{res:04d}-PHE:CD1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.39, "base": True}),
        (f"A{res:04d}-PHE:CG", f"A{res:04d}-PHE:CD2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.39, "base": True}),
        (f"A{res:04d}-PHE:CD1", f"A{res:04d}-PHE:CE1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.39, "base": True}),
        (f"A{res:04d}-PHE:CD2", f"A{res:04d}-PHE:CE2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.39, "base": True}),
        (f"A{res:04d}-PHE:CE1", f"A{res:04d}-PHE:CZ",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.39, "base": True}),
        (f"A{res:04d}-PHE:CE2", f"A{res:04d}-PHE:CZ",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.39, "base": True}),
        (f"A{res:04d}-PHE:CD1", f"A{res:04d}-PHE:HD1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-PHE:CD2", f"A{res:04d}-PHE:HD2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-PHE:CE1", f"A{res:04d}-PHE:HE1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-PHE:CE2", f"A{res:04d}-PHE:HE2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-PHE:CZ", f"A{res:04d}-PHE:HZ",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True})
    ]
    pebble_edges = [
        (f"A{res:04d}-PHE:CG", f"A{res:04d}-PHE:CD1", {"weight": 0}),
        (f"A{res:04d}-PHE:CD1", f"A{res:04d}-PHE:CG", {"weight": 5}),
        (f"A{res:04d}-PHE:CG", f"A{res:04d}-PHE:CD2", {"weight": 0}),
        (f"A{res:04d}-PHE:CD2", f"A{res:04d}-PHE:CG", {"weight": 5}),
        (f"A{res:04d}-PHE:CD1", f"A{res:04d}-PHE:CE1", {"weight": 0}),
        (f"A{res:04d}-PHE:CE1", f"A{res:04d}-PHE:CD1", {"weight": 5}),
        (f"A{res:04d}-PHE:CD2", f"A{res:04d}-PHE:CE2", {"weight": 0}),
        (f"A{res:04d}-PHE:CE2", f"A{res:04d}-PHE:CD2", {"weight": 5}),
        (f"A{res:04d}-PHE:CE1", f"A{res:04d}-PHE:CZ", {"weight": 0}),
        (f"A{res:04d}-PHE:CZ", f"A{res:04d}-PHE:CE1", {"weight": 5})
    ]
    pebble_map = {
        f"A{res:04d}-PHE:CZ": 1,
        f"A{res:04d}-PHE:CE2": 1,
        f"A{res:04d}-PHE:CE1": 1,
        f"A{res:04d}-PHE:CD1": 1,
        f"A{res:04d}-PHE:CD2": 1,
        f"A{res:04d}-PHE:CG": 6
    }
    pebble_bars = [
        (f"A{res:04d}-PHE:CZ", f"A{res:04d}-PHE:CE2", 5)
    ]

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )


def graph_hie_imidazole_ring(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        rotations: Optional[List[Tuple[float, _AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    dist_cg_nd1 = 1.382
    angle_cg = 108.1
    dist_cg_cd2 = 1.359
    angle_nd1 = 106.5
    dist_nd1_ce1 = 1.324
    angle_cd2 = 106.9
    dist_ce1_ne2 = 1.334
    angle_ce1 = 110.7
    dist_cd2_ne2 = 1.367
    angle_ne2 = 107.8

    # NE2 at origin, imidazole ring in -x direction
    half_angle_ne2 = angle_ne2 / 2
    counter_angle_half_ne2 = 90 - half_angle_ne2

    adjusted_angle_cd2 = angle_cd2 - counter_angle_half_ne2
    adjusted_angle_ce1 = angle_ce1 - counter_angle_half_ne2

    ce1_x = -math.cos(math.radians(half_angle_ne2)) * dist_ce1_ne2
    ce1_y = +math.sin(math.radians(half_angle_ne2)) * dist_ce1_ne2

    cd2_x = -math.cos(math.radians(half_angle_ne2)) * dist_cd2_ne2
    cd2_y = -math.sin(math.radians(half_angle_ne2)) * dist_cd2_ne2

    nd1_x = ce1_x - math.sin(math.radians(adjusted_angle_ce1)) * dist_nd1_ce1
    nd1_y = ce1_y - math.cos(math.radians(adjusted_angle_ce1)) * dist_nd1_ce1

    cg_x = cd2_x - math.sin(math.radians(adjusted_angle_cd2)) * dist_cg_cd2
    cg_y = cd2_y + math.cos(math.radians(adjusted_angle_cd2)) * dist_cg_cd2

    adjusted_angle_hd2 = 180 - angle_cd2/2 - half_angle_ne2

    hd2_x = cd2_x + math.cos(math.radians(adjusted_angle_hd2)) * 1.14
    hd2_y = cd2_y - math.sin(math.radians(adjusted_angle_hd2)) * 1.14

    adjusted_angle_he1 = 180 - angle_ce1/2 - half_angle_ne2

    he1_x = ce1_x + math.cos(math.radians(adjusted_angle_he1)) * 1.14
    he1_y = ce1_y + math.sin(math.radians(adjusted_angle_he1)) * 1.14

    coords = _modify_coords(_apply_rotations([
        ["ATOM", start_idx, "A", res, "", "HIE", "CG", f"A{res:04d}-HIE:CG", "C", "",
         cg_x, cg_y, 0, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", "HIE", "ND1", f"A{res:04d}-HIE:ND1", "N", "",
         nd1_x, nd1_y, 0, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", "HIE", "CD2", f"A{res:04d}-HIE:CD2", "C", "",
         cd2_x, cd2_y, 0, start_idx + 2],
        ["ATOM", start_idx + 3, "A", res, "", "HIE", "CE1", f"A{res:04d}-HIE:CE1", "C", "",
         ce1_x, ce1_y, 0, start_idx + 3],
        ["ATOM", start_idx + 4, "A", res, "", "HIE", "NE2", f"A{res:04d}-HIE:NE2", "N", "",
         0, 0, 0, start_idx + 4],
        ["ATOM", start_idx + 5, "A", res, "", "HIE", "HD2", f"A{res:04d}-HIE:HD2", "H", "",
         hd2_x, hd2_y, 0, start_idx + 5],
        ["ATOM", start_idx + 6, "A", res, "", "HIE", "HE1", f"A{res:04d}-HIE:HE1", "H", "",
         he1_x, he1_y, 0, start_idx + 6],
        ["ATOM", start_idx + 7, "A", res, "", "HIE", "HE2", f"A{res:04d}-HIE:HE2", "H", "",
         1.07, 0, 0, start_idx + 7],
    ], rotations), x_offset, y_offset, z_offset)
    hydrogen_mapping = [
        (f"A{res:04d}-HIE:HD2", f"A{res:04d}-HIE:CD2", 1.14),
        (f"A{res:04d}-HIE:HE1", f"A{res:04d}-HIE:CE1", 1.14),
        (f"A{res:04d}-HIE:HE2", f"A{res:04d}-HIE:NE2", 1.07)
    ]
    covalent_edges = [
        (f"A{res:04d}-HIE:CG", f"A{res:04d}-HIE:ND1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": dist_cg_nd1, "base": True}),
        (f"A{res:04d}-HIE:CG", f"A{res:04d}-HIE:CD2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": dist_cg_cd2, "base": True}),
        (f"A{res:04d}-HIE:ND1", f"A{res:04d}-HIE:CE1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": dist_nd1_ce1, "base": True}),
        (f"A{res:04d}-HIE:CD2", f"A{res:04d}-HIE:NE2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": dist_cd2_ne2, "base": True}),
        (f"A{res:04d}-HIE:CE1", f"A{res:04d}-HIE:NE2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": dist_ce1_ne2, "base": True}),
        (f"A{res:04d}-HIE:CD2", f"A{res:04d}-HIE:HD2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-HIE:CE1", f"A{res:04d}-HIE:HE1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-HIE:NE2", f"A{res:04d}-HIE:HE2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True})
    ]
    pebble_edges = [
        (f"A{res:04d}-HIE:CG", f"A{res:04d}-HIE:ND1", {"weight": 0}),
        (f"A{res:04d}-HIE:ND1", f"A{res:04d}-HIE:CG", {"weight": 5}),
        (f"A{res:04d}-HIE:CG", f"A{res:04d}-HIE:CD2", {"weight": 0}),
        (f"A{res:04d}-HIE:CD2", f"A{res:04d}-HIE:CG", {"weight": 5}),
        (f"A{res:04d}-HIE:ND1", f"A{res:04d}-HIE:CE1", {"weight": 0}),
        (f"A{res:04d}-HIE:CE1", f"A{res:04d}-HIE:ND1", {"weight": 5}),
        (f"A{res:04d}-HIE:CD2", f"A{res:04d}-HIE:NE2", {"weight": 0}),
        (f"A{res:04d}-HIE:NE2", f"A{res:04d}-HIE:CD2", {"weight": 5})
    ]
    pebble_map = {
        f"A{res:04d}-HIE:NE2": 1,
        f"A{res:04d}-HIE:CE1": 1,
        f"A{res:04d}-HIE:ND1": 1,
        f"A{res:04d}-HIE:CD2": 1,
        f"A{res:04d}-HIE:CG": 6
    }
    pebble_bars = [
        (f"A{res:04d}-HIE:CE1", f"A{res:04d}-HIE:NE2", 5)
    ]

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )

def graph_methane(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        rotations: Optional[List[Tuple[float, _AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    d = 0.658  # d = 0.658179 => sqrt(3 * d^2) == 1.14, however, PDB only allows 3 digits
    coords = _modify_coords(_apply_rotations([
        ["HETATM", start_idx, "A", res, "", "CH4", "C", f"A{res:04d}-CH4:C", "C", "",
         0, 0, 0, start_idx],
        ["HETATM", start_idx + 1, "A", res, "", "CH4", "H1", f"A{res:04d}-CH4:H1", "H", "",
         d, -d, d, start_idx + 1],
        ["HETATM", start_idx + 2, "A", res, "", "CH4", "H2", f"A{res:04d}-CH4:H2", "H", "",
         d, d, -d, start_idx + 2],
        ["HETATM", start_idx + 3, "A", res, "", "CH4", "H3", f"A{res:04d}-CH4:H3", "H", "",
         -d, d, d, start_idx + 3],
        ["HETATM", start_idx + 4, "A", res, "", "CH4", "H4", f"A{res:04d}-CH4:H4", "H", "",
         -d, -d, -d, start_idx + 4],
    ], rotations), x_offset, y_offset, z_offset)
    hydrogen_mapping = [
        (f"A{res:04d}-CH4:H1", f"A{res:04d}-CH4:C", 1.14),
        (f"A{res:04d}-CH4:H2", f"A{res:04d}-CH4:C", 1.14),
        (f"A{res:04d}-CH4:H3", f"A{res:04d}-CH4:C", 1.14),
        (f"A{res:04d}-CH4:H4", f"A{res:04d}-CH4:C", 1.14)
    ]
    covalent_edges = [
        (f"A{res:04d}-CH4:C", f"A{res:04d}-CH4:H1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-CH4:C", f"A{res:04d}-CH4:H2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-CH4:C", f"A{res:04d}-CH4:H3",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-CH4:C", f"A{res:04d}-CH4:H4",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True})
    ]
    pebble_edges = []
    pebble_map = {
        f"A{res:04d}-CH4:C": 6,
        f"A{res:04d}-CH4:H1": 6,
        f"A{res:04d}-CH4:H2": 6,
        f"A{res:04d}-CH4:H3": 6,
        f"A{res:04d}-CH4:H4": 6
    }
    pebble_bars = [
        (f"A{res:04d}-CH4:C", f"A{res:04d}-CH4:H1", 5),
        (f"A{res:04d}-CH4:C", f"A{res:04d}-CH4:H2", 5),
        (f"A{res:04d}-CH4:C", f"A{res:04d}-CH4:H3", 5),
        (f"A{res:04d}-CH4:C", f"A{res:04d}-CH4:H4", 5)
    ]

    return _format(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )
