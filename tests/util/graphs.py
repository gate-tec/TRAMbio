import math
from typing import Optional, List, Tuple

from TRAMbio.util.constants.interaction import InteractionType
from tests.util.protein_graph_utils import *


def graph_gly_donor(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    coords = modify_coords(apply_rotations([
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

    return format_coords(
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
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    coords = modify_coords(apply_rotations([
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

    return format_coords(
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

    return format_coords(
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

    return format_coords(
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
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    coords = modify_coords(apply_rotations([
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

    return format_coords(
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

    return format_coords(
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

    return format_coords(
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
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    coords = modify_coords(apply_rotations([
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

    return format_coords(
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
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    coords = modify_coords(apply_rotations([
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

    return format_coords(
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
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    c_dist = 1.39
    h_dist = c_dist + 1.14
    coords = modify_coords(apply_rotations([
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

    return format_coords(
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
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
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

    coords = modify_coords(apply_rotations([
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

    return format_coords(
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
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    d = 0.658  # d = 0.658179 => sqrt(3 * d^2) == 1.14, however, PDB only allows 3 digits
    coords = modify_coords(apply_rotations([
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

    return format_coords(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )

def graph_adenine(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        is_dna: bool = True,
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    resi = ("DA" if is_dna else "A").rjust(3, " ")
    coords = modify_coords(apply_rotations([
        ["ATOM", start_idx, "A", res, "", resi, "P", f"A{res:04d}-{resi}:P", "P", "",
         1.861, -8.579, 1.742, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", resi, "OP1", f"A{res:04d}-{resi}:OP1", "O", "",
         1.588, -9.922, 2.372, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", resi, "OP2", f"A{res:04d}-{resi}:OP2", "O", "",
         3.084, -8.504, 0.917, start_idx + 2],
        ["ATOM", start_idx + 3, "A", res, "", resi, "O5'", f"A{res:04d}-{resi}:O5'", "O", "",
         0.435, -8.093, 1.219, start_idx + 3],
        ["ATOM", start_idx + 4, "A", res, "", resi, "C5'", f"A{res:04d}-{resi}:C5'", "C", "",
         -0.730, -8.036, 2.041, start_idx + 4],
        ["ATOM", start_idx + 5, "A", res, "", resi, "C4'", f"A{res:04d}-{resi}:C4'", "C", "",
         -1.602, -6.896, 1.546, start_idx + 5],
        ["ATOM", start_idx + 6, "A", res, "", resi, "O4'", f"A{res:04d}-{resi}:O4'", "O", "",
         -0.904, -5.683, 1.553, start_idx + 6],
        ["ATOM", start_idx + 7, "A", res, "", resi, "C3'", f"A{res:04d}-{resi}:C3'", "C", "",
         -2.233, -7.042, 0.172, start_idx + 7],
        ["ATOM", start_idx + 8, "A", res, "", resi, "O3'", f"A{res:04d}-{resi}:O3'", "O", "",
         -3.658, -7.083, 0.244, start_idx + 8],
        ["ATOM", start_idx + 9, "A", res, "", resi, "C2'", f"A{res:04d}-{resi}:C2'", "C", "",
         -1.653, -5.866, -0.623, start_idx + 9],
        ["ATOM", start_idx + 10, "A", res, "", resi, "C1'", f"A{res:04d}-{resi}:C1'", "C", "",
         -1.281, -4.864, 0.400, start_idx + 10],
        ["ATOM", start_idx + 11, "A", res, "", resi, "N9", f"A{res:04d}-{resi}:N9", "N", "",
         -0.163, -3.934, 0.161, start_idx + 11],
        ["ATOM", start_idx + 12, "A", res, "", resi, "C8", f"A{res:04d}-{resi}:C8", "C", "",
         1.183, -4.183, 0.043, start_idx + 12],
        ["ATOM", start_idx + 13, "A", res, "", resi, "N7", f"A{res:04d}-{resi}:N7", "N", "",
         1.919, -3.104, -0.094, start_idx + 13],
        ["ATOM", start_idx + 14, "A", res, "", resi, "C5", f"A{res:04d}-{resi}:C5", "C", "",
         0.997, -2.064, -0.032, start_idx + 14],
        ["ATOM", start_idx + 15, "A", res, "", resi, "C6", f"A{res:04d}-{resi}:C6", "C", "",
         1.155, -0.652, -0.114, start_idx + 15],
        ["ATOM", start_idx + 16, "A", res, "", resi, "N6", f"A{res:04d}-{resi}:N6", "N", "",
         2.306, -0.014, -0.284, start_idx + 16],
        ["ATOM", start_idx + 17, "A", res, "", resi, "N1", f"A{res:04d}-{resi}:N1", "N", "",
         0.000, 0.082, 0.000, start_idx + 17],
        ["ATOM", start_idx + 18, "A", res, "", resi, "C2", f"A{res:04d}-{resi}:C2", "C", "",
         -1.203, -0.539, 0.168, start_idx + 18],
        ["ATOM", start_idx + 19, "A", res, "", resi, "N3", f"A{res:04d}-{resi}:N3", "N", "",
         -1.421, -1.847, 0.240, start_idx + 19],
        ["ATOM", start_idx + 20, "A", res, "", resi, "C4", f"A{res:04d}-{resi}:C4", "C", "",
         -0.278, -2.554, 0.133, start_idx + 20],
        ["ATOM", start_idx + 21, "A", res, "", resi, "H5'", f"A{res:04d}-{resi}:H5'", "H", "",
         -0.483, -7.896, 2.969, start_idx + 21],
        ["ATOM", start_idx + 22, "A", res, "", resi, "H5''", f"A{res:04d}-{resi}:H5''", "H", "",
         -1.214, -8.876, 2.000, start_idx + 22],
        ["ATOM", start_idx + 23, "A", res, "", resi, "H4'", f"A{res:04d}-{resi}:H4'", "H", "",
         -2.333, -6.917, 2.184, start_idx + 23],
        ["ATOM", start_idx + 24, "A", res, "", resi, "H3'", f"A{res:04d}-{resi}:H3'", "H", "",
         -2.027, -7.882, -0.267, start_idx + 24],
        ["ATOM", start_idx + 25, "A", res, "", resi, "H2'", f"A{res:04d}-{resi}:H2'", "H", "",
         -0.882, -6.140, -1.144, start_idx + 25],
        ["ATOM", start_idx + 26, "A", res, "", resi, "H1'", f"A{res:04d}-{resi}:H1'", "H", "",
         -2.040, -4.263, 0.465, start_idx + 26],
        ["ATOM", start_idx + 27, "A", res, "", resi, "H8", f"A{res:04d}-{resi}:H8", "H", "",
         1.542, -5.040, 0.060, start_idx + 27],
        ["ATOM", start_idx + 28, "A", res, "", resi, "H61", f"A{res:04d}-{resi}:H61", "H", "",
         2.320, 0.845, -0.324, start_idx + 28],
        ["ATOM", start_idx + 29, "A", res, "", resi, "H62", f"A{res:04d}-{resi}:H62", "H", "",
         3.038, -0.460, -0.355, start_idx + 29],
        ["ATOM", start_idx + 30, "A", res, "", resi, "H2", f"A{res:04d}-{resi}:H2", "H", "",
         -1.950, 0.012, 0.240, start_idx + 30],
    ] + ([
        ["ATOM", start_idx + 31, "A", res, "", resi, "H2''", f"A{res:04d}-{resi}:H2''", "H", "",
         -2.305, -5.508, -1.246, start_idx + 31],
    ] if is_dna else [
        ["ATOM", start_idx + 31, "A", res, "", resi, "O2'", f"A{res:04d}-{resi}:O2'", "O", "",
         0, 0, 0, start_idx + 31],
        ["ATOM", start_idx + 32, "A", res, "", resi, "HO2'", f"A{res:04d}-{resi}:HO2'", "H", "",
         0, 0, 0, start_idx + 32],
    ]), rotations), x_offset, y_offset, z_offset)
    hydrogen_mapping = [
        (f"A{res:04d}-{resi}:H8", f"A{res:04d}-{resi}:C8", 0.930),
        (f"A{res:04d}-{resi}:H61", f"A{res:04d}-{resi}:C6", 1.909),
        (f"A{res:04d}-{resi}:H62", f"A{res:04d}-{resi}:C6", 1.908),
        (f"A{res:04d}-{resi}:H2", f"A{res:04d}-{resi}:C2", 0.930),
        # Backbone
        (f"A{res:04d}-{resi}:H5'", f"A{res:04d}-{resi}:C5'", 0.970),
        (f"A{res:04d}-{resi}:H5''", f"A{res:04d}-{resi}:C5'", 0.970),
        (f"A{res:04d}-{resi}:H4'", f"A{res:04d}-{resi}:C4'", 0.970),
        (f"A{res:04d}-{resi}:H3'", f"A{res:04d}-{resi}:C3'", 0.970),
        (f"A{res:04d}-{resi}:H2'", f"A{res:04d}-{resi}:C2'", 0.970),
        (f"A{res:04d}-{resi}:H1'", f"A{res:04d}-{resi}:C1'", 0.970),
    ] + ([
        (f"A{res:04d}-{resi}:H2''", f"A{res:04d}-{resi}:C2'", 0.970),
    ] if is_dna else [
        (f"A{res:04d}-{resi}:HO2'", f"A{res:04d}-{resi}:O2'", 1.04),
    ])
    covalent_edges = [
        (f"A{res:04d}-{resi}:OP1", f"A{res:04d}-{resi}:P",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.508, "base": True}),
        (f"A{res:04d}-{resi}:OP2", f"A{res:04d}-{resi}:P",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.477, "base": True}),
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:O5'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.595, "base": True}),
        (f"A{res:04d}-{resi}:O5'", f"A{res:04d}-{resi}:C5'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.427, "base": True}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:H5'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:H5''",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:C4'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.519, "base": True}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:H4'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C4'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.400, "base": True}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:O3'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 2.441, "base": True}),
        (f"A{res:04d}-{resi}:O3'", f"A{res:04d}-{resi}:C3'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.427, "base": True}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:H3'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:C2'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.534, "base": True}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:H2'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:C1'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.479, "base": True}),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C1'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.464, "base": True}),
        (f"A{res:04d}-{resi}:C1'", f"A{res:04d}-{resi}:H1'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
        (f"A{res:04d}-{resi}:N1", f"A{res:04d}-{resi}:C2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.364, "base": True}),
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:N3",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.328, "base": True}),
        (f"A{res:04d}-{resi}:N3", f"A{res:04d}-{resi}:C4",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.348, "base": True}),
        (f"A{res:04d}-{resi}:N6", f"A{res:04d}-{resi}:C6",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.327, "base": True}),
        (f"A{res:04d}-{resi}:C6", f"A{res:04d}-{resi}:C5",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.424, "base": True}),
        (f"A{res:04d}-{resi}:C5", f"A{res:04d}-{resi}:C4",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.376, "base": True}),
        (f"A{res:04d}-{resi}:C4", f"A{res:04d}-{resi}:N9",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.385, "base": True}),
        (f"A{res:04d}-{resi}:N7", f"A{res:04d}-{resi}:C8",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.312, "base": True}),
        (f"A{res:04d}-{resi}:C8", f"A{res:04d}-{resi}:N9",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.374, "base": True}),
        (f"A{res:04d}-{resi}:N9", f"A{res:04d}-{resi}:C1'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.474, "base": True}),
        (f"A{res:04d}-{resi}:C8", f"A{res:04d}-{resi}:H8",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.930, "base": True}),
        (f"A{res:04d}-{resi}:C6", f"A{res:04d}-{resi}:H61",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.909, "base": True}),
        (f"A{res:04d}-{resi}:C6", f"A{res:04d}-{resi}:H62",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.908, "base": True}),
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:H2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.930, "base": True}),
        (f"A{res:04d}-{resi}:C6", f"A{res:04d}-{resi}:N1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.373, "base": True}),
        (f"A{res:04d}-{resi}:C5", f"A{res:04d}-{resi}:N7",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.391, "base": True}),
    ] + ([
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:H2''",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
    ] if is_dna else [
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:O2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.44, "base": True}),
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:HO2'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.04, "base": True}),
    ])
    pebble_edges = [
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:OP1", {"weight": 0}),
        (f"A{res:04d}-{resi}:OP1", f"A{res:04d}-{resi}:P", {"weight": 5}),
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:OP2", {"weight": 0}),
        (f"A{res:04d}-{resi}:OP2", f"A{res:04d}-{resi}:P", {"weight": 5}),
        (f"A{res:04d}-{resi}:O5'", f"A{res:04d}-{resi}:P", {"weight": 0}),
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:O5'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:O5'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O5'", f"A{res:04d}-{resi}:C5'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:C5'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:C4'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:O4'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C4'", {"weight": 5}),
        (f"A{res:04d}-{resi}:O3'", f"A{res:04d}-{resi}:C4'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:O3'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:O3'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O3'", f"A{res:04d}-{resi}:C3'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:C3'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:C2'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C1'", f"A{res:04d}-{resi}:C2'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:C1'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:N1", {"weight": 0}),
        (f"A{res:04d}-{resi}:N1", f"A{res:04d}-{resi}:C2", {"weight": 5}),
        (f"A{res:04d}-{resi}:N3", f"A{res:04d}-{resi}:C2", {"weight": 0}),
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:N3", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4", f"A{res:04d}-{resi}:N3", {"weight": 0}),
        (f"A{res:04d}-{resi}:N3", f"A{res:04d}-{resi}:C4", {"weight": 5}),
        (f"A{res:04d}-{resi}:C6", f"A{res:04d}-{resi}:N6", {"weight": 0}),
        (f"A{res:04d}-{resi}:N6", f"A{res:04d}-{resi}:C6", {"weight": 5}),
        (f"A{res:04d}-{resi}:C5", f"A{res:04d}-{resi}:C6", {"weight": 0}),
        (f"A{res:04d}-{resi}:C6", f"A{res:04d}-{resi}:C5", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4", f"A{res:04d}-{resi}:C5", {"weight": 0}),
        (f"A{res:04d}-{resi}:C5", f"A{res:04d}-{resi}:C4", {"weight": 5}),
        (f"A{res:04d}-{resi}:N9", f"A{res:04d}-{resi}:C4", {"weight": 0}),
        (f"A{res:04d}-{resi}:C4", f"A{res:04d}-{resi}:N9", {"weight": 5}),
        (f"A{res:04d}-{resi}:C8", f"A{res:04d}-{resi}:N7", {"weight": 0}),
        (f"A{res:04d}-{resi}:N7", f"A{res:04d}-{resi}:C8", {"weight": 5}),
        (f"A{res:04d}-{resi}:N9", f"A{res:04d}-{resi}:C8", {"weight": 0}),
        (f"A{res:04d}-{resi}:C8", f"A{res:04d}-{resi}:N9", {"weight": 5}),
        (f"A{res:04d}-{resi}:C1'", f"A{res:04d}-{resi}:N9", {"weight": 0}),
        (f"A{res:04d}-{resi}:N9", f"A{res:04d}-{resi}:C1'", {"weight": 5}),
    ] + ([
        # not required
    ] if is_dna else [
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:O2'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:C2'", {"weight": 5}),
    ])
    pebble_map = {
        f"A{res:04d}-{resi}:N1": 1,
        # Backbone
        f"A{res:04d}-{resi}:OP1": 1,
        f"A{res:04d}-{resi}:OP2": 1,
        f"A{res:04d}-{resi}:P": 1,
        f"A{res:04d}-{resi}:O5'": 1,
        f"A{res:04d}-{resi}:C5'": 1,
        f"A{res:04d}-{resi}:O4'": 1,
        f"A{res:04d}-{resi}:C4'": 1,
        f"A{res:04d}-{resi}:C3'": 1,
        f"A{res:04d}-{resi}:O3'": 6,
        f"A{res:04d}-{resi}:C2'": 1,
        f"A{res:04d}-{resi}:C1'": 1,
    }
    if not is_dna:
        pebble_map = dict(**pebble_map, **{
            f"A{res:04d}-{resi}:O2'": 1,
        })
    pebble_bars = [
        # Backbone
        (f"A{res:04d}-{resi}:OP1", f"A{res:04d}-{resi}:P", 1),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C1'", 5),
    ] + ([
        # not required
    ] if is_dna else [
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:C2'", 1),
    ])

    return format_coords(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )

def graph_cytosine(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        is_dna: bool = True,
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    resi = ("DC" if is_dna else "C").rjust(3, " ")
    coords = modify_coords(apply_rotations([
        # Nucleotide heavy atoms
        ["ATOM", start_idx, "A", res, "", resi, "N3", f"A{res:04d}-{resi}:N3", "N", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C2", f"A{res:04d}-{resi}:C2", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "O2", f"A{res:04d}-{resi}:O2", "O", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "N1", f"A{res:04d}-{resi}:N1", "N", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C4", f"A{res:04d}-{resi}:C4", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "N4", f"A{res:04d}-{resi}:N4", "N", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C5", f"A{res:04d}-{resi}:C5", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C6", f"A{res:04d}-{resi}:C6", "C", "",
         0, 0, 0, start_idx],
        # Nucleotide hydrogen atoms
        ["ATOM", start_idx, "A", res, "", resi, "H6", f"A{res:04d}-{resi}:H6", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H5", f"A{res:04d}-{resi}:H5", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H41", f"A{res:04d}-{resi}:H41", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H42", f"A{res:04d}-{resi}:H42", "H", "",
         0, 0, 0, start_idx],
        # Backbone
        ["ATOM", start_idx, "A", res, "", resi, "C1'", f"A{res:04d}-{resi}:C1'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C2'", f"A{res:04d}-{resi}:C2'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C3'", f"A{res:04d}-{resi}:C3'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "O3'", f"A{res:04d}-{resi}:O3'", "O", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C4'", f"A{res:04d}-{resi}:C4'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "O4'", f"A{res:04d}-{resi}:O4'", "O", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C5'", f"A{res:04d}-{resi}:C5'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "O5'", f"A{res:04d}-{resi}:O5'", "O", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C4'", f"A{res:04d}-{resi}:C2'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H1'", f"A{res:04d}-{resi}:H1'", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H3'", f"A{res:04d}-{resi}:H3'", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H4'", f"A{res:04d}-{resi}:H4'", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H5'", f"A{res:04d}-{resi}:H5'", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H5''", f"A{res:04d}-{resi}:H5''", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "P", f"A{res:04d}-{resi}:P", "P", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "OP1'", f"A{res:04d}-{resi}:OP1", "O", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "OP2", f"A{res:04d}-{resi}:OP2'", "O", "",
         0, 0, 0, start_idx],
    ] + ([
        ["ATOM", start_idx, "A", res, "", resi, "H2''", f"A{res:04d}-{resi}:H2''", "H", "",
         0, 0, 0, start_idx],
    ] if is_dna else [
        ["ATOM", start_idx, "A", res, "", resi, "O2'", f"A{res:04d}-{resi}:O2'", "O", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "HO2'", f"A{res:04d}-{resi}:HO2'", "H", "",
         0, 0, 0, start_idx],
    ]), rotations), x_offset, y_offset, z_offset)
    hydrogen_mapping = [
        (f"A{res:04d}-{resi}:H6", f"A{res:04d}-{resi}:C6", 1.09),
        (f"A{res:04d}-{resi}:H5", f"A{res:04d}-{resi}:C5", 1.09),
        (f"A{res:04d}-{resi}:H41", f"A{res:04d}-{resi}:N4", 1.07),
        (f"A{res:04d}-{resi}:H42", f"A{res:04d}-{resi}:N4", 1.07),
        # Backbone
        (f"A{res:04d}-{resi}:H5'", f"A{res:04d}-{resi}:C5'", 1.14),
        (f"A{res:04d}-{resi}:H5''", f"A{res:04d}-{resi}:C5'", 1.14),
        (f"A{res:04d}-{resi}:H4'", f"A{res:04d}-{resi}:C4'", 1.14),
        (f"A{res:04d}-{resi}:H3'", f"A{res:04d}-{resi}:C3'", 1.14),
        (f"A{res:04d}-{resi}:H2'", f"A{res:04d}-{resi}:C2'", 1.14),
        (f"A{res:04d}-{resi}:H1'", f"A{res:04d}-{resi}:C1'", 1.14),
    ] + ([
        (f"A{res:04d}-{resi}:H2''", f"A{res:04d}-{resi}:C2'", 1.14),
    ] if is_dna else [
        (f"A{res:04d}-{resi}:HO2'", f"A{res:04d}-{resi}:O2'", 1.04),
    ])
    covalent_edges = [
        (f"A{res:04d}-{resi}:N1", f"A{res:04d}-{resi}:C2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        # Backbone
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:H5''",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:H5'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:H5''",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:H4'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:H3'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:H2'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C1'", f"A{res:04d}-{resi}:H1'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C8", f"A{res:04d}-{resi}:H8",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.09, "base": True}),
        (f"A{res:04d}-{resi}:N6", f"A{res:04d}-{resi}:H61",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-{resi}:N6", f"A{res:04d}-{resi}:H62",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:H2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.04, "base": True}),
    ] + ([
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:H2''",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
    ] if is_dna else [
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:O2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.44, "base": True}),
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:HO2'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.04, "base": True}),
    ])
    pebble_edges = [
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:N1", {"weight": 0}),
        (f"A{res:04d}-{resi}:N1", f"A{res:04d}-{resi}:C2", {"weight": 5}),
        # Backbone
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:OP1", {"weight": 0}),
        (f"A{res:04d}-{resi}:OP1", f"A{res:04d}-{resi}:P", {"weight": 5}),
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:OP2", {"weight": 0}),
        (f"A{res:04d}-{resi}:OP2", f"A{res:04d}-{resi}:P", {"weight": 5}),
        (f"A{res:04d}-{resi}:O5'", f"A{res:04d}-{resi}:P", {"weight": 0}),
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:O5'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:O5'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O5'", f"A{res:04d}-{resi}:C5'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:C5'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:C4'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:O4'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C4'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:C4'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:C3'", {"weight": 5}),
        (f"A{res:04d}-{resi}:O3'", f"A{res:04d}-{resi}:C3'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:O3'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:C2'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:C3'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:C1'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C1'", f"A{res:04d}-{resi}:C2'", {"weight": 5}),
    ] + ([
        # not required
    ] if is_dna else [
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:O2'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:C2'", {"weight": 5}),
    ])
    pebble_map = {
        f"A{res:04d}-{resi}:N1": 1,
        # Backbone
        f"A{res:04d}-{resi}:OP1": 1,
        f"A{res:04d}-{resi}:OP2": 1,
        f"A{res:04d}-{resi}:P": 1,
        f"A{res:04d}-{resi}:O5'": 1,
        f"A{res:04d}-{resi}:C5'": 1,
        f"A{res:04d}-{resi}:O4'": 1,
        f"A{res:04d}-{resi}:C4'": 1,
        f"A{res:04d}-{resi}:C3'": 1,
        f"A{res:04d}-{resi}:O3'": 6,
        f"A{res:04d}-{resi}:C2'": 1,
        f"A{res:04d}-{resi}:C1'": 1,
    }
    if not is_dna:
        pebble_map = dict(**pebble_map, **{
            f"A{res:04d}-{resi}:O2'": 1,
        })
    pebble_bars = [
        # Backbone
        (f"A{res:04d}-{resi}:OP1", f"A{res:04d}-{resi}:P", 1),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C1'", 5),
    ] + ([
        # not required
    ] if is_dna else [
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:C2'", 1),
    ])

    return format_coords(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )

def graph_guanine(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        is_dna: bool = True,
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    resi = ("DG" if is_dna else "G").rjust(3, " ")
    coords = modify_coords(apply_rotations([
        ["ATOM", start_idx, "A", res, "", resi, "P", f"A{res:04d}-{resi}:P", "P", "",
         -1.795, -8.507, 2.382, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", resi, "OP1", f"A{res:04d}-{resi}:OP1", "O", "",
         -1.219, -9.865, 2.670, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", resi, "OP2", f"A{res:04d}-{resi}:OP2", "O", "",
         -3.213, -8.482, 1.933, start_idx + 2],
        ["ATOM", start_idx + 3, "A", res, "", resi, "O5'", f"A{res:04d}-{resi}:O5'", "O", "",
         -0.807, -7.765, 1.316, start_idx + 3],
        ["ATOM", start_idx + 4, "A", res, "", resi, "C5'", f"A{res:04d}-{resi}:C5'", "C", "",
         0.541, -8.373, 1.203, start_idx + 4],
        ["ATOM", start_idx + 5, "A", res, "", resi, "C4'", f"A{res:04d}-{resi}:C4'", "C", "",
         1.383, -7.286, 0.543, start_idx + 5],
        ["ATOM", start_idx + 6, "A", res, "", resi, "O4'", f"A{res:04d}-{resi}:O4'", "O", "",
         0.885, -6.031, 0.997, start_idx + 6],
        ["ATOM", start_idx + 7, "A", res, "", resi, "C3'", f"A{res:04d}-{resi}:C3'", "C", "",
         1.115, -7.262, -0.894, start_idx + 7],
        ["ATOM", start_idx + 8, "A", res, "", resi, "O3'", f"A{res:04d}-{resi}:O3'", "O", "",
         2.281, -7.673, -1.447, start_idx + 8],
        ["ATOM", start_idx + 9, "A", res, "", resi, "C2'", f"A{res:04d}-{resi}:C2'", "C", "",
         0.848, -5.864, -1.184, start_idx + 9],
        ["ATOM", start_idx + 10, "A", res, "", resi, "C1'", f"A{res:04d}-{resi}:C1'", "C", "",
         0.985, -5.083, -0.036, start_idx + 10],
        ["ATOM", start_idx + 11, "A", res, "", resi, "N9", f"A{res:04d}-{resi}:N9", "N", "",
         -0.070, -4.061, -0.119, start_idx + 11],
        ["ATOM", start_idx + 12, "A", res, "", resi, "C8", f"A{res:04d}-{resi}:C8", "C", "",
         -1.404, -4.229, -0.339, start_idx + 12],
        ["ATOM", start_idx + 13, "A", res, "", resi, "N7", f"A{res:04d}-{resi}:N7", "N", "",
         -2.071, -3.103, -0.410, start_idx + 13],
        ["ATOM", start_idx + 14, "A", res, "", resi, "C5", f"A{res:04d}-{resi}:C5", "C", "",
         -1.104, -2.122, -0.233, start_idx + 14],
        ["ATOM", start_idx + 15, "A", res, "", resi, "C6", f"A{res:04d}-{resi}:C6", "C", "",
         -1.200, -0.709, -0.197, start_idx + 15],
        ["ATOM", start_idx + 16, "A", res, "", resi, "O6", f"A{res:04d}-{resi}:O6", "O", "",
         -2.218, -0.005, -0.326, start_idx + 16],
        ["ATOM", start_idx + 17, "A", res, "", resi, "N1", f"A{res:04d}-{resi}:N1", "N", "",
         0.000, -0.064, 0.000, start_idx + 17],
        ["ATOM", start_idx + 18, "A", res, "", resi, "C2", f"A{res:04d}-{resi}:C2", "C", "",
         1.168, -0.757, 0.154, start_idx + 18],
        ["ATOM", start_idx + 19, "A", res, "", resi, "N2", f"A{res:04d}-{resi}:N2", "N", "",
         2.252, 0.005, 0.331, start_idx + 19],
        ["ATOM", start_idx + 20, "A", res, "", resi, "N3", f"A{res:04d}-{resi}:N3", "N", "",
         1.320, -2.081, 0.139, start_idx + 20],
        ["ATOM", start_idx + 21, "A", res, "", resi, "C4", f"A{res:04d}-{resi}:C4", "C", "",
         0.130, -2.698, -0.060, start_idx + 21],
        ["ATOM", start_idx + 22, "A", res, "", resi, "H5'", f"A{res:04d}-{resi}:H5'", "H", "",
         0.894, -8.617, 2.072, start_idx + 22],
        ["ATOM", start_idx + 23, "A", res, "", resi, "H5''", f"A{res:04d}-{resi}:H5''", "H", "",
         0.522, -9.181, 0.667, start_idx + 23],
        ["ATOM", start_idx + 24, "A", res, "", resi, "H4'", f"A{res:04d}-{resi}:H4'", "H", "",
         2.319, -7.445, 0.742, start_idx + 24],
        ["ATOM", start_idx + 25, "A", res, "", resi, "H3'", f"A{res:04d}-{resi}:H3'", "H", "",
         0.382, -7.811, -1.213, start_idx + 25],
        ["ATOM", start_idx + 26, "A", res, "", resi, "H2'", f"A{res:04d}-{resi}:H2'", "H", "",
         -0.050, -5.769, -1.539, start_idx + 26],
        ["ATOM", start_idx + 27, "A", res, "", resi, "HO3'", f"A{res:04d}-{resi}:HO3'", "H", "",
         2.248, -7.567, -2.280, start_idx + 27],
        ["ATOM", start_idx + 28, "A", res, "", resi, "H1'", f"A{res:04d}-{resi}:H1'", "H", "",
         1.805, -4.584, 0.099, start_idx + 28],
        ["ATOM", start_idx + 29, "A", res, "", resi, "H8", f"A{res:04d}-{resi}:H8", "H", "",
         -1.804, -5.064, -0.430, start_idx + 29],
        ["ATOM", start_idx + 30, "A", res, "", resi, "H1", f"A{res:04d}-{resi}:H1", "H", "",
         0.015, 0.796, 0.026, start_idx + 30],
        ["ATOM", start_idx + 31, "A", res, "", resi, "H21", f"A{res:04d}-{resi}:H21", "H", "",
         3.023, -0.361, 0.435, start_idx + 31],
        ["ATOM", start_idx + 32, "A", res, "", resi, "H22", f"A{res:04d}-{resi}:H22", "H", "",
         2.178, 0.862, 0.341, start_idx + 32],
    ] + ([
        ["ATOM", start_idx + 33, "A", res, "", resi, "H2''", f"A{res:04d}-{resi}:H2''", "H", "",
         1.459, -5.550, -1.869, start_idx + 33],
    ] if is_dna else [
        ["ATOM", start_idx + 33, "A", res, "", resi, "O2'", f"A{res:04d}-{resi}:O2'", "O", "",
         0, 0, 0, start_idx + 33],
        ["ATOM", start_idx + 34, "A", res, "", resi, "HO2'", f"A{res:04d}-{resi}:HO2'", "H", "",
         0, 0, 0, start_idx + 34],
    ]), rotations), x_offset, y_offset, z_offset)
    hydrogen_mapping = [
        (f"A{res:04d}-{resi}:H5'", f"A{res:04d}-{resi}:C5'", 0.970),
        (f"A{res:04d}-{resi}:H5''", f"A{res:04d}-{resi}:C5'", 0.970),
        (f"A{res:04d}-{resi}:H4'", f"A{res:04d}-{resi}:C4'", 0.970),
        (f"A{res:04d}-{resi}:H3'", f"A{res:04d}-{resi}:C3'", 0.971),
        (f"A{res:04d}-{resi}:H2'", f"A{res:04d}-{resi}:C2'", 0.970),
        (f"A{res:04d}-{resi}:H1'", f"A{res:04d}-{resi}:C1'", 0.970),
        (f"A{res:04d}-{resi}:H8", f"A{res:04d}-{resi}:C8", 0.930),
        (f"A{res:04d}-{resi}:H21", f"A{res:04d}-{resi}:N2", 0.860),
        (f"A{res:04d}-{resi}:H22", f"A{res:04d}-{resi}:N2", 0.860),
        (f"A{res:04d}-{resi}:H1", f"A{res:04d}-{resi}:N1", 0.860),
    ] + ([
        (f"A{res:04d}-{resi}:H2''", f"A{res:04d}-{resi}:C2'", 0.970),
    ] if is_dna else [
        (f"A{res:04d}-{resi}:HO2'", f"A{res:04d}-{resi}:O2'", 1.04),
    ])
    covalent_edges = [
        (f"A{res:04d}-{resi}:OP1", f"A{res:04d}-{resi}:P",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.503, "base": True}),
        (f"A{res:04d}-{resi}:OP2", f"A{res:04d}-{resi}:P",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.488, "base": True}),
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:O5'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.632, "base": True}),
        (f"A{res:04d}-{resi}:O5'", f"A{res:04d}-{resi}:C5'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.483, "base": True}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:H5'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:H5''",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:C4'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.525, "base": True}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:H4'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C4'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.425, "base": True}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:O3'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 2.217, "base": True}),
        (f"A{res:04d}-{resi}:O3'", f"A{res:04d}-{resi}:C3'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.354, "base": True}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:H3'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.971, "base": True}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:C2'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.453, "base": True}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:H2'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:C1'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.396, "base": True}),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C1'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.406, "base": True}),
        (f"A{res:04d}-{resi}:C1'", f"A{res:04d}-{resi}:H1'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
        (f"A{res:04d}-{resi}:N1", f"A{res:04d}-{resi}:C2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.367, "base": True}),
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:N3",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.333, "base": True}),
        (f"A{res:04d}-{resi}:N2", f"A{res:04d}-{resi}:C2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.337, "base": True}),
        (f"A{res:04d}-{resi}:N3", f"A{res:04d}-{resi}:C4",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.355, "base": True}),
        (f"A{res:04d}-{resi}:O6", f"A{res:04d}-{resi}:C6",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.244, "base": True}),
        (f"A{res:04d}-{resi}:C6", f"A{res:04d}-{resi}:C5",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.417, "base": True}),
        (f"A{res:04d}-{resi}:C5", f"A{res:04d}-{resi}:C4",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.373, "base": True}),
        (f"A{res:04d}-{resi}:C4", f"A{res:04d}-{resi}:N9",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.379, "base": True}),
        (f"A{res:04d}-{resi}:N7", f"A{res:04d}-{resi}:C8",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.311, "base": True}),
        (f"A{res:04d}-{resi}:C8", f"A{res:04d}-{resi}:N9",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.362, "base": True}),
        (f"A{res:04d}-{resi}:N9", f"A{res:04d}-{resi}:C1'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.471, "base": True}),
        (f"A{res:04d}-{resi}:C8", f"A{res:04d}-{resi}:H8",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.930, "base": True}),
        (f"A{res:04d}-{resi}:N2", f"A{res:04d}-{resi}:H21",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.860, "base": True}),
        (f"A{res:04d}-{resi}:N2", f"A{res:04d}-{resi}:H22",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.860, "base": True}),
        (f"A{res:04d}-{resi}:N2", f"A{res:04d}-{resi}:C2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.337, "base": True}),
        (f"A{res:04d}-{resi}:N1", f"A{res:04d}-{resi}:H1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.860, "base": True}),
        (f"A{res:04d}-{resi}:C6", f"A{res:04d}-{resi}:N1",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.376, "base": True}),
        (f"A{res:04d}-{resi}:C5", f"A{res:04d}-{resi}:N7",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.389, "base": True}),
    ] + ([
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:H2''",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 0.970, "base": True}),
    ] if is_dna else [
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:O2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.44, "base": True}),
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:HO2'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.04, "base": True}),
    ])
    pebble_edges = [
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:OP1", {"weight": 0}),
        (f"A{res:04d}-{resi}:OP1", f"A{res:04d}-{resi}:P", {"weight": 5}),
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:OP2", {"weight": 0}),
        (f"A{res:04d}-{resi}:OP2", f"A{res:04d}-{resi}:P", {"weight": 5}),
        (f"A{res:04d}-{resi}:O5'", f"A{res:04d}-{resi}:P", {"weight": 0}),
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:O5'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:O5'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O5'", f"A{res:04d}-{resi}:C5'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:C5'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:C4'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:O4'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C4'", {"weight": 5}),
        (f"A{res:04d}-{resi}:O3'", f"A{res:04d}-{resi}:C4'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:O3'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:O3'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O3'", f"A{res:04d}-{resi}:C3'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:C3'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:C2'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C1'", f"A{res:04d}-{resi}:C2'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:C1'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:N1", {"weight": 0}),
        (f"A{res:04d}-{resi}:N1", f"A{res:04d}-{resi}:C2", {"weight": 5}),
        (f"A{res:04d}-{resi}:N3", f"A{res:04d}-{resi}:C2", {"weight": 0}),
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:N3", {"weight": 5}),
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:N2", {"weight": 0}),
        (f"A{res:04d}-{resi}:N2", f"A{res:04d}-{resi}:C2", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4", f"A{res:04d}-{resi}:N3", {"weight": 0}),
        (f"A{res:04d}-{resi}:N3", f"A{res:04d}-{resi}:C4", {"weight": 5}),
        (f"A{res:04d}-{resi}:C6", f"A{res:04d}-{resi}:O6", {"weight": 0}),
        (f"A{res:04d}-{resi}:O6", f"A{res:04d}-{resi}:C6", {"weight": 5}),
        (f"A{res:04d}-{resi}:C5", f"A{res:04d}-{resi}:C6", {"weight": 0}),
        (f"A{res:04d}-{resi}:C6", f"A{res:04d}-{resi}:C5", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4", f"A{res:04d}-{resi}:C5", {"weight": 0}),
        (f"A{res:04d}-{resi}:C5", f"A{res:04d}-{resi}:C4", {"weight": 5}),
        (f"A{res:04d}-{resi}:N9", f"A{res:04d}-{resi}:C4", {"weight": 0}),
        (f"A{res:04d}-{resi}:C4", f"A{res:04d}-{resi}:N9", {"weight": 5}),
        (f"A{res:04d}-{resi}:C8", f"A{res:04d}-{resi}:N7", {"weight": 0}),
        (f"A{res:04d}-{resi}:N7", f"A{res:04d}-{resi}:C8", {"weight": 5}),
        (f"A{res:04d}-{resi}:N9", f"A{res:04d}-{resi}:C8", {"weight": 0}),
        (f"A{res:04d}-{resi}:C8", f"A{res:04d}-{resi}:N9", {"weight": 5}),
        (f"A{res:04d}-{resi}:C1'", f"A{res:04d}-{resi}:N9", {"weight": 0}),
        (f"A{res:04d}-{resi}:N9", f"A{res:04d}-{resi}:C1'", {"weight": 5}),
    ] + ([
        # not required
    ] if is_dna else [
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:O2'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:C2'", {"weight": 5}),
    ])
    pebble_map = {
        f"A{res:04d}-{resi}:N1": 1,
        # Backbone
        f"A{res:04d}-{resi}:OP1": 1,
        f"A{res:04d}-{resi}:OP2": 1,
        f"A{res:04d}-{resi}:P": 1,
        f"A{res:04d}-{resi}:O5'": 1,
        f"A{res:04d}-{resi}:C5'": 1,
        f"A{res:04d}-{resi}:O4'": 1,
        f"A{res:04d}-{resi}:C4'": 1,
        f"A{res:04d}-{resi}:C3'": 1,
        f"A{res:04d}-{resi}:O3'": 6,
        f"A{res:04d}-{resi}:C2'": 1,
        f"A{res:04d}-{resi}:C1'": 1,
    }
    if not is_dna:
        pebble_map = dict(**pebble_map, **{
            f"A{res:04d}-{resi}:O2'": 1,
        })
    pebble_bars = [
        # Backbone
        (f"A{res:04d}-{resi}:OP1", f"A{res:04d}-{resi}:P", 1),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C1'", 5),
    ] + ([
        # not required
    ] if is_dna else [
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:C2'", 1),
    ])

    return format_coords(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )

def graph_dna_thymine(
        start_idx: int = 1,
        res: int = 1,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        z_offset: float = 0.0,
        is_dna: bool = True,
        rotations: Optional[List[Tuple[float, AXIS_OPTIONS]]] = None,
        exclude: Optional[List[str]] = None
):
    resi = ("DT" if is_dna else "U").rjust(3, " ")
    coords = modify_coords(apply_rotations([
        # Nucleotide heavy atoms
        ["ATOM", start_idx, "A", res, "", resi, "N3", f"A{res:04d}-{resi}:N3", "N", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C2", f"A{res:04d}-{resi}:C2", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "O2", f"A{res:04d}-{resi}:O2", "O", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "N1", f"A{res:04d}-{resi}:N1", "N", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C4", f"A{res:04d}-{resi}:C4", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "O4", f"A{res:04d}-{resi}:O4", "O", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C5", f"A{res:04d}-{resi}:C5", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C6", f"A{res:04d}-{resi}:C6", "C", "",
         0, 0, 0, start_idx],
        # Nucleotide hydrogen atoms
        ["ATOM", start_idx, "A", res, "", resi, "H6", f"A{res:04d}-{resi}:H6", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H3", f"A{res:04d}-{resi}:H3", "H", "",
         0, 0, 0, start_idx],
        # Backbone
        ["ATOM", start_idx, "A", res, "", resi, "C1'", f"A{res:04d}-{resi}:C1'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C2'", f"A{res:04d}-{resi}:C2'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C3'", f"A{res:04d}-{resi}:C3'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "O3'", f"A{res:04d}-{resi}:O3'", "O", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C4'", f"A{res:04d}-{resi}:C4'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "O4'", f"A{res:04d}-{resi}:O4'", "O", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C5'", f"A{res:04d}-{resi}:C5'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "O5'", f"A{res:04d}-{resi}:O5'", "O", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C4'", f"A{res:04d}-{resi}:C2'", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H1'", f"A{res:04d}-{resi}:H1'", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H3'", f"A{res:04d}-{resi}:H3'", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H4'", f"A{res:04d}-{resi}:H4'", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H5'", f"A{res:04d}-{resi}:H5'", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H5''", f"A{res:04d}-{resi}:H5''", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "P", f"A{res:04d}-{resi}:P", "P", "",
         1.861, -8.579, 1.742, start_idx],
        ["ATOM", start_idx + 1, "A", res, "", resi, "OP1", f"A{res:04d}-{resi}:OP1", "O", "",
         1.588, -9.922, 2.372, start_idx + 1],
        ["ATOM", start_idx + 2, "A", res, "", resi, "OP2", f"A{res:04d}-{resi}:OP2", "O", "",
         3.084, -8.504, 0.917, start_idx + 2],
    ] + ([
        ["ATOM", start_idx, "A", res, "", resi, "H2''", f"A{res:04d}-{resi}:H2''", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "C7", f"A{res:04d}-{resi}:C7", "C", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H71", f"A{res:04d}-{resi}:H71", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H72", f"A{res:04d}-{resi}:H72", "H", "",
         0, 0, 0, start_idx],
        ["ATOM", start_idx, "A", res, "", resi, "H73", f"A{res:04d}-{resi}:H73", "H", "",
         0, 0, 0, start_idx],
    ] if is_dna else [
        ["ATOM", start_idx, "A", res, "", resi, "H5", f"A{res:04d}-{resi}:H5", "H", "",
         0, 0, 0, start_idx],
    ]), rotations), x_offset, y_offset, z_offset)
    hydrogen_mapping = [
        (f"A{res:04d}-{resi}:H6", f"A{res:04d}-{resi}:C6", 1.09),
        (f"A{res:04d}-{resi}:H3", f"A{res:04d}-{resi}:N3", 1.07),
        # Backbone
        (f"A{res:04d}-{resi}:H5'", f"A{res:04d}-{resi}:C5'", 1.14),
        (f"A{res:04d}-{resi}:H5''", f"A{res:04d}-{resi}:C5'", 1.14),
        (f"A{res:04d}-{resi}:H4'", f"A{res:04d}-{resi}:C4'", 1.14),
        (f"A{res:04d}-{resi}:H3'", f"A{res:04d}-{resi}:C3'", 1.14),
        (f"A{res:04d}-{resi}:H2'", f"A{res:04d}-{resi}:C2'", 1.14),
        (f"A{res:04d}-{resi}:H1'", f"A{res:04d}-{resi}:C1'", 1.14),
    ] + ([
        (f"A{res:04d}-{resi}:H2''", f"A{res:04d}-{resi}:C2'", 1.14),
        (f"A{res:04d}-{resi}:H71", f"A{res:04d}-{resi}:C7'", 1.14),
        (f"A{res:04d}-{resi}:H72", f"A{res:04d}-{resi}:C7'", 1.14),
        (f"A{res:04d}-{resi}:H73", f"A{res:04d}-{resi}:C7'", 1.14),
    ] if is_dna else [
        (f"A{res:04d}-{resi}:H5", f"A{res:04d}-{resi}:C5'", 1.09),
        (f"A{res:04d}-{resi}:HO2'", f"A{res:04d}-{resi}:O2'", 1.04),
    ])
    covalent_edges = [
        (f"A{res:04d}-{resi}:N1", f"A{res:04d}-{resi}:C2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        # Backbone
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:H5''",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:H5'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:H5''",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:H4'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:H3'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:H2'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C1'", f"A{res:04d}-{resi}:H1'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C8", f"A{res:04d}-{resi}:H8",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.09, "base": True}),
        (f"A{res:04d}-{resi}:N6", f"A{res:04d}-{resi}:H61",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-{resi}:N6", f"A{res:04d}-{resi}:H62",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.07, "base": True}),
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:H2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.04, "base": True}),
    ] + ([
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:H2''",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C7", f"A{res:04d}-{resi}:C5",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.44, "base": True}),
        (f"A{res:04d}-{resi}:C7", f"A{res:04d}-{resi}:H71",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C7", f"A{res:04d}-{resi}:H72",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
        (f"A{res:04d}-{resi}:C7", f"A{res:04d}-{resi}:H73",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.14, "base": True}),
    ] if is_dna else [
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:O2",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.44, "base": True}),
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:HO2'",
         {"kind": {InteractionType.COVALENT.value}, "bond_length": 1.04, "base": True}),
    ])
    pebble_edges = [
        (f"A{res:04d}-{resi}:C2", f"A{res:04d}-{resi}:N1", {"weight": 0}),
        (f"A{res:04d}-{resi}:N1", f"A{res:04d}-{resi}:C2", {"weight": 5}),
        # Backbone
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:OP1", {"weight": 0}),
        (f"A{res:04d}-{resi}:OP1", f"A{res:04d}-{resi}:P", {"weight": 5}),
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:OP2", {"weight": 0}),
        (f"A{res:04d}-{resi}:OP2", f"A{res:04d}-{resi}:P", {"weight": 5}),
        (f"A{res:04d}-{resi}:O5'", f"A{res:04d}-{resi}:P", {"weight": 0}),
        (f"A{res:04d}-{resi}:P", f"A{res:04d}-{resi}:O5'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:O5'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O5'", f"A{res:04d}-{resi}:C5'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:C5'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C5'", f"A{res:04d}-{resi}:C4'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:O4'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C4'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:C4'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C4'", f"A{res:04d}-{resi}:C3'", {"weight": 5}),
        (f"A{res:04d}-{resi}:O3'", f"A{res:04d}-{resi}:C3'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:O3'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C3'", f"A{res:04d}-{resi}:C2'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:C3'", {"weight": 5}),
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:C1'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C1'", f"A{res:04d}-{resi}:C2'", {"weight": 5}),
    ] + ([
        (f"A{res:04d}-{resi}:C5", f"A{res:04d}-{resi}:C7'", {"weight": 0}),
        (f"A{res:04d}-{resi}:C7", f"A{res:04d}-{resi}:C5", {"weight": 5}),
    ] if is_dna else [
        (f"A{res:04d}-{resi}:C2'", f"A{res:04d}-{resi}:O2'", {"weight": 0}),
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:C2'", {"weight": 5}),
    ])
    pebble_map = dict(**{
        f"A{res:04d}-{resi}:N1": 1,
        # Backbone
        f"A{res:04d}-{resi}:OP1": 1,
        f"A{res:04d}-{resi}:OP2": 1,
        f"A{res:04d}-{resi}:P": 1,
        f"A{res:04d}-{resi}:O5'": 1,
        f"A{res:04d}-{resi}:C5'": 1,
        f"A{res:04d}-{resi}:O4'": 1,
        f"A{res:04d}-{resi}:C4'": 1,
        f"A{res:04d}-{resi}:C3'": 1,
        f"A{res:04d}-{resi}:O3'": 6,
        f"A{res:04d}-{resi}:C2'": 1,
        f"A{res:04d}-{resi}:C1'": 1,
    }, **({
        f"A{res:04d}-{resi}:C7": 1,
    } if is_dna else {
        f"A{res:04d}-{resi}:O2'": 1,
    }))
    pebble_bars = [
        # Backbone
        (f"A{res:04d}-{resi}:OP1", f"A{res:04d}-{resi}:P", 1),
        (f"A{res:04d}-{resi}:O4'", f"A{res:04d}-{resi}:C1'", 5),
    ] + ([
        # not required
    ] if is_dna else [
        (f"A{res:04d}-{resi}:O2'", f"A{res:04d}-{resi}:C2'", 1),
    ])

    return format_coords(
        coords=coords, hydrogen_mapping=hydrogen_mapping, covalent_edges=covalent_edges,
        pebble_edges=pebble_edges, pebble_map=pebble_map, pebble_bars=pebble_bars,
        exclude=exclude
    )
