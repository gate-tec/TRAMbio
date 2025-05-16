from typing import List

from TRAMbio.util.structure_library.pdb_entry_map import EntryMap


#####################
# PDB Constants #####
#####################

ATOM_DF_COLUMNS = [
    'record_name', 'atom_number', 'atom_name', 'residue_name', 'chain_id', 'residue_number', 'insertion', 'element_symbol',
    'charge', 'x_coord', 'y_coord', 'z_coord',
    'line_idx'
]

HEADER_RECORDS = [
    'HEADER', 'TITLE',
    'COMPND', 'SOURCE', 'KEYWDS', 'EXPDATA',
    'AUTHOR', 'REVDAT',
    'REMARK',
    'SEQRES',
    'HET', 'HETNAM', 'HETSYN', 'FORMUL',  # HETATM meta info
    'HELIX', 'SHEET',
    'CRYST1', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3'
]

SS_BOND_MAP: List[EntryMap] = [
    {'id': 'ss_entry_number', 'line': [0, 4], "type": int, "required": False},
    {'id': 'residue_name_1', 'line': [5, 8], "type": str, "required": True},
    {'id': 'chain_id_1', 'line': [9, 10], "type": str, "required": True},
    {'id': 'residue_number_1', 'line': [11, 15], "type": int, "required": True},
    {'id': 'residue_name_2', 'line': [19, 22], "type": str, "required": True},
    {'id': 'chain_id_2', 'line': [23, 24], "type": str, "required": True},
    {'id': 'residue_number_2', 'line': [25, 29], "type": int, "required": True},
    {'id': 'bond_length', 'line': [66], "type": float, "required": False},
]

LINK_MAP: List[EntryMap] = [
    {'id': 'atom_name_1', 'line': [6, 10], "type": str, "required": True},
    {'id': 'residue_name_1', 'line': [11, 14], "type": str, "required": True},
    {'id': 'chain_id_1', 'line': [15, 16], "type": str, "required": True},
    {'id': 'residue_number_1', 'line': [16, 20], "type": int, "required": True},
    {'id': 'atom_name_2', 'line': [36, 40], "type": str, "required": True},
    {'id': 'residue_name_2', 'line': [41, 44], "type": str, "required": True},
    {'id': 'chain_id_2', 'line': [45, 46], "type": str, "required": True},
    {'id': 'residue_number_2', 'line': [46, 50], "type": int, "required": True},
    {'id': 'bond_length', 'line': [66], "type": float, "required": False},
]
