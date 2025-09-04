from typing import Dict, List, Tuple, Union

##############################
# Pebble Graph Constants #####
##############################

AA_RESI_GRAPH_TEMPLATES: Dict[str, Dict[str, List[Union[Tuple[str, int], Tuple[str, int, bool]]]]] = {
    'ALA': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)]
    },
    'ARG': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD': [('CG', 5)],
        'NE': [('CD', 5)],
        'CZ': [('NE', 5)],
        'NH1': [('CZ', 5)],
        'NH2': [('CZ', 5), ('NH1', 2, False)]
    },
    'ASN': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'OD1': [('CG', 6)],
        'ND2': [('CG', 5)]
    },
    'ASP': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'OD1': [('CG', 5)],
        'OD2': [('CG', 5), ('OD1', 2, False)]
    },
    'CYH': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'SG': [('CB', 5)]
    },
    'CYS': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'SG': [('CB', 5)]
    },
    'CSS': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'SG': [('CB', 5)]
    },
    'GLU': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD': [('CG', 5)],
        'OE1': [('CD', 5)],
        'OE2': [('CD', 5), ('OE1', 2, False)]
    },
    'GLN': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD': [('CG', 5)],
        'OE1': [('CD', 6)],
        'NE2': [('CD', 5)]
    },
    'GLY': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)]
    },
    'HIS': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD2': [('CG', 5)],
        'ND1': [('CG', 5)],
        'CE1': [('ND1', 5), ('NE2', 5)],
        'NE2': [('CD2', 5)],
    },
    'HID': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD2': [('CG', 5)],
        'ND1': [('CG', 5)],
        'CE1': [('ND1', 5), ('NE2', 5)],
        'NE2': [('CD2', 5)],
    },
    'HIE': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD2': [('CG', 5)],
        'ND1': [('CG', 5)],
        'CE1': [('ND1', 5), ('NE2', 5)],
        'NE2': [('CD2', 5)],
    },
    'HIP': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD2': [('CG', 5)],
        'ND1': [('CG', 5)],
        'CE1': [('ND1', 5), ('NE2', 5)],
        'NE2': [('CD2', 5)],
    },
    'ILE': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG1': [('CB', 5)],
        'CG2': [('CB', 5)],
        'CD1': [('CG1', 5)]
    },
    'LEU': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD1': [('CG', 5)],
        'CD2': [('CG', 5)]
    },
    'LYS': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD': [('CG', 5)],
        'CE': [('CD', 5)],
        'NZ': [('CE', 5)]
    },
    'MET': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'SD': [('CG', 5)],
        'CE': [('SD', 5)]
    },
    'PHE': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD1': [('CG', 5)],
        'CD2': [('CG', 5)],
        'CE1': [('CD1', 5)],
        'CE2': [('CD2', 5)],
        'CZ': [('CE1', 5), ('CE2', 5)]
    },
    'PRO': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD': [('N', 5), ('CG', 5)]
    },
    'SER': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'OG': [('CB', 5)]
    },
    'THR': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'OG1': [('CB', 5)],
        'CG2': [('CB', 5)]
    },
    'TRP': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD1': [('CG', 5)],
        'CD2': [('CG', 5)],
        'NE1': [('CD1', 5), ('CE2', 5)],
        'CE2': [('CD2', 5)],
        'CE3': [('CD2', 5)],
        'CZ2': [('CE2', 5)],
        'CZ3': [('CE3', 5)],
        'CH2': [('CZ2', 5), ('CZ3', 5)]
    },
    'TYR': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG': [('CB', 5)],
        'CD1': [('CG', 5)],
        'CD2': [('CG', 5)],
        'CE1': [('CD1', 5)],
        'CE2': [('CD2', 5)],
        'CZ': [('CE1', 5), ('CE2', 5)],
        'OH': [('CZ', 5)]
    },
    'VAL': {
        'N': [('CA', 5)],
        'CA': [('C', 5)],
        'O': [('C', 6)],
        'CB': [('CA', 5)],
        'CG1': [('CB', 5)],
        'CG2': [('CB', 5)]
    }
}
"""
Templates for mapping well-defined pebble graph edges for each amino acid.

Third value in tuple indicates whether the edge should be counted as regular covalent bond.
If ``False``, the edge only counts towards pebble graph stability. Default ``True``.
"""

NA_RESI_GRAPH_TEMPLATES: Dict[str, Dict[str, List[Union[Tuple[str, int], Tuple[str, int, bool]]]]] = {
    '  A': {
        'OP1': [('P', 5), ('OP2', 2, False)],
        'OP2': [('P', 5)],
        'P': [("O5'", 5)],
        "O5'": [("C5'", 5)],
        "C5'": [("C4'", 5)],
        "O4'": [("C4'", 5), ("C1'", 5)],
        "C4'": [("C3'", 5)],
        "C3'": [("O3'", 5)],
        "O2'": [("C2'", 5)],
        "C2'": [("C3'", 5)],
        "C1'": [("C2'", 5)],
        'N9': [("C1'", 5)],
        'C8': [('N9', 5)],
        'N7': [('C8', 6)],
        'C4': [('N9', 6)],
        'C5': [('C4', 6), ('N7', 5)],
        'C6': [('C5', 5), ('N1', 6)],
        'N6': [('C6', 5)],
        'N3': [('C4', 5)],
        'C2': [('N3', 6)],
        'N1': [('C2', 5)],
    },
    '  C': {
        'OP1': [('P', 5), ('OP2', 2, False)],
        'OP2': [('P', 5)],
        'P': [("O5'", 5)],
        "O5'": [("C5'", 5)],
        "C5'": [("C4'", 5)],
        "O4'": [("C4'", 5), ("C1'", 5)],
        "C4'": [("C3'", 5)],
        "C3'": [("O3'", 5)],
        "O2'": [("C2'", 5)],
        "C2'": [("C3'", 5)],
        "C1'": [("C2'", 5)],
        'N1': [("C1'", 5)],
        'C6': [('N1', 5)],
        'C5': [('C6', 6)],
        'C2': [('N1', 5)],
        'O2': [('C2', 6)],
        'N3': [('C2', 5)],
        'C4': [('N3', 6), ('C5', 5)],
        'N4': [('C4', 5)],
    },
    '  G': {
        'OP1': [('P', 5), ('OP2', 2, False)],
        'OP2': [('P', 5)],
        'P': [("O5'", 5)],
        "O5'": [("C5'", 5)],
        "C5'": [("C4'", 5)],
        "O4'": [("C4'", 5), ("C1'", 5)],
        "C4'": [("C3'", 5)],
        "C3'": [("O3'", 5)],
        "O2'": [("C2'", 5)],
        "C2'": [("C3'", 5)],
        "C1'": [("C2'", 5)],
        'N9': [("C1'", 5)],
        'C8': [('N9', 5)],
        'N7': [('C8', 6)],
        'C4': [('N9', 6)],
        'C5': [('C4', 6), ('N7', 5)],
        'C6': [('C5', 5), ('N1', 5)],
        'O6': [('C6', 6)],
        'N3': [('C4', 5)],
        'C2': [('N3', 6)],
        'N1': [('C2', 5)],
        'N2': [('C2', 5)],
    },
    '  U': {
        'OP1': [('P', 5), ('OP2', 2, False)],
        'OP2': [('P', 5)],
        'P': [("O5'", 5)],
        "O5'": [("C5'", 5)],
        "C5'": [("C4'", 5)],
        "O4'": [("C4'", 5), ("C1'", 5)],
        "C4'": [("C3'", 5)],
        "C3'": [("O3'", 5)],
        "O2'": [("C2'", 5)],
        "C2'": [("C3'", 5)],
        "C1'": [("C2'", 5)],
        'N1': [("C1'", 5)],
        'C6': [('N1', 5)],
        'C5': [('C6', 6)],
        'C2': [('N1', 5)],
        'O2': [('C2', 6)],
        'N3': [('C2', 5)],
        'C4': [('N3', 5), ('C5', 5)],
        'O4': [('C4', 6)],
    },
    ' DA': {
        'OP1': [('P', 5), ('OP2', 2, False)],
        'OP2': [('P', 5)],
        'P': [("O5'", 5)],
        "O5'": [("C5'", 5)],
        "C5'": [("C4'", 5)],
        "O4'": [("C4'", 5), ("C1'", 5)],
        "C4'": [("C3'", 5)],
        "C3'": [("O3'", 5)],
        "C2'": [("C3'", 5)],
        "C1'": [("C2'", 5)],
        'N9': [("C1'", 5)],
        'C8': [('N9', 5)],
        'N7': [('C8', 6)],
        'C4': [('N9', 6)],
        'C5': [('C4', 6), ('N7', 5)],
        'C6': [('C5', 5), ('N1', 6)],
        'N6': [('C6', 5)],
        'N3': [('C4', 5)],
        'C2': [('N3', 6)],
        'N1': [('C2', 5)],
    },
    ' DC': {
        'OP1': [('P', 5), ('OP2', 2, False)],
        'OP2': [('P', 5)],
        'P': [("O5'", 5)],
        "O5'": [("C5'", 5)],
        "C5'": [("C4'", 5)],
        "O4'": [("C4'", 5), ("C1'", 5)],
        "C4'": [("C3'", 5)],
        "C3'": [("O3'", 5)],
        "C2'": [("C3'", 5)],
        "C1'": [("C2'", 5)],
        'N1': [("C1'", 5)],
        'C6': [('N1', 5)],
        'C5': [('C6', 6)],
        'C2': [('N1', 5)],
        'O2': [('C2', 6)],
        'N3': [('C2', 5)],
        'C4': [('N3', 6), ('C5', 5)],
        'N4': [('C4', 5)],
    },
    ' DG': {
        'OP1': [('P', 5), ('OP2', 2, False)],
        'OP2': [('P', 5)],
        'P': [("O5'", 5)],
        "O5'": [("C5'", 5)],
        "C5'": [("C4'", 5)],
        "O4'": [("C4'", 5), ("C1'", 5)],
        "C4'": [("C3'", 5)],
        "C3'": [("O3'", 5)],
        "C2'": [("C3'", 5)],
        "C1'": [("C2'", 5)],
        'N9': [("C1'", 5)],
        'C8': [('N9', 5)],
        'N7': [('C8', 6)],
        'C4': [('N9', 6)],
        'C5': [('C4', 6), ('N7', 5)],
        'C6': [('C5', 5), ('N1', 5)],
        'O6': [('C6', 6)],
        'N3': [('C4', 5)],
        'C2': [('N3', 6)],
        'N1': [('C2', 5)],
        'N2': [('C2', 5)],
    },
    ' DT': {
        'OP1': [('P', 5), ('OP2', 2, False)],
        'OP2': [('P', 5)],
        'P': [("O5'", 5)],
        "O5'": [("C5'", 5)],
        "C5'": [("C4'", 5)],
        "O4'": [("C4'", 5), ("C1'", 5)],
        "C4'": [("C3'", 5)],
        "C3'": [("O3'", 5)],
        "C2'": [("C3'", 5)],
        "C1'": [("C2'", 5)],
        'N1': [("C1'", 5)],
        'C6': [('N1', 5)],
        'C5': [('C6', 6)],
        'C7': [('C5', 5)],
        'C2': [('N1', 5)],
        'O2': [('C2', 6)],
        'N3': [('C2', 5)],
        'C4': [('N3', 5), ('C5', 5)],
        'O4': [('C4', 6)],
    }
}
"""
Templates for mapping well-defined pebble graph edges for nucleotide residue.

Third value in tuple indicates whether the edge should be counted as regular covalent bond.
If ``False``, the edge only counts towards pebble graph stability. Default ``True``.
"""

RESI_GRAPH_TEMPLATES: Dict[str, Dict[str, List[Union[Tuple[str, int], Tuple[str, int, bool]]]]] = \
    dict(**AA_RESI_GRAPH_TEMPLATES, **NA_RESI_GRAPH_TEMPLATES)

AA_RESIDUES = list(AA_RESI_GRAPH_TEMPLATES.keys())

NA_RESIDUES = list(NA_RESI_GRAPH_TEMPLATES.keys())

STANDARD_RESIDUES = AA_RESIDUES + NA_RESIDUES

HETATM_GRAPH_TEMPLATES: Dict[str, Dict[str, List[Tuple[str, int]]]] = {
    'ACE': {
        'O': [('C', 6)],
        'CH3': [('C', 5)]
    }
}

##########################
# Covalent Constants #####
##########################

MAX_NUMBER_OF_BONDS: Dict[str, int] = {
    'C': 4,
    'O': 2,
    'N': 3,
    'H': 1,
    'S': 2
}

########################
# H-Bond Constants #####
########################

RESIDUE_HBOND_DONORS: Dict[str, Dict[str, int]] = {
    "ARG": {"NE": 1, "NH1": 2, "NH2": 2},
    "ASN": {"ND2": 2},
    "CYS": {"SG": 1},
    "CYH": {"SG": 1},
    "GLN": {"NE2": 2},
    "HIS": {"ND1": 1, "NE2": 1},
    "HID": {"ND1": 1, "NE2": 0},
    "HIE": {"ND1": 0, "NE2": 1},
    "HIP": {"ND1": 1, "NE2": 1},
    "LYS": {"NZ": 3},
    "SER": {"OG": 1},
    "THR": {"OG1": 1},
    "TRP": {"NE1": 1},
    "TYR": {"OH": 1},
    "PRO": {"N": 0},  # explicitly exclude Proline main chain N
    "  A": {"N6": 2, "C2": 1},
    "  C": {"N4": 2},
    "  G": {"N1": 1, "N2": 2},
    "  U": {"N3": 1},
    " DA": {"N6": 2, "C2": 1},
    " DC": {"N4": 2},
    " DG": {"N1": 1, "N2": 2},
    " DT": {"N3": 1},
}
"""
Number of hydrogen bonds each donor atom can enter.

https://www.ccpn.ac.uk/manual/v3/NEFAtomNames.html

https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/manual.html

https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/
"""

RESIDUE_HBOND_ACCEPTORS: Dict[str, Dict[str, int]] = {
    "ASN": {"OD1": 2},
    "ASP": {"OD1": 2, "OD2": 2},
    "CYS": {"SG": 1},
    "CSS": {"SG": 1},
    "CYH": {"SG": 1},
    "GLN": {"OE1": 2},
    "GLU": {"OE1": 2, "OE2": 2},
    "HIS": {"ND1": 1, "NE2": 1},  # see HBPLUS manual
    "HID": {"ND1": 0, "NE2": 1},
    "HIE": {"ND1": 1, "NE2": 0},
    "MET": {"SD": 1},
    "SER": {"OG": 2},
    "THR": {"OG1": 2},
    "TYR": {"OH": 1},
    "  A": {"N1": 1},
    "  C": {"O2": 2, "N3": 1},
    "  G": {"O6": 2},
    "  U": {"O2": 2, "O4": 2},
    " DA": {"N1": 1},
    " DC": {"O2": 2, "N3": 1},
    " DG": {"O6": 2},
    " DT": {"O2": 2, "O4": 2},
}
"""
Number of hydrogen bonds each acceptor atom can enter.

https://www.ccpn.ac.uk/manual/v3/NEFAtomNames.html

https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/manual.html

https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/"""

DEFAULT_HYBRIDIZATION: Dict[str, int] = {
    "N": 2,
    "CA": 3,
    "C": 2,
    "O": 2,
    "OXT": 2,
    "CB": 3,
    "P": 3,
    "OP1": 2,
    "OP2": 2,
    "O5'": 3,
    "C5'": 3,
    "O4'": 3,
    "C4'": 3,
    "O3'": 3,
    "C3'": 3,
    "O2'": 3,
    "C1'": 3,
}
"""


Hybridization states for the DNA and RNA phosphate is taken from
`Müller (2017)`_.

.. _Müller (2017):
   https://doi.org/10.1007/s40828-017-0046-8
"""

RESIDUE_ATOM_HYBRIDIZATION: Dict[str, Dict[str, int]] = {
    "VAL": {"CG1": 3, "CG2": 3},
    "LEU": {"CG": 3, "CD1": 3, "CD2": 3},
    "ILE": {"CG1": 3, "CG2": 3, "CD1": 3},
    "MET": {"CG": 3, "SD": 3, "CE": 3},
    "PHE": {"CG": 2, "CD1": 2, "CD2": 2, "CE1": 2, "CE2": 2, "CZ": 2},
    "PRO": {"CG": 3, "CD": 3},
    "SER": {"OG": 3},
    "THR": {"OG1": 3, "CG2": 3},
    "CYS": {"SG": 3},
    "CSS": {"SG": 3},
    "CYH": {"SG": 3},
    "ASN": {"CG": 2, "OD1": 2, "ND2": 2},
    "GLN": {"CG": 3, "CD": 2, "OE1": 2, "NE2": 2},
    "TYR": {"CG": 2, "CD1": 2, "CD2": 2, "CE1": 2, "CE2": 2, "CZ": 2, "OH": 2},
    "TRP": {"CG": 2, "CD1": 2, "CD2": 2, "NE1": 2, "CE2": 2, "CE3": 2, "CZ2": 2, "CZ3": 2, "CH2": 2},
    "ASP": {"CG": 2, "OD1": 2, "OD2": 2},
    "GLU": {"CG": 3, "CD": 2, "OE1": 2, "OE2": 2},
    "HIS": {"CG": 2, "CD2": 2, "ND1": 2, "CE1": 2, "NE2": 2},
    "HID": {"CG": 2, "CD2": 2, "ND1": 3, "CE1": 2, "NE2": 2},
    "HIE": {"CG": 2, "CD2": 2, "ND1": 2, "CE1": 2, "NE2": 2},
    "HIP": {"CG": 2, "CD2": 2, "ND1": 2, "CE1": 2, "NE2": 2},
    "LYS": {"CG": 3, "CD": 3, "CE": 3, "NZ": 3},
    "ARG": {"CG": 3, "CD": 3, "NE": 2, "CZ": 2, "NH1": 2, "NH2": 2},
    "  A": {"N9": 2, "C8": 2, "N7": 2, "C6": 2, "N6": 3, "C5": 2, "C4": 2, "N3": 2, "C2": 2, "N1": 2},
    "  C": {"N1": 2, "C2": 2, "O2": 2, "N3": 2, "C4": 3, "N4": 3, "C5": 2, "C6": 2},
    "  G": {"N9": 2, "C8": 2, "N7": 2, "C6": 2, "O6": 2, "C5": 2, "C4": 2, "N3": 2, "C2": 2, "N1": 2, "N2": 3},
    "  U": {"N1": 2, "C2": 2, "O2": 2, "N3": 2, "C4": 3, "O4": 2, "C5": 2, "C6": 2},
    " DA": {"N9": 2, "C8": 2, "N7": 2, "C6": 2, "N6": 3, "C5": 2, "C4": 2, "N3": 2, "C2": 2, "N1": 2},
    " DC": {"N1": 2, "C2": 2, "O2": 2, "N3": 2, "C4": 3, "N4": 3, "C5": 2, "C6": 2},
    " DG": {"N9": 2, "C8": 2, "N7": 2, "C6": 2, "O6": 2, "C5": 2, "C4": 2, "N3": 2, "C2": 2, "N1": 2, "N2": 3},
    " DT": {"N1": 2, "C2": 2, "O2": 2, "N3": 2, "C4": 3, "O4": 2, "C5": 2, "C6": 2, "C7": 3},
}

#################################
# Disulphide Bond Constants #####
#################################

DISULPHIDE_RESIS: List[str] = ["CYS", "CSS"]

DISULPHIDE_ATOMS: List[str] = ["SG"]

#################################
# (Aromatic) Ring Constants #####
#################################

RING_NORMAL_ATOMS: Dict[str, Dict[str, int]] = {
    "PHE": {"CG": 0, "CE1": 1, "CE2": 2},
    "TRP": {"CD2": 0, "CZ2": 1, "CZ3": 2},
    "TYR": {"CG": 0, "CE1": 1, "CE2": 2},
    "  A": {"N9": 0, "C6": 1, "C2": 2},
    "  C": {"N1": 0, "N3": 1, "C5": 2},
    "  G": {"N9": 0, "C6": 1, "C2": 2},
    "  U": {"N1": 0, "N3": 1, "C5": 2},
    " DA": {"N9": 0, "C6": 1, "C2": 2},
    " DC": {"N1": 0, "N3": 1, "C5": 2},
    " DG": {"N9": 0, "C6": 1, "C2": 2},
    " DT": {"N1": 0, "N3": 1, "C5": 2},
}

RIGID_RING_ATOMS: Dict[str, List[str]] = {
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HID": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HIE": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HIP": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["N", "CA", "CB", "CG", "CD"],
    "  A": ["N9", "C8", "N7", "C6", "C5", "C4", "N3", "C2", "N1"],
    "  C": ["N1", "C2", "O2", "N3", "C4", "C5", "C6"],
    "  G": ["N9", "C8", "N7", "C6", "O6", "C5", "C4", "N3", "C2", "N1"],
    "  U": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"],
    " DA": ["N9", "C8", "N7", "C6", "C5", "C4", "N3", "C2", "N1"],
    " DC": ["N1", "C2", "O2", "N3", "C4", "C5", "C6"],
    " DG": ["N9", "C8", "N7", "C6", "O6", "C5", "C4", "N3", "C2", "N1"],
    " DT": ["N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"],
}
