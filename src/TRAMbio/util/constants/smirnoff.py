from typing import TypedDict, List, Dict


#############
# Types #####
#############

class SmirnoffTemplate(TypedDict):
    epsilon: float
    rmin_half: float


class SmirnoffParameter(TypedDict):
    smirks: str
    epsilon: float
    rmin_half: float


#################
# Constants #####
#################


DEFAULT_POTENTIAL: Dict[str, SmirnoffTemplate] = {
    "N": {"epsilon": 0.17, "rmin_half": 1.824},
    "CA": {"epsilon": 0.1094, "rmin_half": 1.908},
    "C": {"epsilon": 0.086, "rmin_half": 1.908},
    "O": {"epsilon": 0.21, "rmin_half": 1.6612},
    "CB": {"epsilon": 0.1094, "rmin_half": 1.908},
    "H": {"epsilon": 0.0157, "rmin_half": 0.6},
    "HA": {"epsilon": 0.0157, "rmin_half": 1.387},
}


RESIDUE_ATOM_POTENTIAL: Dict[str, Dict[str, SmirnoffTemplate]] = {
    "VAL": {
        "CG1": {"epsilon": 0.1094, "rmin_half": 1.908},
        "CG2": {"epsilon": 0.1094, "rmin_half": 1.908},
        "HB": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG11": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG12": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG13": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG21": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG22": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG23": {"epsilon": 0.0157, "rmin_half": 1.487}
    },
    "LEU": {
        "CG": {"epsilon": 0.1094, "rmin_half": 1.908},
        "CD1": {"epsilon": 0.1094, "rmin_half": 1.908},
        "CD2": {"epsilon": 0.1094, "rmin_half": 1.908},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD11": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD12": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD13": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD21": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD22": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD23": {"epsilon": 0.0157, "rmin_half": 1.487}
    },
    "ILE": {
        "CG1": {"epsilon": 0.1094, "rmin_half": 1.908},
        "CG2": {"epsilon": 0.1094, "rmin_half": 1.908},
        "CD1": {"epsilon": 0.1094, "rmin_half": 1.908},
        "HB": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG12": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG13": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG21": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG22": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG23": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD11": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD12": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD13": {"epsilon": 0.0157, "rmin_half": 1.487}
    },
    "MET": {
        "CG": {"epsilon": 0.1094, "rmin_half": 1.908},
        "SD": {"epsilon": 0.25, "rmin_half": 2.0},
        "CE": {"epsilon": 0.1094, "rmin_half": 1.908},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG2": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HG3": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HE1": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HE2": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HE3": {"epsilon": 0.0157, "rmin_half": 1.387}
    },
    "PHE": {
        "CG": {"epsilon": 0.086, "rmin_half": 1.908},
        "CD1": {"epsilon": 0.086, "rmin_half": 1.908},
        "CD2": {"epsilon": 0.086, "rmin_half": 1.908},
        "CE1": {"epsilon": 0.086, "rmin_half": 1.908},
        "CE2": {"epsilon": 0.086, "rmin_half": 1.908},
        "CZ": {"epsilon": 0.086, "rmin_half": 1.908},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD1": {"epsilon": 0.015, "rmin_half": 1.459},
        "HD2": {"epsilon": 0.015, "rmin_half": 1.459},
        "HE1": {"epsilon": 0.015, "rmin_half": 1.459},
        "HE2": {"epsilon": 0.015, "rmin_half": 1.459},
        "HZ": {"epsilon": 0.015, "rmin_half": 1.459}
    },
    "PRO": {
        "CG": {"epsilon": 0.1094, "rmin_half": 1.908},
        "CD": {"epsilon": 0.1094, "rmin_half": 1.908},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD1": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HD2": {"epsilon": 0.0157, "rmin_half": 1.387}
    },
    "SER": {
        "OG": {"epsilon": 0.2104, "rmin_half": 1.721},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HG": {"epsilon": 5.27e-05, "rmin_half": 0.3}
    },
    "THR": {
        "OG1": {"epsilon": 0.2104, "rmin_half": 1.721},
        "CG2": {"epsilon": 0.1094, "rmin_half": 1.908},
        "HB": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HG1": {"epsilon": 5.27e-05, "rmin_half": 0.3},
        "HG21": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG22": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG23": {"epsilon": 0.0157, "rmin_half": 1.487}
    },
    "CYS": {
        "SG": {"epsilon": 0.25, "rmin_half": 2.0},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HG": {"epsilon": 0.0157, "rmin_half": 0.6}
    },
    "CSS": {
        "SG": {"epsilon": 0.25, "rmin_half": 2.0},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.387}
    },
    "CYH": {
        "SG": {"epsilon": 0.25, "rmin_half": 2.0},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HG": {"epsilon": 0.0157, "rmin_half": 0.6}
    },
    "ASN": {
        "CG": {"epsilon": 0.086, "rmin_half": 1.908},
        "OD1": {"epsilon": 0.21, "rmin_half": 1.6612},
        "ND2": {"epsilon": 0.17, "rmin_half": 1.824},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD21": {"epsilon": 0.0157, "rmin_half": 0.6},
        "HD22": {"epsilon": 0.0157, "rmin_half": 0.6}
    },
    "GLN": {
        "CG": {"epsilon": 0.1094, "rmin_half": 1.908},
        "CD": {"epsilon": 0.086, "rmin_half": 1.908},
        "OE1": {"epsilon": 0.21, "rmin_half": 1.6612},
        "NE2": {"epsilon": 0.17, "rmin_half": 1.824},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HE21": {"epsilon": 0.0157, "rmin_half": 0.6},
        "HE22": {"epsilon": 0.0157, "rmin_half": 0.6}
    },
    "TYR": {
        "CG": {"epsilon": 0.086, "rmin_half": 1.908},
        "CD1": {"epsilon": 0.086, "rmin_half": 1.908},
        "CD2": {"epsilon": 0.086, "rmin_half": 1.908},
        "CE1": {"epsilon": 0.086, "rmin_half": 1.908},
        "CE2": {"epsilon": 0.086, "rmin_half": 1.908},
        "CZ": {"epsilon": 0.086, "rmin_half": 1.908},
        "OH": {"epsilon": 0.2104, "rmin_half": 1.721},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD1": {"epsilon": 0.015, "rmin_half": 1.459},
        "HD2": {"epsilon": 0.015, "rmin_half": 1.459},
        "HE1": {"epsilon": 0.015, "rmin_half": 1.459},
        "HE2": {"epsilon": 0.015, "rmin_half": 1.459},
        "HH": {"epsilon": 5.27e-05, "rmin_half": 0.3}
    },
    "TRP": {
        "CG": {"epsilon": 0.086, "rmin_half": 1.908},
        "CD1": {"epsilon": 0.086, "rmin_half": 1.908},
        "CD2": {"epsilon": 0.086, "rmin_half": 1.908},
        "NE1": {"epsilon": 0.17, "rmin_half": 1.824},
        "CE2": {"epsilon": 0.086, "rmin_half": 1.908},
        "CE3": {"epsilon": 0.086, "rmin_half": 1.908},
        "CZ2": {"epsilon": 0.086, "rmin_half": 1.908},
        "CZ3": {"epsilon": 0.086, "rmin_half": 1.908},
        "CH2": {"epsilon": 0.086, "rmin_half": 1.908},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD1": {"epsilon": 0.015, "rmin_half": 1.409},
        "HE1": {"epsilon": 0.0157, "rmin_half": 0.6},
        "HE3": {"epsilon": 0.015, "rmin_half": 1.459},
        "HZ2": {"epsilon": 0.015, "rmin_half": 1.459},
        "HZ3": {"epsilon": 0.015, "rmin_half": 1.459},
        "HH2": {"epsilon": 0.015, "rmin_half": 1.459},
    },
    "ASP": {
        "CG": {"epsilon": 0.086, "rmin_half": 1.908},
        "OD1": {"epsilon": 0.21, "rmin_half": 1.6612},
        "OD2": {"epsilon": 0.21, "rmin_half": 1.6612},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
    },
    "GLU": {
        "CG": {"epsilon": 0.1094, "rmin_half": 1.908},
        "CD": {"epsilon": 0.086, "rmin_half": 1.908},
        "OE1": {"epsilon": 0.21, "rmin_half": 1.6612},
        "OE2": {"epsilon": 0.21, "rmin_half": 1.6612},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG3": {"epsilon": 0.0157, "rmin_half": 1.487},
    },
    "HIS": {
        "CG": {"epsilon": 0.086, "rmin_half": 1.908},
        "CD2": {"epsilon": 0.086, "rmin_half": 1.908},
        "ND1": {"epsilon": 0.17, "rmin_half": 1.824},
        "CE1": {"epsilon": 0.086, "rmin_half": 1.908},
        "NE2": {"epsilon": 0.17, "rmin_half": 1.824},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD1": {"epsilon": 0.0157, "rmin_half": 0.6},
        "HD2": {"epsilon": 0.015, "rmin_half": 1.409},
        "HE1": {"epsilon": 0.015, "rmin_half": 1.409},
        "HE2": {"epsilon": 0.0157, "rmin_half": 0.6}
    },
    "HID": {
        "CG": {"epsilon": 0.086, "rmin_half": 1.908},
        "CD2": {"epsilon": 0.086, "rmin_half": 1.908},
        "ND1": {"epsilon": 0.17, "rmin_half": 1.824},
        "CE1": {"epsilon": 0.086, "rmin_half": 1.908},
        "NE2": {"epsilon": 0.17, "rmin_half": 1.824},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD1": {"epsilon": 0.0157, "rmin_half": 0.6},
        "HD2": {"epsilon": 0.015, "rmin_half": 1.409},
        "HE1": {"epsilon": 0.015, "rmin_half": 1.409}
    },
    "HIE": {
        "CG": {"epsilon": 0.086, "rmin_half": 1.908},
        "CD2": {"epsilon": 0.086, "rmin_half": 1.908},
        "ND1": {"epsilon": 0.17, "rmin_half": 1.824},
        "CE1": {"epsilon": 0.086, "rmin_half": 1.908},
        "NE2": {"epsilon": 0.17, "rmin_half": 1.824},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD2": {"epsilon": 0.015, "rmin_half": 1.409},
        "HE1": {"epsilon": 0.015, "rmin_half": 1.409},
        "HE2": {"epsilon": 0.0157, "rmin_half": 0.6}
    },
    "HIP": {
        "CG": {"epsilon": 0.086, "rmin_half": 1.908},
        "CD2": {"epsilon": 0.086, "rmin_half": 1.908},
        "ND1": {"epsilon": 0.17, "rmin_half": 1.824},
        "CE1": {"epsilon": 0.086, "rmin_half": 1.908},
        "NE2": {"epsilon": 0.17, "rmin_half": 1.824},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD1": {"epsilon": 0.0157, "rmin_half": 0.6},
        "HD2": {"epsilon": 0.015, "rmin_half": 1.409},
        "HE1": {"epsilon": 0.015, "rmin_half": 1.409},
        "HE2": {"epsilon": 0.0157, "rmin_half": 0.6}
    },
    "LYS": {
        "CG": {"epsilon": 0.1094, "rmin_half": 1.908},
        "CD": {"epsilon": 0.1094, "rmin_half": 1.908},
        "CE": {"epsilon": 0.1094, "rmin_half": 1.908},
        "NZ": {"epsilon": 0.17, "rmin_half": 1.824},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HE2": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HE3": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HZ1": {"epsilon": 0.17, "rmin_half": 1.824},
        "HZ2": {"epsilon": 0.17, "rmin_half": 1.824},
        "HZ3": {"epsilon": 0.17, "rmin_half": 1.824}
    },
    "ARG": {
        "CG": {"epsilon": 0.1094, "rmin_half": 1.908},
        "CD": {"epsilon": 0.1094, "rmin_half": 1.908},
        "NE": {"epsilon": 0.17, "rmin_half": 1.824},
        "CZ": {"epsilon": 0.086, "rmin_half": 1.908},
        "NH1": {"epsilon": 0.17, "rmin_half": 1.824},
        "NH2": {"epsilon": 0.17, "rmin_half": 1.824},
        "HB2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HB3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG2": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HG3": {"epsilon": 0.0157, "rmin_half": 1.487},
        "HD2": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HD3": {"epsilon": 0.0157, "rmin_half": 1.387},
        "HE": {"epsilon": 0.17, "rmin_half": 1.824},
        "HH11": {"epsilon": 0.17, "rmin_half": 1.824},
        "HH12": {"epsilon": 0.17, "rmin_half": 1.824},
        "HH21": {"epsilon": 0.17, "rmin_half": 1.824},
        "HH22": {"epsilon": 0.17, "rmin_half": 1.824}
    },
}


LENNARD_JONES_12_6: List[SmirnoffParameter] = [
    {"smirks": "[#1:1]", "epsilon": 0.0157, "rmin_half": 0.6},
    {"smirks": "[#1:1]-[#6X4]", "epsilon": 0.0157, "rmin_half": 1.487},
    {"smirks": "[#1:1]-[#6X4]-[#7,#8,#9,#16,#17,#35]", "epsilon": 0.0157, "rmin_half": 1.387},
    {"smirks": "[#1:1]-[#6X4](-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]", "epsilon": 0.0157, "rmin_half": 1.287},
    {"smirks": "[#1:1]-[#6X4](-[#7,#8,#9,#16,#17,#35])(-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]", "epsilon": 0.0157, "rmin_half": 1.187},
    {"smirks": "[#1:1]-[#6X4]~[*+1,*+2]", "epsilon": 0.0157, "rmin_half": 1.1},
    {"smirks": "[#1:1]-[#6X3]", "epsilon": 0.015, "rmin_half": 1.459},
    {"smirks": "[#1:1]-[#6X3]~[#7,#8,#9,#16,#17,#35]", "epsilon": 0.015, "rmin_half": 1.409},
    {"smirks": "[#1:1]-[#6X3](~[#7,#8,#9,#16,#17,#35])~[#7,#8,#9,#16,#17,#35]", "epsilon": 0.015, "rmin_half": 1.359},
    {"smirks": "[#1:1]-[#6X2]", "epsilon": 0.015, "rmin_half": 1.459},
    {"smirks": "[#1:1]-[#7]", "epsilon": 0.0157, "rmin_half": 0.6},
    {"smirks": "[#1:1]-[#8]", "epsilon": 5.27e-05, "rmin_half": 0.3},
    {"smirks": "[#1:1]-[#16]", "epsilon": 0.0157, "rmin_half": 0.6},
    {"smirks": "[#3+1:1]", "epsilon": 0.0279896, "rmin_half": 1.025},
    {"smirks": "[#6:1]", "epsilon": 0.086, "rmin_half": 1.908},
    {"smirks": "[#6X2:1]", "epsilon": 0.21, "rmin_half": 1.908},
    {"smirks": "[#6X4:1]", "epsilon": 0.1094, "rmin_half": 1.908},
    {"smirks": "[#7:1]", "epsilon": 0.17, "rmin_half": 1.824},
    {"smirks": "[#8:1]", "epsilon": 0.21, "rmin_half": 1.6612},
    {"smirks": "[#8X2H0+0:1]", "epsilon": 0.17, "rmin_half": 1.6837},
    {"smirks": "[#8X2H1+0:1]", "epsilon": 0.2104, "rmin_half": 1.721},
    {"smirks": "[#9:1]", "epsilon": 0.061, "rmin_half": 1.75},
    {"smirks": "[#9X0-1:1]", "epsilon": 0.003364, "rmin_half": 2.303},
    {"smirks": "[#11+1:1]", "epsilon": 0.0874393, "rmin_half": 1.369},
    {"smirks": "[#15:1]", "epsilon": 0.2, "rmin_half": 2.1},
    {"smirks": "[#16:1]", "epsilon": 0.25, "rmin_half": 2.0},
    {"smirks": "[#17:1]", "epsilon": 0.265, "rmin_half": 1.948},
    {"smirks": "[#17X0-1:1]", "epsilon": 0.035591, "rmin_half": 2.513},
    {"smirks": "[#19+1:1]", "epsilon": 0.1936829, "rmin_half": 1.705},
    {"smirks": "[#35:1]", "epsilon": 0.32, "rmin_half": 2.22},
    {"smirks": "[#35X0-1:1]", "epsilon": 0.0586554, "rmin_half": 2.608},
    {"smirks": "[#37+1:1]", "epsilon": 0.3278219, "rmin_half": 1.813},
    {"smirks": "[#53:1]", "epsilon": 0.4, "rmin_half": 2.35},
    {"smirks": "[#53X0-1:1]", "epsilon": 0.0536816, "rmin_half": 2.86},
    {"smirks": "[#55+1:1]", "epsilon": 0.4065394, "rmin_half": 1.976},
]
"""
Lennard-Jones-12-6 potential.

Values from
https://github.com/openforcefield/smirnoff99Frosst/blob/master/smirnoff99frosst/offxml/smirnoff99Frosst-1.1.0.offxml#L305
"""
