from typing import Dict


#################
# Constants #####
#################

COVALENT_RADII: Dict[str, float] = {
    "Csb": 0.77,
    "Cres": 0.72,
    "Cdb": 0.67,
    "Osb": 0.67,
    "Ores": 0.635,
    "Odb": 0.60,
    "Nsb": 0.70,
    "Nres": 0.66,
    "Ndb": 0.62,
    "Hsb": 0.37,
    "Ssb": 1.04,
    "Psb": 1.11,
    "Pres": 1.06,
    "Pdb": 1.02
}
"""
Covalent radii for OpenSCAD output.
Adding ``Ores`` between ``Osb`` and ``Odb`` for ``Asp`` and ``Glu``, ``Nres`` between ``Nsb`` and ``Ndb`` for ``Arg``, as PDB does not specify.
Adding ``Pres`` between ``Psb`` and ``Pdb``.

Covalent radii from:

    Heyrovska (2008) : 'Atomic Structures of all the Twenty Essential Amino Acids and a Tripeptide, with Bond Lengths as Sums of Atomic Covalent Radii'

Paper: https://arxiv.org/pdf/0804.2488.pdf

Covalent radii for ``Psb`` and ``Pdb`` from:

    Pyykkö and Atsumi (2009) : 'Molecular Double-Bond Covalent Radii for Elements Li-E112'

Paper: https://doi.org/10.1002/chem.200901472
"""

DEFAULT_BOND_STATE: Dict[str, str] = {
    "N": "Nsb",
    "CA": "Csb",
    "C": "Cdb",
    "O": "Odb",
    "OXT": "Osb",
    "CB": "Csb",
    "H": "Hsb",
    # Not sure about these - assuming they're all standard Hydrogen. Won't make much difference given
    # the tolerance is larger than Hs covalent radius
    "HG1": "Hsb",
    "HE": "Hsb",
    "1HH1": "Hsb",
    "1HH2": "Hsb",
    "2HH1": "Hsb",
    "2HH2": "Hsb",
    "HG": "Hsb",
    "HH": "Hsb",
    "1HD2": "Hsb",
    "2HD2": "Hsb",
    "HZ1": "Hsb",
    "HZ2": "Hsb",
    "HZ3": "Hsb",
    "P": "Pres",
    "OP1": "Ores",
    "OP2": "Ores",
    "O5'": "Osb",
    "C5'": "Csb",
    "O4'": "Osb",
    "C4'": "Csb",
    "O3'": "Osb",
    "C3'": "Csb",
    "C1'": "Csb",
    "H5'": "Hsb",
    "H5''": "Hsb",
    "H4'": "Hsb",
    "H3'": "Hsb",
    "H1'": "Hsb",
}
"""Assignment of atom classes to atomic radii.

Covalent radii from:

    Heyrovska, Raji : 'Atomic Structures of all the Twenty Essential Amino Acids and a Tripeptide, with Bond Lengths as Sums of Atomic Covalent Radii'

Paper: https://arxiv.org/pdf/0804.2488.pdf
"""

RESIDUE_ATOM_BOND_STATE: Dict[str, Dict[str, str]] = {
    "XXX": {
        "N": "Nsb",
        "CA": "Csb",
        "C": "Cdb",
        "O": "Odb",
        "OXT": "Osb",
        "CB": "Csb",
        "H": "Hsb",
    },
    "VAL": {"CG1": "Csb", "CG2": "Csb"},
    "LEU": {"CG": "Csb", "CD1": "Csb", "CD2": "Csb"},
    "ILE": {"CG1": "Csb", "CG2": "Csb", "CD1": "Csb"},
    "MET": {"CG": "Csb", "SD": "Ssb", "CE": "Csb"},
    "PHE": {
        "CG": "Cdb",
        "CD1": "Cres",
        "CD2": "Cres",
        "CE1": "Cdb",
        "CE2": "Cdb",
        "CZ": "Cres",
    },
    "PRO": {"CG": "Csb", "CD": "Csb"},
    "SER": {"OG": "Osb"},
    "THR": {"OG1": "Osb", "CG2": "Csb"},
    "CYS": {"SG": "Ssb"},
    "ASN": {"CG": "Csb", "OD1": "Odb", "ND2": "Ndb"},
    "GLN": {"CG": "Csb", "CD": "Csb", "OE1": "Odb", "NE2": "Ndb"},
    "TYR": {
        "CG": "Cdb",
        "CD1": "Cres",
        "CD2": "Cres",
        "CE1": "Cdb",
        "CE2": "Cdb",
        "CZ": "Cres",
        "OH": "Osb",
    },
    "TRP": {
        "CG": "Cdb",
        "CD1": "Cdb",
        "CD2": "Cres",
        "NE1": "Nsb",
        "CE2": "Cdb",
        "CE3": "Cdb",
        "CZ2": "Cres",
        "CZ3": "Cres",
        "CH2": "Cdb",
    },
    "ASP": {"CG": "Csb", "OD1": "Ores", "OD2": "Ores"},
    "GLU": {"CG": "Csb", "CD": "Csb", "OE1": "Ores", "OE2": "Ores"},
    "HIS": {
        "CG": "Cdb",
        "CD2": "Cdb",
        "ND1": "Nsb",
        "CE1": "Cdb",
        "NE2": "Ndb",
    },
    "LYS": {"CG": "Csb", "CD": "Csb", "CE": "Csb", "NZ": "Nsb"},
    "ARG": {
        "CG": "Csb",
        "CD": "Csb",
        "NE": "Nsb",
        "CZ": "Cdb",
        "NH1": "Nres",
        "NH2": "Nres",
    },
    "A": {
        "C2'": "Csb",
        "H2'": "Hsb",
        "H2''": "Hsb",
        # Nucleo-base
        "N1": "Ndb",
        "C2": "Cdb",
        "N3": "Ndb",
        "C4": "Cdb",
        "C5": "Cdb",
        "C6": "Cdb",
        "N6": "Nsb",
        "N7": "Ndb",
        "C8": "Cdb",
        "N9": "Nsb"
    },
    "C": {
        "C2'": "Csb",
        "H2'": "Hsb",
        "H2''": "Hsb",
        # Nucleo-base
        "N1": "Nsb",
        "C2": "Cdb",
        "O2": "Odb",
        "N3": "Ndb",
        "C4": "Cdb",
        "N4": "Nsb",
        "C5": "Cdb",
        "C6": "Cdb"
    },
    "G": {
        "C2'": "Csb",
        "H2'": "Hsb",
        "H2''": "Hsb",
        # Nucleo-base
        "N1": "Nsb",
        "C2": "Cdb",
        "N2": "Nsb",
        "N3": "Ndb",
        "C4": "Cdb",
        "C5": "Cdb",
        "C6": "Cdb",
        "O6": "Odb",
        "N7": "Ndb",
        "C8": "Cdb",
        "N9": "Nsb"
    },
    "U": {
        "C2'": "Csb",
        "H2'": "Hsb",
        "H2''": "Hsb",
        # Nucleo-base
        "N1": "Nsb",
        "C2": "Cdb",
        "O2": "Odb",
        "N3": "Nsb",
        "C4": "Cdb",
        "O4": "Odb",
        "C5": "Cdb",
        "C6": "Cdb"
    },
    "DA": {
        "O2'": "Odb",
        "C2'": "Cdb",
        "HO2'": "Hsb",
        # Nucleo-base
        "N1": "Ndb",
        "C2": "Cdb",
        "N3": "Ndb",
        "C4": "Cdb",
        "C5": "Cdb",
        "C6": "Cdb",
        "N6": "Nsb",
        "N7": "Ndb",
        "C8": "Cdb",
        "N9": "Nsb"
    },
    "DC": {
        "O2'": "Odb",
        "C2'": "Cdb",
        "HO2'": "Hsb",
        # Nucleo-base
        "N1": "Nsb",
        "C2": "Cdb",
        "O2": "Odb",
        "N3": "Ndb",
        "C4": "Cdb",
        "N4": "Nsb",
        "C5": "Cdb",
        "C6": "Cdb"
    },
    "DG": {
        "O2'": "Odb",
        "C2'": "Cdb",
        "HO2'": "Hsb",
        # Nucleo-base
        "N1": "Nsb",
        "C2": "Cdb",
        "N2": "Nsb",
        "N3": "Ndb",
        "C4": "Cdb",
        "C5": "Cdb",
        "C6": "Cdb",
        "O6": "Odb",
        "N7": "Ndb",
        "C8": "Cdb",
        "N9": "Nsb"
    },
    "DT": {
        "O2'": "Odb",
        "C2'": "Cdb",
        "HO2'": "Hsb",
        # Nucleo-base
        "N1": "Nsb",
        "C2": "Cdb",
        "O2": "Odb",
        "N3": "Nsb",
        "C4": "Cdb",
        "O4": "Odb",
        "C5": "Cdb",
        "C6": "Cdb",
        "C7": "Csb"
    }
}
"""Assignment of consituent atom classes with each standard residue to atomic radii.

Covalent radii from:

    Heyrovska, Raji : 'Atomic Structures of all the Twenty Essential Amino Acids and a Tripeptide, with Bond Lengths as Sums of Atomic Covalent Radii'

Paper: https://arxiv.org/pdf/0804.2488.pdf
"""

VDW_RADII: Dict[str, float] = {
    "H": 1.2,  # 1.09
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.8,
    "S": 1.8,
    "Cl": 1.75,
    "Cu": 1.4,
}
"""van der Waals radii of the most common atoms. Taken from:

> Bondi, A. (1964). "van der Waals Volumes and Radii".
> J. Phys. Chem. 68 (3): 441–451.

https://pubs.acs.org/doi/10.1021/j100785a001
"""
