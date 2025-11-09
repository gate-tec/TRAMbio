import os
from TRAMbio.services.parameter._base_parameter import BaseParameter
from TRAMbio.util.errors import SafeDependencyTest, DependencyWithSafeTest


class GeneralWorkflowParameter(BaseParameter):
    VERBOSE = ("TRAM_VERBOSE",
               "Verbosity switch for various outputs like progress bars (default: True)",
               True, None)

class PebbleGameParameter(BaseParameter):
    """Parameter for Pebble Game workflow"""
    K = ("TRAM_PEBBLE_GAME_K",
         "Parameter k for the pebble game",
         2, lambda x: x >= 1)
    L = ("TRAM_PEBBLE_GAME_L",
         "Parameter l for the pebble game",
         3, lambda x: x >= 0)
    THREADS = ("TRAM_PEBBLE_GAME_THREADS",
               "Number of threads for multiprocessing. A value of 1 specifies no multiprocessing (default: 1)",
               1, lambda x: 1 <= x <= os.cpu_count())

class XtcParameter(BaseParameter):
    """Parameter for XTC workflow"""
    MODULE = ("TRAM_XTC_MODULE",
              "The selected third-party module for loading XTC files",
              DependencyWithSafeTest.MDAnalysis.value
              if SafeDependencyTest.is_present(DependencyWithSafeTest.MDAnalysis) else
              (DependencyWithSafeTest.mdtraj.value
               if SafeDependencyTest.is_present(DependencyWithSafeTest.mdtraj) else
               DependencyWithSafeTest.MDAnalysis.value),
              None)
    STRIDE = ("TRAM_XTC_STRIDE",
              "Selected stride for XTC frames (default: 50 frames)",
              50, lambda x: x >= 1)
    DYNAMIC_SCALING = ("TRAM_XTC_DYNAMIC_SCALING",
                       "Whether to activate dynamic core allocation during multiprocessing setup. If False, 2/3 of the available cores will be used for graph construction and 1/3 for the pebble game runs (default: True)",
                       True, None)

class ResidueParameter(BaseParameter):
    """Parameter for Residue workflow"""
    MIN_KEY = ("TRAM_RESIDUE_MIN_KEY",
               "Minimum value for state keys. Either a float for single frame PDBs, indicating the minimum strength for present hydrogen bonds, or a starting frame number for trajectories (default: no limit)",
               None, None)
    MAX_STATES = ("TRAM_RESIDUE_MAX_STATES",
                  "Maximum number of states to visualize. Values lower than 1 indicate no limit. (default: 0)",
                  0, lambda x: x >= 0)
    THRESHOLD = ("TRAM_RESIDUE_THRESHOLD",
                 "Percentage of required residue atoms within a component to count it as present on the residue level. Requires PDB data to be provided (default: 0.8)",
                 0.8, lambda x: 0 <= x <= 1)
    USE_MAIN_CHAIN = ("TRAM_RESIDUE_USE_MAIN_CHAIN",
                      "Whether to consider a residue present in the residue-level component if the N-CA-C main chain atoms are contained (default: True)",
                      True, None)

class PyMolParameter(BaseParameter):
    """Parameter for PyMol workflow"""
    ALL_WEIGHTED_BONDS = ("TRAM_PYMOL_ALL_WEIGHTED_BONDS",
                          "Whether to include all weighted interactions (from a provided bond file) in the PyMol frames or only hydrogen bonds and salt-bridges (PDB to PyMol only, default: False)",
                          False, None)
