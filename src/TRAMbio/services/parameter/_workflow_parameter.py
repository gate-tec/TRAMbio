from enum import Enum


class GeneralWorkflowParameter(Enum):
    VERBOSE = "TRAM_VERBOSE"

class PebbleGameParameter(Enum):
    """Parameter for Pebble Game workflow"""
    K = "TRAM_PEBBLE_GAME_K"
    L = "TRAM_PEBBLE_GAME_L"
    THREADS = "TRAM_PEBBLE_GAME_THREADS"

class XtcParameter(Enum):
    """Parameter for XTC workflow"""
    MODULE = "TRAM_XTC_MODULE"
    STRIDE = "TRAM_XTC_STRIDE"
    DYNAMIC_SCALING = "TRAM_XTC_DYNAMIC_SCALING"

class ResidueParameter(Enum):
    """Parameter for Residue workflow"""
    MIN_KEY = "TRAM_RESIDUE_MIN_KEY"
    MAX_STATES = "TRAM_RESIDUE_MAX_STATES"
    THRESHOLD = "TRAM_RESIDUE_THRESHOLD"
    USE_MAIN_CHAIN = "TRAM_RESIDUE_USE_MAIN_CHAIN"

class PyMolParameter(Enum):
    """Parameter for PyMol workflow"""
    ALL_WEIGHTED_BONDS = "TRAM_PYMOL_ALL_WEIGHTED_BONDS"
