from TRAMbio.services.parameter._base_parameter import BaseParameter


class HydrogenBondParameter(BaseParameter):
    INCLUDE = ("TRAM_HYDROGEN_INCLUDE",
               "Whether to insert hydrogen bonds (default: True)",
               True, None)
    MINIMAL_LENGTH = ("TRAM_HYDROGEN_MINIMAL_LENGTH",
                      "Minimum distance between atoms to consider for a hydrogen bond (in Angstroms, default: 2.6)",
                      2.6, lambda x: x >= 0.0)
    ENERGY_THRESHOLD = ("TRAM_HYDROGEN_ENERGY_THRESHOLD",
                        "Energy threshold for inclusion of hydrogen bonds. All bonds with energy lower or equal to this threshold are included (default: -0.1)",
                        -0.1, None)
    CUTOFF_DISTANCE = ("TRAM_HYDROGEN_CUTOFF_DISTANCE",
                       "Maximum distance between atoms to consider for a hydrogen bond (in Angstroms, default: 3.0)",
                       3.0, lambda x: x > 0.0)
    STRONG_ENERGY_THRESHOLD = ("TRAM_HYDROGEN_STRONG_ENERGY_THRESHOLD",
                               "Separate energy threshold for 'strong' hydrogen bonds. Every bond with energies greater than this will be considered 'weak' and modeled with a gradually fewer number of bars (default: 0.0)",
                               0.0, None)
    BAR_COUNT = ("TRAM_HYDROGEN_BAR_COUNT",
                 "Number of bars modeling a hydrogen bond in the pebble game (default: 5)",
                 5, lambda x: 1 <= x <= 6)


class HydrophobicInteractionParameter(BaseParameter):
    INCLUDE = ("TRAM_HYDROPHOBIC_INCLUDE",
               "Whether to insert hydrophobic interactions (default: True)",
               True, None)
    MINIMAL_LENGTH = ("TRAM_HYDROPHOBIC_MINIMAL_LENGTH",
                      "Whether to only keep the hydrophobic interaction with the shortest length for each atom (default: True)",
                      True, None)
    POTENTIAL = ("TRAM_HYDROPHOBIC_USE_POTENTIAL",
                 "Whether to use hydrophobic potential calculation instead of surface calculation (default: False)",
                 False, None)
    SURFACE_CUTOFF_DISTANCE = ("TRAM_HYDROPHOBIC_SURFACE_CUTOFF_DISTANCE",
                               "Maximum distance between the VDW surfaces of atoms in hydrophobic interactions (Surface only, in Angstroms, default: 0.25)",
                               0.25, lambda x: x > 0.0)
    POTENTIAL_CUTOFF_DISTANCE = ("TRAM_HYDROPHOBIC_POTENTIAL_CUTOFF_DISTANCE",
                                 "Maximum distance between atoms to consider for hydrophobic interactions (Potential only, in Angstroms, default: 9.0)",
                                 9.0, lambda x: x > 0.0)
    SCALE_14 = ("TRAM_HYDROPHOBIC_SCALE_14",
                "Scaling factor for energy potential in hydrophobic interaction between 3rd-degree (1-4) covalent neighbors (Potential only, default: 0.5)",
                0.5, lambda x: x > 0.0)
    SCALE_15 = ("TRAM_HYDROPHOBIC_SCALE_15",
                "Scaling factor for energy potential in hydrophobic interaction between 4th-degree (1-5) covalent neighbors (Potential only, default: 1.0)",
                1.0, lambda x: x > 0.0)
    SCALE_UNBOUNDED = ("TRAM_HYDROPHOBIC_SCALE_UNBOUNDED",
                       "Scaling factor for energy potential in hydrophobic interaction between 5th or higher degree covalent neighbors (Potential only, default: 0.5)",
                       1.0, lambda x: x > 0.0)
    ENERGY_THRESHOLD = ("TRAM_HYDROPHOBIC_ENERGY_THRESHOLD",
                        "Energy threshold for inclusion of hydrophobic interactions. All interactions with energy lower or equal to this threshold are included (Potential only, default: -0.1)",
                        -0.1, None)
    BAR_COUNT = ("TRAM_HYDROPHOBIC_BAR_COUNT",
                 "Number of bars modeling a hydrophobic interaction in the pebble game (default: 3)",
                 3, lambda x: 1 <= x <= 5)


class DisulphideBridgeParameter(BaseParameter):
    INCLUDE = ("TRAM_DISULPHIDE_INCLUDE",
               "Whether to insert disulphide bridges (default: True)",
               True, None)
    CUTOFF_DISTANCE = ("TRAM_DISULPHIDE_CUTOFF_DISTANCE",
                       "Maximum distance between sulphur atoms to consider for a disulphide bridge(in Angstroms, default: 3.0)",
                       3.0, lambda x: x > 0.0)


class CationPiInteractionParameter(BaseParameter):
    INCLUDE = ("TRAM_CATION_PI_INCLUDE",
               "Whether to insert cation-pi interactions (default: True)",
               True, None)
    CUTOFF_DISTANCE = ("TRAM_CATION_PI_CUTOFF_DISTANCE",
                       "Maximum distance between atoms to consider for a cation-pi interaction (in Angstroms, default: 3.0)",
                       6.0, lambda x: x > 0.0)
    BAR_COUNT = ("TRAM_CATION_PI_BAR_COUNT",
                 "Number of bars modeling a cation-pi interaction in the pebble game (default: 3)",
                 3, lambda x: 1 <= x <= 5)


class AromaticInteractionParameter(BaseParameter):
    INCLUDE = ("TRAM_AROMATIC_INCLUDE",
               "Whether to insert aromatic interactions (default: True)",
               True, None)
    CUTOFF_DISTANCE_PI = ("TRAM_AROMATIC_CUTOFF_DISTANCE_PI",
                          "Maximum distance between aromatic centers to consider for a pi-stacking (in Angstroms, default: 7.0)",
                          7.0, lambda x: x > 0.0)
    CUTOFF_DISTANCE_T = ("TRAM_AROMATIC_CUTOFF_DISTANCE_T",
                         "Maximum distance between aromatic centers to consider for a t-stacking (in Angstroms, default: 5.0)",
                         5.0, lambda x: x > 0.0)
    ANGLE_VARIANCE = ("TRAM_AROMATIC_ANGLE_VARIANCE",
                      "Maximum allowed variance of the interaction angle in aromatic interactions (in degrees, default: 5.0)",
                      5.0, lambda x: x > 0.0)
    BAR_COUNT = ("TRAM_AROMATIC_BAR_COUNT",
                 "Number of bars modeling an aromatic interaction in the pebble game (default: 3)",
                 3, lambda x: 1 <= x <= 5)


class PdbEntryInteractionParameter(BaseParameter):
    SSBOND_INCLUDE = ("TRAM_SSBOND_INCLUDE",
               "Whether to insert SSBOND records from PDB data (default: True)",
               True, None)
    LINK_INCLUDE = ("TRAM_LINK_INCLUDE",
               "Whether to insert LINK records from PDB data (default: True)",
               True, None)
    CONECT_INCLUDE = ("TRAM_CONECT_INCLUDE",
               "Whether to insert CONECT records from PDB data (default: True)",
               True, None)
