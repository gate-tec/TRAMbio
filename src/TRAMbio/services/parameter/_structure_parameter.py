from TRAMbio.services.parameter._base_parameter import BaseParameter


class PdbParameter(BaseParameter):
    UNIQUE_BONDS = ("TRAM_PDB_UNIQUE_BONDS",
                    "Whether to limit annotations for atomic edges to a single, unique label (default: False)",
                    False, None)
    KEEP_HETS = ("TRAM_PDB_KEEP_HETS",
                 "Whether to include HETATM records from PDB data (default: True)",
                 True, None)
