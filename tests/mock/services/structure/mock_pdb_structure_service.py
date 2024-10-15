from io import StringIO
from typing import Dict, Optional, Union

import pandas as pd

from TRAMbio.services.structure import IPdbStructureService
from TRAMbio.util.structure_library.graph_struct import ProteinGraph


class MockHydrogenTestPdbStructureService(IPdbStructureService):

    def __init__(self, has_hydrogen: bool = True):
        self.has_hydrogen = has_hydrogen

    @property
    def name(self):
        return "PdbStructureService"

    def export_atom_df(self, raw_df: pd.DataFrame, check_ids: bool = False, parameter_id: str = '') -> pd.DataFrame:
        raise NotImplementedError

    def has_hydrogen_atoms(self, raw_or_atom_df: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
                           parameter_id: str = '') -> bool:
        return self.has_hydrogen

    def export_others_df(self, raw_df: Dict[str, pd.DataFrame], ter_only: bool = False,
                         parameter_id: str = '') -> pd.DataFrame:
        raise NotImplementedError

    def export_header_stream(self, raw_df: Dict[str, pd.DataFrame], pdb_name: Optional[str] = None,
                             parameter_id: str = '') -> StringIO:
        raise NotImplementedError

    def create_graph_struct(self, atom_df: pd.DataFrame, others_df: pd.DataFrame,
                            parameter_id: str = '') -> ProteinGraph:
        raise NotImplementedError

    def copy_graph_for_frame(self, atom_df: pd.DataFrame, others_df: pd.DataFrame, protein_graph: ProteinGraph,
                             parameter_id: str = '') -> ProteinGraph:
        raise NotImplementedError

    def apply_non_covalent_interactions(self, protein_graph: ProteinGraph, parameter_id: str = '') -> None:
        raise NotImplementedError

