import re
from io import StringIO
from typing import Union, Optional, Dict

import pandas as pd

from TRAMbio.services import InteractionServiceRegistry, InteractionServiceException
from TRAMbio.services.parameter import ParameterRegistry, lock_registry, PdbParameter
from TRAMbio.services.structure import IPdbStructureService, StructureServiceRegistry
from TRAMbio.services.structure.util import create_base_graphs, export_ter_flagged_residues
from TRAMbio.util.constants.interaction import INTERACTION_RANKING, map_interaction_ranking
from TRAMbio.util.constants.pdb import ATOM_DF_COLUMNS, HEADER_RECORDS
from TRAMbio.util.wrapper.biopandas.pandas_pdb import CustomPandasPdb
from TRAMbio.util.structure_library.graph_struct import ProteinGraph as ProteinGraph, copy_graphs_for_dataframe, GraphKey
from loguru import logger as logger


__all__ = []


class PdbStructureService(IPdbStructureService):

    @property
    def name(self):
        return 'PdbStructureService'

    @lock_registry(kwargs_name="parameter_id")
    def export_atom_df(
            self,
            raw_df: Dict[str, pd.DataFrame],
            check_ids: bool = False,
            parameter_id: str = ''
    ) -> pd.DataFrame:
        keep_hets = ParameterRegistry.get_parameter_set(parameter_id)(PdbParameter.KEEP_HETS.value)
        atom_df = pd.concat([raw_df['ATOM'], raw_df['HETATM']], ignore_index=True) if keep_hets else raw_df['ATOM']

        if len(atom_df) < 1:
            logger.warning("No atom coordinates in data frame.")
            return pd.DataFrame([], columns=ATOM_DF_COLUMNS + ["node_id"])

        # limit columns
        atom_df = atom_df.loc[:, ATOM_DF_COLUMNS].reset_index(drop=True)

        # label node ids
        atom_df["node_id"] = (
                atom_df["chain_id"].apply(str)
                + atom_df["residue_number"].apply(lambda x: str(x).rjust(4, '0'))
                + atom_df["insertion"].apply(lambda x: str(x) if len(x) > 0 else '-')
                + atom_df["residue_name"]
                + ":"
                + atom_df["atom_name"]
        )
        if check_ids:
            self._check_for_duplicated_node_ids(atom_df=atom_df)

        return atom_df

    @lock_registry(kwargs_name="parameter_id")
    def has_hydrogen_atoms(
            self,
            raw_or_atom_df: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
            parameter_id: str = ''
    ) -> bool:
        if isinstance(raw_or_atom_df, pd.DataFrame):
            df = raw_or_atom_df
        else:
            df = raw_or_atom_df['ATOM']
        return len(df) > 0 and any(df.element_symbol == 'H')  # noqa comparing against a pandas series returns a series

    @lock_registry(kwargs_name="parameter_id")
    def export_others_df(
            self,
            raw_df: Dict[str, pd.DataFrame],
            ter_only: bool = False,
            parameter_id: str = ''
    ) -> pd.DataFrame:
        others_df = raw_df['OTHERS']
        if ter_only:
            others_df = others_df.loc[others_df["record_name"] == 'TER', :]
        return others_df.reset_index(drop=True)

    @lock_registry(kwargs_name="parameter_id")
    def export_header_stream(
            self,
            raw_df: Dict[str, pd.DataFrame],
            pdb_name: Optional[str] = None,
            parameter_id: str = ''
    ) -> StringIO:
        header_frame = None
        if len(raw_df['OTHERS']) > 0:
            header_frame = raw_df['OTHERS'].loc[
                           raw_df['OTHERS']["record_name"].isin(HEADER_RECORDS), :].reset_index(drop=True)

        if header_frame is None or len(header_frame) == 0:
            new_title = re.sub("[^A-Z0-9 ]", "", pdb_name.upper().replace("_", " ") if pdb_name is not None else "PROTEIN")

            parts = []
            while len(new_title) > 70:
                parts.append(new_title[:70])
                new_title = " " + new_title[70:]
            parts.append(new_title)

            header_lines = [
                ["TITLE", f"    {parts[0]}", 0]
            ]
            for i in range(1, len(parts)):
                header_lines.append(
                    ["TITLE", f" {i+1:3d}{parts[i]}", i]
                )

            header_frame = pd.DataFrame(header_lines, columns=['record_name', 'entry', 'line_idx'])

        header_pdb_df = CustomPandasPdb()
        header_pdb_df.df['OTHERS'] = header_frame
        header_stream = header_pdb_df.to_pdb_stream(records=['OTHERS'])

        return header_stream

    @staticmethod
    def _check_for_duplicated_node_ids(
            atom_df: pd.DataFrame,
    ):
        duplicate_ids = atom_df.loc[atom_df.duplicated(subset=['node_id'], keep=False), :]
        if len(duplicate_ids) > 0:
            logger.error(
                f"Duplicate atom names within residues:\n{duplicate_ids[['atom_number', 'chain_id', 'residue_number', 'residue_name', 'atom_name', 'node_id']].to_string()}")
            raise KeyError("Duplicate atom names within residues.")

    @lock_registry(kwargs_name="parameter_id")
    def create_graph_struct(
            self,
            atom_df: pd.DataFrame,
            others_df: pd.DataFrame,
            parameter_id: str = ''
    ) -> ProteinGraph:
        # subset heavy atoms
        hydrogen_df = atom_df.copy().loc[atom_df.element_symbol == "H", :].reset_index(drop=True)
        heavy_atom_df = atom_df.copy().loc[
                        (atom_df["record_name"].str.strip() == "HETATM") | (atom_df["element_symbol"] != "H"),
                        :].reset_index(drop=True)

        # get TER flagged residues
        ter_flags = export_ter_flagged_residues(atom_df=atom_df, others_df=others_df)

        # construct graphs and insert covalent edges within residues and peptide bonds
        graphs, hydrogen_mapping = create_base_graphs(heavy_atom_df=heavy_atom_df, all_atom_df=atom_df,
                                                      ter_flags=ter_flags)

        protein_graph = ProteinGraph(
            graphs=graphs,
            atom_df=atom_df, others_df=others_df,
            heavy_atom_df=heavy_atom_df, hydrogen_df=hydrogen_df,
            hydrogen_mapping=hydrogen_mapping
        )

        for interaction_service in sorted(
                InteractionServiceRegistry.COV.list_services(),
                key=lambda x: map_interaction_ranking(x.interaction_types[0])
        ):
            try:
                interaction_service.apply_interactions(protein_graph, parameter_id)
            except InteractionServiceException as iE:
                logger.error(str(iE))

        return protein_graph

    @lock_registry(kwargs_name="parameter_id")
    def copy_graph_for_frame(
            self,
            atom_df: pd.DataFrame,
            others_df: pd.DataFrame,
            protein_graph: ProteinGraph,
            parameter_id: str = ''
    ) -> ProteinGraph:
        # subset heavy atoms
        hydrogen_df = atom_df.copy().loc[atom_df.element_symbol == "H", :].reset_index(drop=True)
        heavy_atom_df = atom_df.copy().loc[
                        (atom_df["record_name"].str.strip() == "HETATM") | (atom_df["element_symbol"] != "H"),
                        :].reset_index(drop=True)

        return ProteinGraph(
            graphs=copy_graphs_for_dataframe(protein_graph.graphs, atom_df),
            atom_df=atom_df, others_df=others_df,
            heavy_atom_df=heavy_atom_df, hydrogen_df=hydrogen_df,
            hydrogen_mapping=protein_graph.hydrogen_mapping
        )

    @lock_registry(kwargs_name="parameter_id")
    def apply_non_covalent_interactions(
            self,
            protein_graph: ProteinGraph,
            parameter_id: str = ''
    ) -> None:
        for interaction_service in sorted(
                InteractionServiceRegistry.NON_COV.list_services(),
                key=lambda x: map_interaction_ranking(x.interaction_types[0])
        ):
            try:
                interaction_service.apply_interactions(protein_graph, parameter_id)
            except InteractionServiceException as iE:
                logger.error(str(iE))

        protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value] = \
            sorted(protein_graph.graphs['pebble'].graph[GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value], key=lambda x: x[3])

        if ParameterRegistry.get_parameter_set(parameter_id)(PdbParameter.UNIQUE_BONDS.value):
            for _, _, d in protein_graph.graphs['atom'].edges(data=True):
                d['kind'] = {sorted([(INTERACTION_RANKING[bond], bond) for bond in d['kind']])[0][1]}
            for _, _, d in protein_graph.graphs['full'].edges(data=True):
                d['kind'] = {sorted([(INTERACTION_RANKING[bond], bond) for bond in d['kind']])[0][1]}


StructureServiceRegistry.PDB.register_service(PdbStructureService())
