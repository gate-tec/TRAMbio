from typing import Tuple, Optional, List

import pandas as pd

from TRAMbio.services.parameter import HydrogenBondParameter
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.services.workflow import WorkflowServiceRegistry, BaseWorkflowService, IPdbWorkflowService
from TRAMbio.services import StructureServiceRegistry, IOServiceRegistry, ParameterRegistry, lock_registry
from TRAMbio.services.workflow.util import apply_pebble_game
from TRAMbio.util.structure_library.components import PebbleGameResult
from TRAMbio.util.structure_library.graph_struct import export_bond_frame, ProteinGraph
from TRAMbio.util.structure_library.generator import CustomGenerator, as_custom_generator


__all__ = []


class PdbWorkflowService(BaseWorkflowService, IPdbWorkflowService):

    @property
    def name(self):
        return "PdbWorkflowService"

    @lock_registry(kwargs_name="parameter_id")
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def pdb_to_components(
            self,
            pdb_path: str,
            edge_data_file: Optional[str] = None,
            verbose: bool = False,
            parameter_id: str = ''
    ) -> CustomGenerator[Tuple[str, List[PebbleGameResult]], pd.DataFrame]:

        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id)
        pdb_service = StructureServiceRegistry.PDB.single_service()
        pdb_io_service = IOServiceRegistry.PDB.single_service()

        ################
        # Read PDB #####
        ################

        raw_df = pdb_io_service.read(pdb_path).export_first_model()

        if parameter_registry.get_parameter(HydrogenBondParameter.INCLUDE.value) and not pdb_service.has_hydrogen_atoms(raw_df, parameter_id=parameter_id):
            raise ValueError(f"No hydrogen atoms in ATOM records of {pdb_path}")

        ###################
        # Build Graph #####
        ###################

        atom_df = pdb_service.export_atom_df(raw_df=raw_df, check_ids=True, parameter_id=parameter_id)
        others_df = pdb_service.export_others_df(raw_df=raw_df, ter_only=False, parameter_id=parameter_id)

        protein_graph: ProteinGraph = pdb_service.create_graph_struct(atom_df=atom_df, others_df=others_df,
                                                        parameter_id=parameter_id)

        # insert interactions into first frame
        pdb_service.apply_non_covalent_interactions(protein_graph=protein_graph, parameter_id=parameter_id)

        if edge_data_file is not None:
            bond_frame = export_bond_frame(protein_graph.graphs['full'], include_base_edges=True)
            IOServiceRegistry.BND.single_service().store_bonds(edge_data_file, bond_frame, mode='w')

        @as_custom_generator(Tuple[str, List[PebbleGameResult]], pd.DataFrame)
        def generator():
            yield from apply_pebble_game(protein_graph, verbose=verbose)

            return protein_graph.hydrogen_mapping

        return generator()  # noqa - decorator converts function to correct type


WorkflowServiceRegistry.PDB.register_service(PdbWorkflowService())
