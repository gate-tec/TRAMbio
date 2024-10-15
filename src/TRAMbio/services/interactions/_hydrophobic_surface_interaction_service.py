from typing import List

from TRAMbio.services.parameter import ParameterRegistry as ParameterRegistry, lock_registry, \
    HydrophobicInteractionParameter
from TRAMbio.services.interactions.registry import IInteractionService, InteractionServiceRegistry
from TRAMbio.services.interactions.util import calculate_hydrophobic_interactions
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import ProteinGraph as ProteinGraph, GraphKey, add_bonds_from_frame

from loguru import logger as logger


__all__ = []


class HydrophobicSurfaceInteractionService(IInteractionService):

    @property
    def name(self):
        return "HydrophobicSurfaceInteractionService"

    @property
    def interaction_types(self) -> List[InteractionType]:
        return [InteractionType.HYDROPHOBIC]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id)
        if not parameter_registry(HydrophobicInteractionParameter.INCLUDE.value):
            return
        if parameter_registry(HydrophobicInteractionParameter.POTENTIAL.value):
            # this service does not use the Lennard-Jones-12-6 potential
            return
        cutoff_distance = parameter_registry(HydrophobicInteractionParameter.SURFACE_CUTOFF_DISTANCE.value)
        minimal_length = parameter_registry(HydrophobicInteractionParameter.MINIMAL_LENGTH.value)
        bar_count = parameter_registry(HydrophobicInteractionParameter.BAR_COUNT.value)
        # a hydrophobic interaction is
        # placed for every carbon–carbon, carbon–sulfur, or
        # sulfur–sulfur pair, when their van der Waals surfaces are
        # within a cutoff distance of 0.25 Å.

        hydrophobic_frame = calculate_hydrophobic_interactions(
            graphs=protein_graph.graphs, heavy_atom_df=protein_graph.heavy_atom_df, cutoff_distance=cutoff_distance,
            minimal_length=minimal_length, verbose=verbose
        )

        if hydrophobic_frame is None:
            return

        new_edges = add_bonds_from_frame(
            graphs=protein_graph.graphs, bond_frame=hydrophobic_frame,
            bond_attributes=None,  # defaults to just length
            pebble_graph_key=GraphKey.NON_COVALENT_EDGES.value,
            pebble_graph_weight=bar_count,
            pebble_graph_quantified_keys=None
        )
        if verbose:
            if new_edges > 1:
                logger.info(f"Found {new_edges} hydrophobic interactions")
            if new_edges == 1:
                logger.info(f"Found 1 hydrophobic interaction")


InteractionServiceRegistry.NON_COV.register_service(HydrophobicSurfaceInteractionService())
