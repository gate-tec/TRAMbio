import operator
from typing import Optional

import pytest

from TRAMbio.services import lock_registry
from TRAMbio.services.interactions.registry import IInteractionService
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import ProteinGraph
from tests.mock import MockServiceEvalException


class MockEmptyHydrophobicPotentialInteractionService(IInteractionService):
    """
    Parameter-independent service ensuring no interactions.
    """

    def __init__(self, tested_edge_count: Optional[int] = None):
        self.tested_edge_count = tested_edge_count

    @property
    def name(self):
        return "HydrophobicPotentialInteractionService"

    @property
    def interaction_types(self):
        return [InteractionType.HYDROPHOBIC]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        __tracebackhide__ = operator.methodcaller("errisinstance", MockServiceEvalException)
        if self.tested_edge_count and len(protein_graph.graphs['atom'].edges) != self.tested_edge_count:
            raise MockServiceEvalException(f"wrong number of edges detected. Expected {self.tested_edge_count}, got {len(protein_graph.graphs['atom'].edges)} instead.")
        pass