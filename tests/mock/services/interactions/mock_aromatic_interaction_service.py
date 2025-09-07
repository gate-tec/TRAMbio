import operator
from typing import Optional

import pytest

import pandas as pd

from TRAMbio.services import lock_registry
from TRAMbio.services.interactions.registry import IInteractionService
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import ProteinGraph, add_bonds_from_frame, GraphKey
from tests.mock import MockServiceEvalException


class MockEmptyAromaticInteractionService(IInteractionService):
    """
    Parameter-independent service ensuring no interactions.
    """

    def __init__(self, tested_edge_count: Optional[int] = None):
        self.tested_edge_count = tested_edge_count

    @property
    def name(self):
        return "AromaticInteractionService"

    @property
    def interaction_types(self):
        return [InteractionType.PI_STACKING, InteractionType.T_STACKING]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        __tracebackhide__ = operator.methodcaller("errisinstance", MockServiceEvalException)
        if self.tested_edge_count and len(protein_graph.graphs['atom'].edges) != self.tested_edge_count:
            raise MockServiceEvalException(f"wrong number of edges detected. Expected {self.tested_edge_count}, got {len(protein_graph.graphs['atom'].edges)} instead.")
        pass


class MockAromaticInteractionService(IInteractionService):
    """
    Parameter-independent service ensuring no interactions.
    """

    def __init__(self, tested_edge_count: Optional[int] = None):
        self.tested_edge_count = tested_edge_count

    @property
    def name(self):
        return "AromaticInteractionService"

    @property
    def interaction_types(self):
        return [InteractionType.PI_STACKING, InteractionType.T_STACKING]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        __tracebackhide__ = operator.methodcaller("errisinstance", MockServiceEvalException)
        if self.tested_edge_count and len(protein_graph.graphs['atom'].edges) != self.tested_edge_count:
            raise MockServiceEvalException(f"wrong number of edges detected. Expected {self.tested_edge_count}, got {len(protein_graph.graphs['atom'].edges)} instead.")

        pdb_df = protein_graph.atom_df.copy()

        pivots = pdb_df \
            .loc[pdb_df.atom_name.isin(["CG", "CE1", "CE2"]), ["node_id", "atom_name", "element_symbol"]] \
            .pivot(index="atom_name", columns="node_id", values="element_symbol")

        interactions = []

        for _, row in pivots.iterrows():
            node_ids = [x for x in row.loc[~row.isnull()].index]
            for i in range(len(node_ids) - 1):
                for j in range(i + 1, len(node_ids)):
                    interactions.append({
                        'node_1': node_ids[i], 'node_2': node_ids[j], 'bond_length': 5.0,
                        'bond_type': InteractionType.PI_STACKING.value
                    })

        add_bonds_from_frame(
            graphs=protein_graph.graphs, bond_frame=pd.DataFrame.from_records(interactions),
            bond_attributes=None,  # defaults to just length
            pebble_graph_key=GraphKey.NON_COVALENT_EDGES.value,
            pebble_graph_weight=3,
            pebble_graph_quantified_keys=None
        )
