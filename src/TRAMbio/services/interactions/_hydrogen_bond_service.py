from typing import List

from TRAMbio.services.parameter import ParameterRegistry as ParameterRegistry, lock_registry, HydrogenBondParameter
from TRAMbio.services.interactions.registry import IInteractionService, InteractionServiceRegistry
from TRAMbio.services.interactions.util import calculate_hydrogen_and_salt_bridge_bonds
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import ProteinGraph as ProteinGraph, GraphKey, add_bonds_from_frame, \
    add_missing_nodes

from loguru import logger as logger


__all__ = []


class HydrogenBondService(IInteractionService):

    @property
    def name(self):
        return "HydrogenBondService"

    @property
    def interaction_types(self) -> List[InteractionType]:
        return [InteractionType.H_BOND, InteractionType.SALT_BRIDGE]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        if not ParameterRegistry.get_parameter_set(parameter_id)(HydrogenBondParameter.INCLUDE.value):
            return
        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id)
        energy_threshold = parameter_registry(HydrogenBondParameter.ENERGY_THRESHOLD.value)
        minimal_length = parameter_registry(HydrogenBondParameter.MINIMAL_LENGTH.value)

        strong_energy_threshold = parameter_registry(HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value)
        default_bar_count = parameter_registry(HydrogenBondParameter.BAR_COUNT.value)

        bond_frame, covalent_frame = calculate_hydrogen_and_salt_bridge_bonds(
            graphs=protein_graph.graphs, heavy_atom_df=protein_graph.heavy_atom_df,
            hydrogen_mapping=protein_graph.hydrogen_mapping, hydrogen_df=protein_graph.hydrogen_df,
            minimum_distance=minimal_length, energy_threshold=energy_threshold,
            verbose=verbose
        )

        if bond_frame is None:
            if verbose:
                logger.info("Found 0 hydrogen bonds")
                logger.info("Found 0 salt-bridges")
            return

        if covalent_frame is not None and len(covalent_frame) > 0:
            # hydrogen are node_2
            hydrogen_subset = protein_graph.hydrogen_df.loc[protein_graph.hydrogen_df['node_id'].isin(covalent_frame['node_2']), :] \
                .copy().reset_index(drop=True)
            # add relevant hydrogen atoms as nodes
            add_missing_nodes(graphs=protein_graph.graphs, atom_df=hydrogen_subset, clear_key_lists=False)
            add_bonds_from_frame(
                graphs=protein_graph.graphs, bond_frame=covalent_frame,
                bond_attributes=None,  # defaults to just length
                pebble_graph_key=GraphKey.COVALENT_EDGES.value,
                pebble_graph_weight=5,
                pebble_graph_quantified_keys=None
            )

        if verbose:
            logger.info(
                f"Found {len(bond_frame.loc[bond_frame['bond_type'] == InteractionType.H_BOND.value, :])} hydrogen bonds")
            logger.info(
                f"Found {len(bond_frame.loc[bond_frame['bond_type'] == InteractionType.SALT_BRIDGE.value, :])} salt-bridges")

        if len(bond_frame) == 0:
            # sanity check
            return

        if strong_energy_threshold < energy_threshold:
            strong_bond_frame = bond_frame.loc[bond_frame.energy <= strong_energy_threshold, :].copy().reset_index(drop=True)

            energy_difference = energy_threshold - strong_energy_threshold
            bar_difference = default_bar_count - 1

            if bar_difference >= 1:
                weak_bond_frame = bond_frame.loc[bond_frame.energy > strong_energy_threshold, :].reset_index(drop=True)
                interval_length = energy_difference / bar_difference
                for i in range(bar_difference):
                    interval_frame = weak_bond_frame.loc[weak_bond_frame.energy.between(
                        energy_threshold - (i + 1) * interval_length + 1e-5,
                        energy_threshold - i * interval_length
                    )].copy().reset_index(drop=True)

                    add_bonds_from_frame(
                        graphs=protein_graph.graphs, bond_frame=interval_frame,
                        bond_attributes={"bond_length": 'bond_length', "energy": 'energy', "extra": 'extra'},
                        # annotate edges with length and energy metadata
                        pebble_graph_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value,
                        pebble_graph_weight=i + 1,
                        pebble_graph_quantified_keys=['energy']  # adds energy value to edge tuple for later sorting
                    )
        else:
            strong_bond_frame = bond_frame

        add_bonds_from_frame(
            graphs=protein_graph.graphs, bond_frame=strong_bond_frame,
            bond_attributes={"bond_length": 'bond_length', "energy": 'energy', "extra": 'extra'},
            # annotate edges with length and energy metadata
            pebble_graph_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value,
            pebble_graph_weight=default_bar_count,
            pebble_graph_quantified_keys=['energy']  # adds energy value to edge tuple for later sorting
        )


InteractionServiceRegistry.NON_COV.register_service(HydrogenBondService())
