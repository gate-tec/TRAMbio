from typing import List

import pandas as pd
import numpy as np

from TRAMbio.services.interactions.util import PotentialEvaluator
from TRAMbio.services.parameter import ParameterRegistry as ParameterRegistry, lock_registry, \
    HydrophobicInteractionParameter
from TRAMbio.services.interactions.registry import IInteractionService, InteractionServiceRegistry
from TRAMbio.services.parameter.registry import verbosity_from_parameter

from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.constants.smirnoff import DEFAULT_POTENTIAL, RESIDUE_ATOM_POTENTIAL
from TRAMbio.util.functions.graph_neighborhood import NeighborType, AtomNeighborhood
from TRAMbio.util.graphein.functions import compute_distmat
from TRAMbio.util.structure_library.graph_struct import ProteinGraph as ProteinGraph, GraphKey, add_bonds_from_frame, \
    add_missing_nodes

from loguru import logger as logger
from tqdm import tqdm


__all__ = []


class HydrophobicPotentialInteractionService(IInteractionService):

    _SCALE_MAPPING = {
        NeighborType.P_14: HydrophobicInteractionParameter.SCALE_14,
        NeighborType.P_15: HydrophobicInteractionParameter.SCALE_15,
        NeighborType.UNBOUNDED: HydrophobicInteractionParameter.SCALE_UNBOUNDED
    }

    @property
    def name(self):
        return "HydrophobicPotentialInteractionService"

    @property
    def interaction_types(self) -> List[InteractionType]:
        return [InteractionType.HYDROPHOBIC]

    @lock_registry(kwargs_name='parameter_id')
    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def apply_interactions(self, protein_graph: ProteinGraph, parameter_id: str, verbose: bool = False):
        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id)
        if not parameter_registry(HydrophobicInteractionParameter.INCLUDE.value):
            return
        if not parameter_registry(HydrophobicInteractionParameter.POTENTIAL.value):
            # this service does not use the Lennard-Jones-12-6 potential
            return
        cutoff_distance = parameter_registry(HydrophobicInteractionParameter.POTENTIAL_CUTOFF_DISTANCE.value)
        energy_threshold = parameter_registry(HydrophobicInteractionParameter.ENERGY_THRESHOLD.value)
        minimal_length = parameter_registry(HydrophobicInteractionParameter.MINIMAL_LENGTH.value)
        bar_count = parameter_registry(HydrophobicInteractionParameter.BAR_COUNT.value)
        # a hydrophobic interaction is
        # placed for every carbon–carbon, carbon–sulfur, or
        # sulfur–sulfur pair, when their van der Waals surfaces are
        # within a cutoff distance of 0.25 Å.

        pdb_df = protein_graph.atom_df.copy()
        hydrogen_mapping = protein_graph.hydrogen_mapping

        naive_vdw = pdb_df['atom_name'].map(DEFAULT_POTENTIAL).apply(pd.Series)

        vdw = (
            pd.DataFrame(RESIDUE_ATOM_POTENTIAL)
            .unstack()
            .rename_axis(("residue_name", "atom_name"))
            .apply(pd.Series)
            .drop([0], axis=1)
        )

        pdb_df = pdb_df.join(vdw, on=["residue_name", "atom_name"])

        # Fill the NaNs with the naive states
        if "epsilon" in naive_vdw.columns:
            pdb_df = pdb_df.fillna(value={"epsilon": naive_vdw['epsilon'], "rmin_half": naive_vdw['rmin_half']})

        missing_indices = pdb_df['epsilon'].index[pdb_df['epsilon'].apply(np.isnan)]

        if len(missing_indices) > 0:
            evaluator = PotentialEvaluator(protein_graph.graphs['full'], pdb_df, hydrogen_mapping)
            for index in missing_indices:
                node_id = pdb_df.loc[index, 'node_id']
                epsilon, rmin_half = evaluator.evaluate_node(node_id)
                pdb_df.loc[index, 'epsilon'] = epsilon
                pdb_df.loc[index, 'rmin_half'] = rmin_half

        pdb_df['sigma'] = pdb_df['rmin_half'].apply(lambda x: 2.0 * x / 2 ** (1 / 6))

        # compute eucl distances between valid atoms
        dist_mat = compute_distmat(pdb_df)

        # Threshold Distance Matrix to entries where the eucl distance is less than or equal to
        # the cutoff distance and larger than or equal to 2.2 (in order to exclude covalent bonds)
        dist_mat = dist_mat[dist_mat >= 2.2]
        t_distmat = dist_mat[dist_mat <= cutoff_distance]

        hydrophobic_records = []
        num_potential_interactions = 0
        inds = list(i for i in zip(*np.where(~np.isnan(t_distmat))) if i[0] < i[1])
        with AtomNeighborhood(protein_graph.graphs['full']) as neighborhood, tqdm(
                iterable=enumerate(inds), total=len(inds), desc="Test for Hydrophobics",
                disable=not verbose
        ) as progress_bar:
            for i, ind in progress_bar:
                length = t_distmat[ind[0]][ind[1]]
                node_1 = pdb_df["node_id"][ind[0]]
                node_2 = pdb_df["node_id"][ind[1]]

                # exclude 1-2 and 1-3 neighbors
                neighbor_type = neighborhood.distance_type(node_1, node_2)
                if neighbor_type not in [NeighborType.P_14, NeighborType.P_15, NeighborType.UNBOUNDED]:
                    continue

                # get scaling factor for neighbor type
                scale = parameter_registry.get_parameter(self._SCALE_MAPPING[neighbor_type].value)

                sigma_1 = pdb_df.loc[ind[0], 'sigma']
                sigma_2 = pdb_df.loc[ind[1], 'sigma']
                epsilon_1 = pdb_df.loc[ind[0], 'epsilon']
                epsilon_2 = pdb_df.loc[ind[1], 'epsilon']

                # Combine via Lorentz-Berthelot
                sigma = (sigma_1 + sigma_2) * 0.5
                epsilon = (epsilon_1 * epsilon_2) ** 0.5 * scale

                # Calculate Lennard-Jones-12-6 potential
                r = 1 / length
                energy = 4 * epsilon * ((sigma * r) ** 12 - (sigma * r) ** 6)
                if energy <= energy_threshold:
                    hydrophobic_records.append({
                        'node_1': node_1, 'node_2': node_2, 'bond_length': length, 'energy': energy,
                        'bond_type': InteractionType.HYDROPHOBIC.value
                    })
                    num_potential_interactions += 1
                    progress_bar.set_postfix({'HP': num_potential_interactions})

        hydrophobic_frame = pd.DataFrame.from_records(hydrophobic_records)

        if len(hydrophobic_frame) < 1:
            if verbose:
                logger.debug(f"Found 0 hydrophobic interactions within {cutoff_distance} Angstroms.")
            return

        if minimal_length:
            hydrophobic_frame = hydrophobic_frame.sort_values(by=['energy', 'bond_length'], ignore_index=True).reset_index(
                drop=True)

            minimal_indices = []
            node_ids = []
            for i, row in hydrophobic_frame.iterrows():
                node_1, node_2 = row["node_1"], row["node_2"]
                if node_1 not in node_ids and node_2 not in node_ids:
                    node_ids.append(node_1)
                    node_ids.append(node_2)
                    minimal_indices.append(i)

            hydrophobic_frame = hydrophobic_frame.loc[minimal_indices, :].reset_index(drop=True)

        # check for hydrogen bonds
        hydrogen_atoms = protein_graph.hydrogen_df.loc[
            (protein_graph.hydrogen_df['node_id'].isin(hydrophobic_frame['node_1'])) |
            (protein_graph.hydrogen_df['node_id'].isin(hydrophobic_frame['node_2']))
        ]
        if len(hydrogen_atoms) > 0:
            add_missing_nodes(graphs=protein_graph.graphs, atom_df=hydrogen_atoms, clear_key_lists=False)
            covalent_frame = protein_graph.hydrogen_mapping\
                .loc[protein_graph.hydrogen_mapping['h_id'].isin(hydrogen_atoms['node_id']), :]\
                .copy().reset_index(drop=True)\
                .rename(columns={'node_id': 'node_1', 'h_id': 'node_2', 'length': 'bond_length'})
            covalent_frame['bond_type'] = InteractionType.COVALENT.value

            add_bonds_from_frame(
                graphs=protein_graph.graphs, bond_frame=covalent_frame,
                bond_attributes=None,  # defaults to just length
                pebble_graph_key=GraphKey.COVALENT_EDGES.value,
                pebble_graph_weight=5,
                pebble_graph_quantified_keys=None
            )

        new_edges = add_bonds_from_frame(
            graphs=protein_graph.graphs, bond_frame=hydrophobic_frame,
            bond_attributes={"bond_length": 'bond_length', "energy": 'energy'},
            # annotate edges with length and energy metadata
            pebble_graph_key=GraphKey.QUANTIFIED_NON_COVALENT_EDGES.value,
            pebble_graph_weight=bar_count,
            pebble_graph_quantified_keys=['energy']  # adds energy value to edge tuple for later sorting
        )

        if verbose:
            if new_edges > 1:
                logger.info(f"Found {new_edges} hydrophobic interactions")
            if new_edges == 1:
                logger.info(f"Found 1 hydrophobic interaction")


InteractionServiceRegistry.NON_COV.register_service(HydrophobicPotentialInteractionService())
