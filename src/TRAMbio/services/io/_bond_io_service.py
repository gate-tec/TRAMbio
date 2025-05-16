from typing import Tuple, Generator, Set, Literal
import pandas as pd
import numpy as np

from TRAMbio.services.io import IBondIOService, IOServiceRegistry
from TRAMbio.util.constants.interaction import InteractionType


__all__ = []


def _convert_node_id_to_pymol_sel(node_id: str):
    return node_id[0] + '/' + node_id[6:9].strip() + '`' + str(int(node_id[1:5])) + '/' + node_id[10:]


class BondIOService(IBondIOService):

    @property
    def name(self):
        return "BondIOService"

    def read(self, bond_path: str) -> pd.DataFrame:
        with open(bond_path) as bond_file:
            bond_frame = pd.read_csv(bond_file, sep='\t', header=0, index_col=None)
        return bond_frame

    def store_bonds(self, bond_path: str, bond_data: pd.DataFrame, mode: Literal['w', 'a'] = 'w') -> None:
        with open(bond_path, mode) as bnd_f:
            bond_data.to_csv(bnd_f, sep='\t', header=(mode == 'w'), index=False)

    def get_bonds_for_key(
            self,
            bond_path: str,
            all_weighted_bonds: bool = False
    ) -> Generator[Set[Tuple[str, str]], str, None]:
        if all_weighted_bonds:
            bond_types = [InteractionType.H_BOND.value, InteractionType.SALT_BRIDGE.value, InteractionType.HYDROPHOBIC.value]
        else:
            bond_types = [InteractionType.H_BOND.value, InteractionType.SALT_BRIDGE.value]

        bond_frame = self.read(bond_path=bond_path)

        bond_frame: pd.DataFrame = bond_frame.loc[bond_frame['type'].isin(bond_types), :].reset_index(drop=True)
        if '_merge' in bond_frame.columns:
            # Trajectory

            def _generator() -> Generator[Set[Tuple[str, str]], str, None]:
                bond_set: Set[Tuple[str, str]] = set()
                iterator = bond_frame.iterrows()

                prev_row = None
                while True:
                    next_key = (yield bond_set)
                    if not next_key:
                        continue
                    try:
                        next_key = int(next_key)
                    except ValueError:
                        continue

                    state = 0
                    while state < 2:
                        if prev_row is not None:
                            row = prev_row
                            prev_row = None
                        else:
                            try:
                                row = next(iter(iterator))[1]
                            except StopIteration:
                                yield bond_set
                                return

                        if state == 1 and row['frame'] != next_key:
                            prev_row = row
                            state = 2
                            continue
                        if state == 0 and row['frame'] == next_key:
                            state = 1

                        atom1, atom2 = row['node1'], row['node2']

                        bond = (_convert_node_id_to_pymol_sel(atom1), _convert_node_id_to_pymol_sel(atom2))

                        if row['_merge'] == 'right_only':
                            bond_set.add(bond)
                        else:
                            bond_set.remove(bond)

        else:
            # single frame PDB
            def _generator() -> Generator[Set[Tuple[str, str]], str, None]:
                bond_set: Set[Tuple[str, str]] = set()
                while True:
                    next_key = (yield bond_set)
                    if next_key and next_key != '-INF':
                        next_key = float(next_key) + 0.001  # epsilon for rounding in XML
                        for i, row in bond_frame.iterrows():
                            if row['type'] not in bond_types:
                                # sanity check; should never reach
                                continue

                            if not np.isnan(row['key']) and row['key'] <= next_key:
                                atom1, atom2 = row['node1'], row['node2']

                                bond_set.add((_convert_node_id_to_pymol_sel(atom1), _convert_node_id_to_pymol_sel(atom2)))

        gen = _generator()
        next(gen)  # advance to starting position
        return gen


IOServiceRegistry.BND.register_service(BondIOService())
