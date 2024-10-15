import textwrap
from typing import Optional

from TRAMbio.services.io import IPyMolIOService, IOServiceRegistry

__all__ = []


class PyMolIOService(IPyMolIOService):

    @property
    def name(self):
        return "PyMolIOService"

    def write_pymol_template(
            self,
            pml_path: str,
            out_prefix: str,
            pdb_path: str,
            num_states: Optional[int],
            max_color_value: int,
            bond_commands: Optional[str]
    ) -> None:
        with open(pml_path, 'w') as pml_f:
            pml_f.write(
                textwrap.dedent("""
                ###########################LENGTH#######
                ###### PyMol Commands for #PREFIX# #####
                ###########################LENGTH#######
                reinitialize
                load #FILE#, #PREFIX#, discrete=1
                hide all
                """)
                .strip()
                .replace('#PREFIX#', out_prefix)
                .replace('#FILE#', pdb_path)
                .replace('#LENGTH#', '#' * len(out_prefix))
                + '\n'
                )

            if num_states is not None:
                # Structural alignment for trajectories
                pml_f.write(
                    textwrap.dedent("""
                    load #FILE#, #PREFIX#_ref, discrete=0
                    intra_fit #PREFIX#_ref
                    """)
                    .replace('#PREFIX#', out_prefix)
                    .replace('#FILE#', pdb_path)
                    + '\n'
                )
                for state in range(1, num_states + 1):
                    pml_f.write(
                        "stored.coords = {}\n"
                        f"iterate_state {state}, {out_prefix}_ref, stored.coords[segi+'/'+chain+'/'+resn+'`'+resi+'/'+name] = (x,y,z)\n"
                        f"alter_state {state}, {out_prefix}, (x,y,z) = stored.coords.get(segi+'/'+chain+'/'+resn+'`'+resi+ '/'+name)\n"
                    )
                pml_f.write(f'dele {out_prefix}_ref\n')

            pml_f.write(f"\nspectrum b, rainbow, {out_prefix}, minimum=0.0, maximum={max_color_value:.2f}\n\n")

            pml_f.write(
                textwrap.dedent("""
                set_bond stick_color, white, #PREFIX#
                set_bond stick_radius, 0.15, #PREFIX#
                show sticks, #PREFIX#

                set sphere_scale, 0.5, #PREFIX#
                set sphere_transparency, 0.5
                show spheres, #PREFIX# & b>0
                orient
                """)
                .strip()
                .replace('#PREFIX#', out_prefix)
                + '\n'
            )

            if bond_commands is not None:
                pml_f.write(f'\n\n{bond_commands}\n')


IOServiceRegistry.PYMOL.register_service(PyMolIOService())
