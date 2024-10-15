from io import StringIO
from typing import Tuple, Optional, List, Set, Generator

import pandas as pd

from TRAMbio.services import StructureServiceRegistry, IOServiceRegistry, lock_registry, ParameterRegistry
from TRAMbio.services.io import AbstractPdbIOContext
from TRAMbio.services.parameter.registry import verbosity_from_parameter
from TRAMbio.services.workflow import WorkflowServiceRegistry, IPyMolWorkflowService
from TRAMbio.services.parameter import XtcParameter, PyMolParameter
from TRAMbio.util.constants.xml import XMLConstants
from TRAMbio.util.structure_library.generator import CustomGenerator as CustomGenerator
from TRAMbio.util.wrapper.biopandas.pandas_pdb import CustomPandasPdb

from tqdm import tqdm


__all__ = []


class PyMolWorkflowService(IPyMolWorkflowService):

    @property
    def name(self):
        return "PyMolWorkflowService"

    @lock_registry(kwargs_name='parameter_id')
    def load_states_for_sequence(
            self,
            pdb_path: str,
            xtc_path: Optional[str],
            key_sequence: List[str],
            parameter_id: str = ''
    ) -> Generator[Tuple[str, pd.DataFrame], None, None]:
        pdb_service = StructureServiceRegistry.PDB.single_service()
        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id=parameter_id)

        if xtc_path is not None:
            module = parameter_registry.get_parameter(XtcParameter.MODULE.value)
            stride = parameter_registry.get_parameter(XtcParameter.STRIDE.value)
            xtc_io_service = IOServiceRegistry.XTC.query_service(module)

            _, frame_generator = xtc_io_service.read(xtc_path=xtc_path, pdb_path=pdb_path, stride=stride)

            for frame_number, raw_df in frame_generator:
                frame_number = str(frame_number)
                if frame_number in key_sequence:
                    atom_df = pdb_service.export_atom_df(raw_df=raw_df, parameter_id=parameter_id)
                    yield frame_number, atom_df
        else:
            pdb_io_service = IOServiceRegistry.PDB.query_service('PdbIOService')
            raw_df = pdb_io_service.read(input_data=pdb_path, verbose=False).export_first_model()
            atom_df = pdb_service.export_atom_df(raw_df=raw_df, parameter_id=parameter_id)

            for state_key in key_sequence:
                yield state_key, atom_df.copy()

    @lock_registry(kwargs_name='parameter_id')
    def export_header_stream_and_ter_frame(
            self,
            pdb_path: str,
            out_prefix: str,
            parameter_id: str = ''
    ) -> Tuple[StringIO, pd.DataFrame]:
        pdb_io_service = IOServiceRegistry.PDB.single_service()
        pdb_structure_service = StructureServiceRegistry.PDB.single_service()

        model_frame = pdb_io_service.read(pdb_path, verbose=False)
        raw_df = model_frame.export_first_model()

        others_frame = pdb_structure_service.export_others_df(raw_df=raw_df, ter_only=True, parameter_id=parameter_id)
        header_stream = pdb_structure_service.export_header_stream(raw_df=raw_df, pdb_name=out_prefix,
                                                                   parameter_id=parameter_id)

        return header_stream, others_frame

    @lock_registry(kwargs_name='parameter_id')
    def construct_pdb_file(
            self,
            out_pdb_path: str,
            in_pdb_path: str,
            out_prefix: str,
            color_map_generator: CustomGenerator,
            num_states: int,
            custom_frame_generator: Generator[Tuple[str, pd.DataFrame], None, None],
            bond_generator: Optional[Generator[Set[Tuple[str, str]], str, None]],
            verbose: bool = True,
            parameter_id: str = ''
    ) -> Tuple[int, Optional[str]]:

        header_stream, ter_frame = self.export_header_stream_and_ter_frame(
            pdb_path=in_pdb_path, out_prefix=out_prefix, parameter_id=parameter_id
        )

        bond_commands_1 = []
        bond_commands_2 = set()
        b_factor_shift = 10.0

        model_idx = 1
        with IOServiceRegistry.PDB.single_service().pdb_file_context(out_pdb_path, header_stream) as pdb_out:
            pdb_out: AbstractPdbIOContext
            for key, component_nodes, color_mapping in tqdm(color_map_generator, total=num_states, desc='States',
                                                            disable=not verbose):
                state_key = ''
                while key != state_key:
                    try:
                        state_key, base_frame = next(custom_frame_generator)
                    except StopIteration:
                        exc = KeyError(f'Unable to find next state key: {key}')
                        try:
                            raise exc
                        finally:
                            exc.__context__ = None

                base_frame.loc[:, 'b_factor'] = 0.0  # default value

                for i, component in enumerate(component_nodes):
                    color = float(color_mapping[i]) / b_factor_shift  # offset b-factor value
                    base_frame.loc[base_frame['node_id'].isin(component), 'b_factor'] = color

                pdb_df = CustomPandasPdb()
                pdb_df.df['ATOM'] = base_frame.loc[base_frame['record_name'] == 'ATOM', :]
                pdb_df.df['HETATM'] = base_frame.loc[base_frame['record_name'] == 'HETATM', :]
                pdb_df.df['OTHERS'] = ter_frame
                pdb_out.write_model(pdb_df, model_idx=model_idx)

                if bond_generator is not None:
                    for (atom1, atom2) in bond_generator.send(key):
                        bond_commands_1.append(
                            f"bond state {model_idx} & /{out_prefix}/*/{atom1}, state {model_idx} & /{out_prefix}/*/{atom2}")
                        bond_commands_2.add(f"/{out_prefix}/*/{atom1}, /{out_prefix}/*/{atom2}")
                model_idx += 1

        max_color_value = (color_map_generator.stop() - 1) / b_factor_shift

        pml_text = None
        if len(bond_commands_1) > 0:
            pml_text = '\n'.join(bond_commands_1) + '\n\n' + '\n'.join(
                [f"set_bond stick_color, yellow, {entry}\nset_bond stick_radius, 0.1, {entry}" for entry in
                 bond_commands_2]
            )

        return max_color_value, pml_text

    @verbosity_from_parameter(parameter_name="parameter_id", verbose_name="verbose")
    def convert_to_pymol_files(
            self,
            in_pdb_path: str,
            in_xtc_path: Optional[str],
            in_xml_path: str,
            in_bond_path: Optional[str],
            out_pdb_path: str,
            out_pml_path: str,
            out_prefix: str,
            rel_out_pdb_path: str,
            verbose: bool = True,
            parameter_id: str = ''
    ) -> None:
        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id=parameter_id)

        xml_io_service = IOServiceRegistry.XML.single_service()
        xml_structure_service = StructureServiceRegistry.XML.single_service()
        pymol_io_service = IOServiceRegistry.PYMOL.single_service()

        if not xml_io_service.validate_xml(xml_path=in_xml_path):
            raise ValueError(f"Provided file {in_xml_path} does not match expected XML schema.")

        base_components, states = xml_io_service.read(xml_path=in_xml_path)

        is_trajectory = False
        stride = states.get(XMLConstants.KEY_ATTRIBUTE_NAME.value)
        if (stride is None) != (in_xtc_path is None):
            if in_xtc_path is None:
                raise ValueError("A single-state PDB file is incompatible with a trajectory-type XML.")
            else:
                raise ValueError("A trajectory XTC file is incompatible with a single-frame protein XML.")
        if stride is not None:
            is_trajectory = True
            parameter_registry.set_parameter(XtcParameter.STRIDE.value, int(stride))

        key_sequence = [state.get(XMLConstants.KEY_ATTRIBUTE_NAME.value) for state in states]
        num_states = len(key_sequence)

        custom_frame_generator = self.load_states_for_sequence(
            pdb_path=in_pdb_path,
            xtc_path=in_xtc_path,
            key_sequence=key_sequence,
            parameter_id=parameter_id
        )
        color_map_generator = xml_structure_service.consistent_color_components(base_components, states)

        all_weighted_bonds = parameter_registry.get_parameter(PyMolParameter.ALL_WEIGHTED_BONDS.value)

        bond_generator = None if in_bond_path is None else \
            IOServiceRegistry.BND.single_service().get_bonds_for_key(
                bond_path=in_bond_path, all_weighted_bonds=all_weighted_bonds
            )

        ############################

        max_color_value, bond_commands = self.construct_pdb_file(
            out_pdb_path=out_pdb_path,
            in_pdb_path=in_pdb_path,
            out_prefix=out_prefix,
            color_map_generator=color_map_generator,
            num_states=num_states,
            custom_frame_generator=custom_frame_generator,
            bond_generator=bond_generator,
            verbose=verbose,
            parameter_id=parameter_id
        )

        pymol_io_service.write_pymol_template(
            pml_path=out_pml_path,
            out_prefix=out_prefix,
            pdb_path=rel_out_pdb_path,
            num_states=num_states if is_trajectory else None,
            max_color_value=max_color_value,
            bond_commands=bond_commands
        )


WorkflowServiceRegistry.PYMOL.register_service(PyMolWorkflowService())
