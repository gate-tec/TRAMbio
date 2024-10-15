import json
from typing import Optional

from TRAMbio.services import StructureServiceRegistry, IOServiceRegistry, ParameterRegistry, lock_registry
from TRAMbio.services.workflow import WorkflowServiceRegistry, IResidueWorkflowService
from TRAMbio.services.parameter import ResidueParameter
from TRAMbio.services.workflow.util import export_base_map, export_residue_states


__all__ = []

from TRAMbio.util.constants.xml import XMLConstants


class ResidueWorkflowService(IResidueWorkflowService):

    @property
    def name(self):
        return "ResidueWorkflowService"

    @lock_registry(kwargs_name='parameter_id')
    def convert_to_residue_level(
            self,
            xml_path: str,
            pdb_path: Optional[str],
            out_file: str,
            parameter_id: str = ''
    ) -> None:
        xml_io_service = IOServiceRegistry.XML.single_service()
        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id=parameter_id)

        if not xml_io_service.validate_xml(xml_path):
            raise ValueError(f"Provided file {xml_path} does not match expected XML schema.")

        base_components, states = xml_io_service.read(xml_path)

        stride_key = states.get(XMLConstants.KEY_ATTRIBUTE_NAME.value)
        is_trajectory = stride_key is not None

        min_key = parameter_registry(ResidueParameter.MIN_KEY.value)
        max_states = parameter_registry(ResidueParameter.MAX_STATES.value)
        # clamp threshold to percentage
        threshold = min(max(parameter_registry(ResidueParameter.THRESHOLD.value), 1e-5), 1 - 1e-5)
        use_main_chain = parameter_registry(ResidueParameter.USE_MAIN_CHAIN.value)

        ############################
        # Decode min_key_value #####
        ############################
        min_key_val = None
        if min_key is not None:
            # Check validity of key
            if is_trajectory:
                # trajectory file
                try:
                    min_key_val = int(min_key)
                except ValueError:
                    raise ValueError(f"Unable to use minimal key '{min_key}' for trajectory-type XML.")
            else:
                # single frame file
                if min_key != "-INF":
                    try:
                        min_key_val = float(min_key)
                    except ValueError:
                        raise ValueError(f"Unable to use minimal key '{min_key}' for single-frame protein XML.")
        ############################

        resi_counts = None
        if pdb_path:
            pdb_io_service = IOServiceRegistry.PDB.single_service()
            pdb_structure_service = StructureServiceRegistry.PDB.single_service()

            atom_df = pdb_structure_service.export_atom_df(
                pdb_io_service.read(pdb_path, verbose=False).export_first_model(),
                parameter_id=parameter_id
            )
            resi_counts = atom_df.loc[:, ['chain_id', 'residue_number', 'node_id']].groupby(by=['chain_id', 'residue_number']).count()
        elif not use_main_chain:
            raise ValueError("Unable to calculate rigidity percentages without a provided PDB file.")

        base_map = export_base_map(base_components, resi_counts, threshold, use_main_chain)

        key_list, components = export_residue_states(
            states=states, base_map=base_map, min_key=min_key, min_key_val=min_key_val, is_trajectory=is_trajectory
        )

        if max_states > 0:
            key_list = key_list[-min(max_states, len(key_list)):]
            components = {key: components[key] for key in key_list}

        with open(out_file, 'w') as json_file:
            json.dump(components, json_file, indent=2)


WorkflowServiceRegistry.RESIDUE.register_service(ResidueWorkflowService())
