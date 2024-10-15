import os
from io import StringIO
from typing import Tuple, Optional, List, Generator, Set, Dict, Any
import abc

from TRAMbio.services import ParameterRegistry, IOServiceRegistry
from TRAMbio.services.core import IBaseService, BaseServiceRegistry, ExtendedServiceRegistry
from TRAMbio.util.errors import MissingDependencyError

import pandas as pd

from TRAMbio.services.parameter import XtcParameter, ResidueParameter, PebbleGameParameter, PyMolParameter, \
    GeneralWorkflowParameter
from TRAMbio.util.structure_library.components import StructureRef, PebbleGameResult, IntermediateComponents
from TRAMbio.util.structure_library.generator import CustomGenerator, as_custom_generator


__all__ = [
    "WorkflowServiceRegistry",
    "IBaseWorkflowService", "IPdbWorkflowService", "IXtcWorkflowService",
    "IResidueWorkflowService", "IPyMolWorkflowService", "IPebbleGameWorkflowService"
]

ParameterRegistry.register_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

ParameterRegistry.register_parameter(PebbleGameParameter.K.value, 2, lambda x: x >= 1)
ParameterRegistry.register_parameter(PebbleGameParameter.L.value, 3, lambda x: x >= 0)
ParameterRegistry.register_parameter(PebbleGameParameter.THREADS.value, 2, lambda x: 1 <= x <= os.cpu_count())

ParameterRegistry.register_parameter(
    XtcParameter.MODULE.value,
    'MDAnalysis'
    if 'MDAnalysis' in IOServiceRegistry.XTC.list_service_names() else
    (IOServiceRegistry.XTC.single_service(error_on_none=False) or 'MDAnalysis')
)
ParameterRegistry.register_parameter(XtcParameter.STRIDE.value, 50, lambda x: x >= 1)
ParameterRegistry.register_parameter(XtcParameter.DYNAMIC_SCALING.value, True)

ParameterRegistry.register_parameter(ResidueParameter.MIN_KEY.value, None)
ParameterRegistry.register_parameter(ResidueParameter.MAX_STATES.value, 0, lambda x: x >= 0)
ParameterRegistry.register_parameter(ResidueParameter.THRESHOLD.value, 0.8, lambda x: 0 <= x <= 1)
ParameterRegistry.register_parameter(ResidueParameter.USE_MAIN_CHAIN.value, True)

ParameterRegistry.register_parameter(PyMolParameter.ALL_WEIGHTED_BONDS.value, False)


class IBaseWorkflowService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'calculate_components_from_generator') and
                callable(subclass.calculate_components_from_generator) and
                hasattr(subclass, 'convert_component_archive_to_mapping') and
                callable(subclass.convert_component_archive_to_mapping) and
                hasattr(subclass, 'run_pipeline_on_generator') and
                callable(subclass.run_pipeline_on_generator)):
            return True
        return NotImplemented

    @abc.abstractmethod
    def calculate_components_from_generator(
            self,
            structure_generator: CustomGenerator[Tuple[str, List[PebbleGameResult]], Any],
            temp_path
    ) -> Tuple[Dict[str, IntermediateComponents], List[str]]:
        raise NotImplementedError

    @abc.abstractmethod
    def convert_component_archive_to_mapping(
            self,
            comp_dict: Dict[str, IntermediateComponents],
            hydrogen_dict
    ) -> Tuple[List[str], int, Dict[str, StructureRef]]:
        raise NotImplementedError

    @abc.abstractmethod
    def run_pipeline_on_generator(
            self,
            generator: CustomGenerator[Tuple[str, List[PebbleGameResult]], pd.DataFrame],
            out_file: str,
            temp_file: str,
            verbose: bool = True,
            parameter_id: str = ''
    ) -> None:
        raise NotImplementedError


class IPdbWorkflowService(IBaseWorkflowService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseWorkflowService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'pdb_to_components') and
                callable(subclass.pdb_to_components)):
            return True
        return NotImplemented

    @as_custom_generator(Tuple[str, List[PebbleGameResult]], pd.DataFrame)
    @abc.abstractmethod
    def pdb_to_components(
            self,
            pdb_path: str,
            log_level: str = 'INFO',
            edge_data_file: Optional[str] = None,
            parameter_id: str = ''
    ) -> CustomGenerator[Tuple[str, List[PebbleGameResult]], pd.DataFrame]:
        raise NotImplementedError


class IXtcWorkflowService(IBaseWorkflowService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseWorkflowService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'trajectory_to_components') and
                callable(subclass.trajectory_to_components)):
            return True
        return NotImplemented

    @as_custom_generator(Tuple[str, List[PebbleGameResult]], pd.DataFrame)
    @abc.abstractmethod
    def trajectory_to_components(
            self,
            xtc_path: str,
            pdb_path: str,
            edge_data_file: Optional[str] = None,
            cores: int = os.cpu_count() - 1,
            verbose: bool = False,
            parameter_id: str = ''
    ) -> CustomGenerator[Tuple[str, List[PebbleGameResult]], pd.DataFrame]:
        raise NotImplementedError


class IResidueWorkflowService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'convert_to_residue_level') and
                callable(subclass.convert_to_residue_level)):
            return True
        return NotImplemented

    @abc.abstractmethod
    def convert_to_residue_level(
            self,
            xml_path: str,
            pdb_path: Optional[str],
            out_file: str,
            parameter_id: str = ''
    ) -> None:
        raise NotImplementedError


class IPyMolWorkflowService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'load_states_for_sequence') and
                callable(subclass.load_states_for_sequence) and
                hasattr(subclass, 'export_header_stream_and_ter_frame') and
                callable(subclass.export_header_stream_and_ter_frame) and
                hasattr(subclass, 'construct_pdb_file') and
                callable(subclass.constructpdb_file) ):
            return True
        return NotImplemented

    @abc.abstractmethod
    def load_states_for_sequence(
            self,
            pdb_path: str,
            xtc_path: Optional[str],
            key_sequence: List[str]
    ) -> Generator[Tuple[str, pd.DataFrame], None, None]:
        raise NotImplementedError

    @abc.abstractmethod
    def export_header_stream_and_ter_frame(
            self,
            pdb_path: str,
            out_prefix: str,
            parameter_id: str = ''
    ) -> Tuple[StringIO, pd.DataFrame]:
        raise NotImplementedError

    @abc.abstractmethod
    def construct_pdb_file(
            self,
            outpdb_path: str,
            inpdb_path: str,
            out_prefix: str,
            color_map_generator: CustomGenerator,
            num_states: int,
            custom_frame_generator: Generator[Tuple[str, pd.DataFrame], None, None],
            bond_generator: Optional[Generator[Set[Tuple[str, str]], str, None]],
            verbose: bool = True,
            parameter_id: str = ''
    ) -> Tuple[int, Optional[str]]:
        raise NotImplementedError

    @abc.abstractmethod
    def convert_to_pymol_files(
            self,
            inpdb_path: str,
            in_xtc_path: Optional[str],
            in_xml_path: str,
            in_bond_path: Optional[str],
            outpdb_path: str,
            out_pml_path: str,
            out_prefix: str,
            rel_outpdb_path: str,
            verbose: bool = True,
            parameter_id: str = ''
    ) -> None:
        raise NotImplementedError


class IPebbleGameWorkflowService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'analyze_graph') and
                callable(subclass.analyze_graph)):
            return True
        return NotImplemented

    @abc.abstractmethod
    def analyze_graph(
            self,
            graph_ml_path: str,
            out_path: str,
            verbose: bool = False,
            parameter_id: str = ''
    ) -> None:
        raise NotImplementedError


class _WorkflowServiceRegistry:

    __PDB = BaseServiceRegistry[IPdbWorkflowService]()
    __XTC = BaseServiceRegistry[IXtcWorkflowService]()
    __RESIDUE = BaseServiceRegistry[IResidueWorkflowService]()
    __PYMOL = ExtendedServiceRegistry[IPyMolWorkflowService, MissingDependencyError]()
    __PEBBLE = BaseServiceRegistry[IPebbleGameWorkflowService]()

    @property
    def PDB(self) -> BaseServiceRegistry[IPdbWorkflowService]:
        return self.__PDB

    @property
    def XTC(self) -> BaseServiceRegistry[IXtcWorkflowService]:
        return self.__XTC

    @property
    def RESIDUE(self) -> BaseServiceRegistry[IResidueWorkflowService]:
        return self.__RESIDUE

    @property
    def PYMOL(self) -> ExtendedServiceRegistry[IPyMolWorkflowService, MissingDependencyError]:
        return self.__PYMOL

    @property
    def PEBBLE(self) -> BaseServiceRegistry[IPebbleGameWorkflowService]:
        return self.__PEBBLE


WorkflowServiceRegistry = _WorkflowServiceRegistry()
