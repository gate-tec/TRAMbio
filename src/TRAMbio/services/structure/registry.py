from io import StringIO
from typing import Tuple, Union, Optional, Dict, Set, List
import abc

import xml.etree.ElementTree as ET
from TRAMbio.services.core import IBaseService, BaseServiceRegistry, ExtendedServiceRegistry
from TRAMbio.util.structure_library.generator import as_custom_generator, CustomGenerator
from TRAMbio.util.structure_library.graph_struct import ProteinGraph
from TRAMbio.util.errors.dependencies import MissingDependencyError
from TRAMbio.services.parameter import ParameterRegistry, PdbParameter

import pandas as pd


__all__ = [
    "StructureServiceRegistry",
    "IPdbStructureService", "IXmlStructureService"
]


ParameterRegistry.register_parameter(PdbParameter.UNIQUE_BONDS.value, False)
ParameterRegistry.register_parameter(PdbParameter.KEEP_HETS.value, True)


class IPdbStructureService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'export_atom_df') and
                callable(subclass.export_atom_df) and
                hasattr(subclass, 'has_hydrogen_atoms') and
                callable(subclass.has_hydrogen_atoms) and
                hasattr(subclass, 'export_others_df') and
                callable(subclass.export_others_df) and
                hasattr(subclass, 'export_header_stream') and
                callable(subclass.export_header_stream) and
                hasattr(subclass, 'create_graph_struct') and
                callable(subclass.create_graph_struct) and
                hasattr(subclass, 'copy_graph_for_frame') and
                callable(subclass.copy_graph_for_frame) and
                hasattr(subclass, 'apply_non_covalent_interactions') and
                callable(subclass.apply_non_covalent_interactions)):
            return True
        return NotImplemented

    @abc.abstractmethod
    def export_atom_df(
            self,
            raw_df: pd.DataFrame,
            check_ids: bool = False,
            parameter_id: str = ''
    ) -> pd.DataFrame:
        raise NotImplementedError

    @abc.abstractmethod
    def has_hydrogen_atoms(
            self,
            raw_or_atom_df: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
            parameter_id: str = ''
    ) -> bool:
        raise NotImplementedError

    @abc.abstractmethod
    def export_others_df(
            self,
            raw_df: Dict[str, pd.DataFrame],
            ter_only: bool = False,
            parameter_id: str = ''
    ) -> pd.DataFrame:
        raise NotImplementedError

    @abc.abstractmethod
    def export_header_stream(
            self,
            raw_df: Dict[str, pd.DataFrame],
            pdb_name: Optional[str] = None,
            parameter_id: str = ''
    ) -> StringIO:
        raise NotImplementedError

    @abc.abstractmethod
    def create_graph_struct(
            self,
            atom_df: pd.DataFrame,
            others_df: pd.DataFrame,
            parameter_id: str = ''
    ) -> ProteinGraph:
        raise NotImplementedError

    @abc.abstractmethod
    def copy_graph_for_frame(
            self,
            atom_df: pd.DataFrame,
            others_df: pd.DataFrame,
            protein_graph: ProteinGraph,
            parameter_id: str = ''
    ) -> ProteinGraph:
        raise NotImplementedError

    @abc.abstractmethod
    def apply_non_covalent_interactions(
            self,
            protein_graph: ProteinGraph,
            parameter_id: str = ''
    ) -> None:
        raise NotImplementedError


class IXmlStructureService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'resolve_base_mapping') and
                callable(subclass.resolve_base_mapping) and
                hasattr(subclass, 'create_list_from_state') and
                callable(subclass.create_list_from_state) and
                hasattr(subclass, 'consistent_color_neighbor_states') and
                callable(subclass.consistent_color_neighbor_states) and
                hasattr(subclass, 'consistent_color_components') and
                callable(subclass.consistent_color_components)):
            return True
        return NotImplemented

    @abc.abstractmethod
    def resolve_base_mapping(
            self,
            base_components: ET.Element
    ) -> Tuple[Dict[str, Set[str]], Dict[str, int]]:
        raise NotImplementedError

    @abc.abstractmethod
    def create_list_from_state(
            self,
            xml_state_element: ET.Element,
            node_map: Dict[str, Set]
    ) -> Tuple[List[Set[str]], List[Set[str]]]:
        raise NotImplementedError

    @abc.abstractmethod
    def consistent_color_neighbor_states(
            self,
            state1: List[Set[str]],
            state2: List[Set[str]],
            base_lengths: Dict[str, int],
            color_map: Dict[int, int],
            next_color: int
    ) -> Tuple[Dict[int, int], int]:
        raise NotImplementedError

    @as_custom_generator(Tuple[str, List[Set[str]], Dict[int, int]], int)
    @abc.abstractmethod
    def consistent_color_components(
            self,
            base_components: ET.Element,
            states: ET.Element
    ) -> CustomGenerator[Tuple[str, List[Set[str]], Dict[int, int]], int]:
        raise NotImplementedError


class _StructureServiceRegistry:

    __PDB = ExtendedServiceRegistry[IPdbStructureService, MissingDependencyError]()
    __XML = BaseServiceRegistry[IXmlStructureService]()

    @property
    def PDB(self) -> ExtendedServiceRegistry[IPdbStructureService, MissingDependencyError]:
        return self.__PDB

    @property
    def XML(self) -> BaseServiceRegistry[IXmlStructureService]:
        return self.__XML


StructureServiceRegistry = _StructureServiceRegistry()
