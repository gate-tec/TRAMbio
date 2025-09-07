from __future__ import annotations
from contextlib import AbstractContextManager, contextmanager
from typing import Tuple, Union, Generator, TextIO, List, Dict, Set, Optional, Literal
import abc

import xml.etree.ElementTree as ET
from io import StringIO
import pandas as pd
import networkx as nx

from TRAMbio.services.core import IBaseService, BaseServiceRegistry, ExtendedServiceRegistry
from TRAMbio.util.structure_library.components import StructureRef
from TRAMbio.util.errors import MissingDependencyError
try:
    from TRAMbio.util.wrapper.biopandas.pandas_pdb import CustomPandasPdb
except MissingDependencyError:
    # ignore dependency errors in interface definition
    CustomPandasPdb = None
    pass


__all__ = [
    "IOServiceRegistry",
    "IPdbIOService", "IXtcIOService", "IXmlIOService", "IBondIOService", "IPyMolIOService",
    "AbstractPdbIOContext"
]


class AbstractPdbIOContext(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def write_model(self, model: CustomPandasPdb, model_idx: int) -> None:
        raise NotImplementedError


class IPdbIOService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'read') and
                callable(subclass.read) and
                hasattr(subclass, 'pdb_file_context') and
                callable(subclass.pdb_file_context)):
            return True
        return NotImplemented

    @abc.abstractmethod
    def read(self, input_data: str | StringIO, verbose: bool = True) -> CustomPandasPdb:
        raise NotImplementedError

    @contextmanager
    @abc.abstractmethod
    def pdb_file_context(
            self,
            pdb_path: str,
            header_stream: StringIO
    ) -> Generator[AbstractPdbIOContext]:
        raise NotImplementedError


class IXtcIOService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'read') and
                callable(subclass.read)):
            return True
        return NotImplemented

    @abc.abstractmethod
    def read(self, xtc_path: str, pdb_path: str, stride: int) -> tuple[int, Generator[tuple[int, pd.DataFrame]]]:
        raise NotImplementedError


class IXmlIOService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'read') and
                callable(subclass.read) and
                hasattr(subclass, 'validate_components_xml') and
                callable(subclass.validate_xml) and
                hasattr(subclass, 'write_temp_xml_fragment') and
                callable(subclass.write_temp_xml_fragment) and
                hasattr(subclass, 'convert_temp_to_xml') and
                callable(subclass.convert_temp_to_xml)):
            return True
        return NotImplemented

    @abc.abstractmethod
    def read(self, xml_path: str) -> tuple[ET.Element, ET.Element]:
        raise NotImplementedError

    @abc.abstractmethod
    def read_graphml(self, graphml_path: str) -> nx.Graph | nx.MultiGraph:
        raise NotImplementedError

    @abc.abstractmethod
    def validate_xml(self, xml_path: str) -> bool:
        raise NotImplementedError

    @abc.abstractmethod
    def state_context(
            self,
            xml_file: TextIO,
            structure_key: str
    ) -> AbstractContextManager:
        raise NotImplementedError

    @abc.abstractmethod
    def write_temp_xml_fragment(self, xml_file: TextIO, halo: list[str] | None, sub_components: list[str]) -> None:
        raise NotImplementedError

    @abc.abstractmethod
    def convert_temp_to_xml(
            self,
            xml_path: str,
            temp_path: str,
            base_components: list[str],
            num_base_components: int,
            component_mapping: dict[str, StructureRef],
            is_trajectory: bool,
            discarded_keys: list[str] | None = None,
            parameter_id: str = ''
    ) -> None:
        raise NotImplementedError

    @abc.abstractmethod
    def write_pebble_game_results(
            self,
            xml_out_path: str,
            category: str,
            components: list[list[str]],
            parameter_id: str = ''
    ) -> None:
        raise NotImplementedError


class IBondIOService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'read') and
                callable(subclass.read) and
                hasattr(subclass, 'store_bonds') and
                callable(subclass.store_bonds) and
                hasattr(subclass, 'get_bonds_for_key') and
                callable(subclass.get_bonds_for_key)):
            return True
        return NotImplemented

    @abc.abstractmethod
    def read(self, bond_path: str) -> pd.DataFrame:
        raise NotImplementedError

    @abc.abstractmethod
    def store_bonds(self, bond_path: str, bond_data: pd.DataFrame, mode: Literal['w', 'a'] = 'w') -> None:
        raise NotImplementedError

    @abc.abstractmethod
    def get_bonds_for_key(
            self,
            bond_path: str,
            all_weighted_bonds: bool = False
    ) -> Generator[set[tuple[str, str]], str, None]:
        raise NotImplementedError


class IPyMolIOService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'write_pymol_template') and
                callable(subclass.write_pymol_template)):
            return True
        return NotImplemented

    @abc.abstractmethod
    def write_pymol_template(
            self,
            pml_path: str,
            out_prefix: str,
            pdb_path: str,
            num_states: int | None,
            max_color_value: int,
            bond_commands: str
    ) -> None:
        raise NotImplementedError


class _IOServiceRegistry:

    __PDB = ExtendedServiceRegistry[IPdbIOService, MissingDependencyError]()
    __XTC = ExtendedServiceRegistry[IXtcIOService, MissingDependencyError]()
    __XML = BaseServiceRegistry[IXmlIOService]()
    __BND = BaseServiceRegistry[IBondIOService]()
    __PYMOL = BaseServiceRegistry[IPyMolIOService]()

    @property
    def PDB(self) -> ExtendedServiceRegistry[IPdbIOService, MissingDependencyError]:
        return self.__PDB

    @property
    def XTC(self) -> ExtendedServiceRegistry[IXtcIOService, MissingDependencyError]:
        return self.__XTC

    @property
    def XML(self) -> BaseServiceRegistry[IXmlIOService]:
        return self.__XML

    @property
    def BND(self) -> BaseServiceRegistry[IBondIOService]:
        return self.__BND

    @property
    def PYMOL(self) -> BaseServiceRegistry[IPyMolIOService]:
        return self.__PYMOL


IOServiceRegistry = _IOServiceRegistry()
