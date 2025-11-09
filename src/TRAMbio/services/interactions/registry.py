from typing import List
import abc

from TRAMbio.services.core import IBaseService, BaseServiceRegistry
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import ProteinGraph
from TRAMbio.services.parameter import ParameterRegistry, HydrogenBondParameter, HydrophobicInteractionParameter, \
    DisulphideBridgeParameter, CationPiInteractionParameter, AromaticInteractionParameter, PdbEntryInteractionParameter


__all__ = ["InteractionServiceRegistry", "IInteractionService", "InteractionServiceException"]


for parameter_set in [
    HydrogenBondParameter,
    HydrophobicInteractionParameter,
    DisulphideBridgeParameter,
    CationPiInteractionParameter,
    AromaticInteractionParameter,
    PdbEntryInteractionParameter
]:
    for parameter in parameter_set:
        ParameterRegistry.register_parameter(parameter)


class IInteractionService(IBaseService, metaclass=abc.ABCMeta):

    @classmethod
    def _subclasshook__(cls, subclass):
        if IBaseService.__subclasshook__(subclass) is NotImplemented:
            return NotImplemented
        if (hasattr(subclass, 'apply_interactions') and
                callable(subclass.apply_interactions) and
                hasattr(subclass, 'interaction_types') and
                callable(subclass.interaction_types)):
            return True
        return NotImplemented

    @property
    @abc.abstractmethod
    def interaction_types(self) -> List[InteractionType]:
        raise NotImplementedError

    @abc.abstractmethod
    def apply_interactions(
            self,
            protein_graph: ProteinGraph,
            parameter_id: str,
            verbose: bool = False
    ) -> None:
        raise NotImplementedError


class _InteractionServiceRegistry:

    __COV = BaseServiceRegistry[IInteractionService]()
    __NON_COV = BaseServiceRegistry[IInteractionService]()

    @property
    def COV(self) -> BaseServiceRegistry[IInteractionService]:
        return self.__COV

    @property
    def NON_COV(self) -> BaseServiceRegistry[IInteractionService]:
        return self.__NON_COV


InteractionServiceRegistry = _InteractionServiceRegistry()


class InteractionServiceException(BaseException):
    pass
