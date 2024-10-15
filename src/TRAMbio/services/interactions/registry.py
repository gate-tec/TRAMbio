from typing import List
import abc

from TRAMbio.services.core import IBaseService, BaseServiceRegistry
from TRAMbio.util.constants.interaction import InteractionType
from TRAMbio.util.structure_library.graph_struct import ProteinGraph
from TRAMbio.services.parameter import ParameterRegistry, HydrogenBondParameter, HydrophobicInteractionParameter, \
    DisulphideBridgeParameter, CationPiInteractionParameter, AromaticInteractionParameter, PdbEntryInteractionParameter


__all__ = ["InteractionServiceRegistry", "IInteractionService", "InteractionServiceException"]


ParameterRegistry.register_parameter(HydrogenBondParameter.INCLUDE.value, True)
ParameterRegistry.register_parameter(HydrogenBondParameter.MINIMAL_LENGTH.value, 2.6, lambda x: x >= 0.0)
ParameterRegistry.register_parameter(HydrogenBondParameter.ENERGY_THRESHOLD.value, -0.1)
ParameterRegistry.register_parameter(HydrogenBondParameter.CUTOFF_DISTANCE.value, 3.0, lambda x: x > 0.0)
ParameterRegistry.register_parameter(HydrogenBondParameter.STRONG_ENERGY_THRESHOLD.value, 0.0)
ParameterRegistry.register_parameter(HydrogenBondParameter.BAR_COUNT.value, 5, lambda x: 1 <= x <= 6)

ParameterRegistry.register_parameter(HydrophobicInteractionParameter.INCLUDE.value, True)
ParameterRegistry.register_parameter(HydrophobicInteractionParameter.MINIMAL_LENGTH.value, True)
ParameterRegistry.register_parameter(HydrophobicInteractionParameter.POTENTIAL.value, False)
ParameterRegistry.register_parameter(HydrophobicInteractionParameter.SURFACE_CUTOFF_DISTANCE.value, 0.25, lambda x: x > 0.0)
ParameterRegistry.register_parameter(HydrophobicInteractionParameter.POTENTIAL_CUTOFF_DISTANCE.value, 9.0, lambda x: x > 0.0)
ParameterRegistry.register_parameter(HydrophobicInteractionParameter.SCALE_14.value, 0.5, lambda x: x > 0.0)
ParameterRegistry.register_parameter(HydrophobicInteractionParameter.SCALE_15.value, 1.0, lambda x: x > 0.0)
ParameterRegistry.register_parameter(HydrophobicInteractionParameter.SCALE_UNBOUNDED.value, 1.0, lambda x: x > 0.0)
ParameterRegistry.register_parameter(HydrophobicInteractionParameter.ENERGY_THRESHOLD.value, -0.1)
ParameterRegistry.register_parameter(HydrophobicInteractionParameter.BAR_COUNT.value, 3, lambda x: 1 <= x <= 5)

ParameterRegistry.register_parameter(DisulphideBridgeParameter.INCLUDE.value, True)
ParameterRegistry.register_parameter(DisulphideBridgeParameter.CUTOFF_DISTANCE.value, 3.0, lambda x: x > 0.0)

ParameterRegistry.register_parameter(CationPiInteractionParameter.INCLUDE.value, True)
ParameterRegistry.register_parameter(CationPiInteractionParameter.CUTOFF_DISTANCE.value, 6.0, lambda x: x > 0.0)
ParameterRegistry.register_parameter(CationPiInteractionParameter.BAR_COUNT.value, 3, lambda x: 1 <= x <= 5)

ParameterRegistry.register_parameter(AromaticInteractionParameter.INCLUDE.value, True)
ParameterRegistry.register_parameter(AromaticInteractionParameter.CUTOFF_DISTANCE_PI.value, 7.0, lambda x: x > 0.0)
ParameterRegistry.register_parameter(AromaticInteractionParameter.CUTOFF_DISTANCE_T.value, 5.0, lambda x: x > 0.0)
ParameterRegistry.register_parameter(AromaticInteractionParameter.ANGLE_VARIANCE.value, 5.0, lambda x: x > 0.0)
ParameterRegistry.register_parameter(AromaticInteractionParameter.BAR_COUNT.value, 3, lambda x: 1 <= x <= 5)

ParameterRegistry.register_parameter(PdbEntryInteractionParameter.SSBOND_INCLUDE.value, True)
ParameterRegistry.register_parameter(PdbEntryInteractionParameter.LINK_INCLUDE.value, True)
ParameterRegistry.register_parameter(PdbEntryInteractionParameter.CONECT_INCLUDE.value, True)


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
