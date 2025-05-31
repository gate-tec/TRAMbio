from TRAMbio.services.parameter import HydrogenBondParameter, HydrophobicInteractionParameter, \
    DisulphideBridgeParameter, AromaticInteractionParameter, CationPiInteractionParameter, \
    PdbEntryInteractionParameter
from TRAMbio.services import ParameterRegistry


def create_registry_with_bindable_properties():
    for param_group in [
        HydrogenBondParameter, HydrophobicInteractionParameter,
        DisulphideBridgeParameter, AromaticInteractionParameter, CationPiInteractionParameter,
        PdbEntryInteractionParameter
    ]:
        for entry in param_group:
            prop = property(lambda self, k=entry.value: self.get_parameter(k))
            prop = prop.setter(lambda self, v, k=entry.value: self.set_parameter(k, v))
            setattr(ParameterRegistry, entry.value, prop)