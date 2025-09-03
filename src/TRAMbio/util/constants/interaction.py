from enum import Enum


__all__ = ["InteractionType", "INTERACTION_RANKING", "map_interaction_ranking"]


#######################################
# Protein Graph Interaction Types #####
#######################################


class InteractionType(Enum):
    COVALENT = 'covalent'
    PEPTIDE_BOND = 'peptide_bond'
    PHOSPHODIESTER_BOND = 'phosphodiester_bond'
    SS_BOND = 'disulphide'
    SALT_BRIDGE = 'salt_bridge'
    H_BOND = 'hbond'
    PI_STACKING = 'pi_stacking'
    T_STACKING = 't_stacking'
    HYDROPHOBIC = 'hydrophobic'
    CATION_PI = 'cation_pi'
    CONECT = 'CONECT'
    LINK = 'LINK'
    UNKNOWN = 'UNKNOWN'


INTERACTION_RANKING = {
    InteractionType.COVALENT.value: 0,
    InteractionType.PEPTIDE_BOND.value: 1,
    InteractionType.PHOSPHODIESTER_BOND.value: 1,  # can have PEPTIDE_BOND value, since a double assignment is impossible
    InteractionType.SS_BOND.value: 2,
    InteractionType.SALT_BRIDGE.value: 3,
    InteractionType.H_BOND.value: 4,
    InteractionType.PI_STACKING.value: 5,
    InteractionType.T_STACKING.value: 6,
    InteractionType.HYDROPHOBIC.value: 7,
    InteractionType.CATION_PI.value: 8,
    InteractionType.CONECT.value: 9,
    InteractionType.LINK.value: 10,
    InteractionType.UNKNOWN.value: 11
}


def map_interaction_ranking(x: str):
    try:
        return INTERACTION_RANKING[x]
    finally:
        return INTERACTION_RANKING[InteractionType.UNKNOWN.value]
