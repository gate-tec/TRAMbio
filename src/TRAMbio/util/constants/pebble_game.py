from enum import Enum


class PebbleGameCategory(Enum):

    WELL_CONSTRAINED = "well-constrained"
    OVER_CONSTRAINED = "over-constrained"
    UNDER_CONSTRAINED = "under-constrained"
    OTHER = "other"
