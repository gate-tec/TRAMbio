from typing import TypedDict, List, Iterable, Optional, Set


class PebbleGameComponent(TypedDict):
    nodes: List[str]
    halo: Set[str]


class PebbleGameResult(TypedDict):
    size: int
    halo_size: int
    nodes: List[str]
    halo: List[str]


class IntermediateComponents(TypedDict):
    size: int
    nodes: Optional[List[str]]
    components: Optional[List[str]]


class StructureRef(TypedDict):
    out: List[str]
    stack: Optional[Iterable[str]]
