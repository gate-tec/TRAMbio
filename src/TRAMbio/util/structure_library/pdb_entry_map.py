from typing import TypedDict, List, Type, Union


class EntryMap(TypedDict):
    id: str
    line: List[int]
    type: Union[Type[str], Type[int], Type[float]]
    required: bool
