from typing import List, Optional, Dict, TypedDict, Callable


class OptionsDictionary(TypedDict):
    id: List[str]
    args: Dict
    default: Optional[Callable]
