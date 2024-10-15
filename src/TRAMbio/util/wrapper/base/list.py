from typing import Iterable, Optional, Any


def first(iterable: Optional[Iterable]) -> Optional[Any]:
    """
    Get the first value from an iterable (if exists)
    """
    return next(iter(iterable or []), None)
