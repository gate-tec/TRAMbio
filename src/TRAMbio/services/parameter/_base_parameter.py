from __future__ import annotations
from typing import TypeVar, Callable
from enum import Enum

_Parameter = TypeVar('_Parameter', str, int, float, bool)

class BaseParameter(str, Enum):

    description: str
    default_value: _Parameter | None
    validator: Callable[[_Parameter], bool] | None

    def __new__(
            cls, title: str,
            description: str = "",
            default_value: _Parameter | None = None,
            validator: Callable[[_Parameter], bool] | None = None
    ) -> BaseParameter:
        obj = str.__new__(cls, title)
        obj._value_ = title

        obj.description = description
        obj.default_value = default_value
        obj.validator = validator if validator is not None else lambda x: True
        return obj