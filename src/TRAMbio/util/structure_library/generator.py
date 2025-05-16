from __future__ import annotations
import functools
from typing import Iterable, Generic, TypeVar, Generator, Union, Type


__all__ = ["CustomGenerator", "as_custom_generator"]


_I = TypeVar('_I')
"""Generic Iterator Type"""
_R = TypeVar('_R')
"""Generic Return Type"""


class CustomGenerator(Iterable, Generic[_I, _R]):
    """Generator wrapper with access to return value"""
    def __init__(self, gen: CustomGenerator[_I, _R] | Generator[_I, None, _R]):
        self._gen = gen
        self._value = None

    def __iter__(self) -> _I:
        self._value = yield from self._gen
        return self._value

    def stop(self) -> _I | _R | None:
        return self._value


def as_custom_generator(iter_type: type, return_type: type):

    def as_custom_generator_decorator(f):
        @functools.wraps(f)
        def wrapped_func(*args, **kwargs):
            generator = f(*args, **kwargs)
            if isinstance(generator, CustomGenerator) or isinstance(generator, Generator):
                return CustomGenerator[iter_type, return_type](generator)
            raise ValueError(f"Function {f.__name__} is not a generator.")
        return wrapped_func
    return as_custom_generator_decorator
