from __future__ import annotations
import os
from typing import Optional, Dict, TypeVar, Tuple, Callable

import functools
import inspect
import time
from loguru import logger

import random

from TRAMbio.services.parameter import GeneralWorkflowParameter

__all__ = ["ParameterRegistry", "DefaultParameterRegistry", "lock_registry", "verbosity_from_parameter"]

_Parameter = TypeVar('_Parameter', str, int, float, bool)


def _bool_env_var_parser(x: str) -> bool:
    return x.lower() in ('true', '1', 't')


class ParameterRegistry:
    __default_values: dict[str, tuple[_Parameter, type, Callable[[_Parameter], bool]]] = {}
    __instances = {}

    def __init__(self, parameter_id: str):
        self.__parameter_id = parameter_id
        self.__instance_values: dict[str, _Parameter] = {}
        self.__acquired_lock: str | None = None

    @classmethod
    def register_parameter(
            cls,
            parameter_key: str,
            default_value: _Parameter | None,
            validator: Callable[[_Parameter], bool] | None = None
    ):
        param_type = type(default_value) if default_value is not None else str
        env_var_parser = param_type
        if param_type == bool:
            env_var_parser = _bool_env_var_parser
        env_var = os.getenv(parameter_key)
        if env_var is not None:
            try:
                value = env_var_parser(env_var)
                if validator is not None and not validator(value):
                    logger.warning(f"Unable to set {parameter_key} to {value} using {default_value} instead.")
                else:
                    default_value = value
            except ValueError:
                pass
        if validator is None:
            cls.__default_values[parameter_key] = (default_value, param_type, lambda x: True)
        else:
            cls.__default_values[parameter_key] = (default_value, param_type, validator)

    @classmethod
    def get_parameter_set(cls, parameter_id: str) -> ParameterRegistry:
        if parameter_id not in cls.__instances.keys():
            cls.__instances[parameter_id] = cls(parameter_id=parameter_id)
        return cls.__instances[parameter_id]

    @classmethod
    def load_registry(cls, registry: ParameterRegistry):
        cls.__instances[registry.__parameter_id] = registry

    def __call__(self, *args, **kwargs) -> _Parameter | None:
        if len(args) == 1 and isinstance(args[0], str):
            parameter_key = str(args[0])
            if parameter_key in self.__instance_values.keys():
                return self.__instance_values[parameter_key]
            if parameter_key in self.__default_values.keys():
                return self.__default_values[parameter_key][0]
        return None

    def get_parameter(self, parameter_key: str) -> _Parameter | None:
        return self.__call__(parameter_key)

    def set_parameter(self, parameter_key: str, value: _Parameter):
        if self.__acquired_lock:
            logger.warning(f'Attempt to set "{parameter_key}" in parameter set "{self.__parameter_id}" during lock.')
            return
        if parameter_key not in self.__default_values.keys():
            raise KeyError(f"Unknown parameter: {parameter_key}")
        reference = self.__default_values[parameter_key]
        typing: type = reference[1]
        # PyCharm may flag a potential misuse of the typing library (ignored with noqa)
        if not isinstance(value, typing) and not (typing == float and isinstance(value, int)):  # noqa
            raise ValueError(f'Mismatching type for parameter "{parameter_key}". '
                             f'Expected {typing}, got {type(value)} instead.')
        if not reference[2](value):
            raise ValueError(f'Unable to pass check for parameter "{parameter_key}".')
        self.__instance_values[parameter_key] = value

    def clone(self, new_parameter_key: str):
        cloned_registry = self.__class__.get_parameter_set(new_parameter_key)
        cloned_registry.__instance_values = self.__instance_values.copy()
        return cloned_registry

    def acquire_lock(self, call_id: str):
        if self.__acquired_lock is None:
            self.__acquired_lock = call_id

    def release_lock(self, call_id: str):
        if self.__acquired_lock is not None and self.__acquired_lock == call_id:
            self.__acquired_lock = None

    def __str__(self):
        return f"ParameterRegistry(id={self.__parameter_id})"


DefaultParameterRegistry: ParameterRegistry = ParameterRegistry.get_parameter_set('')


def __value_by_signature(f, kwargs_name, args, kwargs, default = None, signature_default: bool = False):
    # passed as keyword
    if kwargs_name in kwargs.keys():
        return kwargs[kwargs_name]
    # passed as argument (get position from signature)
    signature = inspect.signature(f)
    if kwargs_name in signature.parameters.keys():
        args_pos = list(signature.parameters.keys()).index(kwargs_name)
        if len(args) > args_pos:
            return args[args_pos]
        elif signature_default and signature.parameters[kwargs_name].default is not inspect.Parameter.empty:
            return signature.parameters[kwargs_name].default
    return default

def lock_registry(kwargs_name: str = 'parameter_id'):
    def lock_registry_decorator(f):
        @functools.wraps(f)
        def wrapped_func(*args, **kwargs):
            param_id = __value_by_signature(f, kwargs_name, args, kwargs, signature_default=True)

            if param_id is None:
                logger.warning(f'Unable to extract argument {kwargs_name} from function {f.__name__}')
                return f(*args, **kwargs)
            if not isinstance(param_id, str):
                logger.warning(f'Argument {kwargs_name} from function {f.__name__} is incompatible to parameter_id. '
                               f'Expected <class \'str\'>, got {type(param_id)} instead.')
                return f(*args, **kwargs)

            call_id = f"{random.randint(0, 100)}{time.monotonic_ns()}"
            registry: ParameterRegistry = ParameterRegistry.get_parameter_set(parameter_id=param_id)

            registry.acquire_lock(call_id=call_id)
            return_value = f(*args, **kwargs)
            registry.release_lock(call_id=call_id)

            return return_value

        return wrapped_func

    return lock_registry_decorator


def verbosity_from_parameter(parameter_name: str = "parameter_id", verbose_name: str = "verbose"):
    def verbosity_from_parameter_decorator(f):
        @functools.wraps(f)
        def wrapped_func(*args, **kwargs):
            verbosity = __value_by_signature(f, verbose_name, args, kwargs)

            if verbosity is None:
                # if no parameter for verbosity is passed, use registry value
                param_id = __value_by_signature(f, parameter_name, args, kwargs, signature_default=True)

                if param_id is None:
                    logger.warning(f'Unable to extract argument {parameter_name} from function {f.__name__}')
                elif not isinstance(param_id, str):
                    logger.warning(f'Argument {parameter_name} from function {f.__name__} is incompatible to parameter_id. '
                                   f'Expected <class \'str\'>, got {type(param_id)} instead.')
                else:
                    registry: ParameterRegistry = ParameterRegistry.get_parameter_set(parameter_id=param_id)
                    verbosity = registry(GeneralWorkflowParameter.VERBOSE.value)

            if verbosity is None:
                logger.warning(f'Unable to extract argument {verbose_name} from function {f.__name__}, defaulting to False.')
                verbosity = False

            try:
                return_value = f(*args, **dict(kwargs, **{verbose_name: verbosity}))
            except TypeError:
                logger.warning(f'Incompatible argument {verbose_name} for function {f.__name__}')
                return_value = f(*args, **kwargs)

            return return_value

        return wrapped_func

    return verbosity_from_parameter_decorator
