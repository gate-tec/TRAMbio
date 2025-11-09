from typing import List, Optional, Dict, TypedDict, Callable
import argparse
import os

from TRAMbio.services import ParameterRegistry
from TRAMbio.services.parameter import BaseParameter


class OptionsDictionary(TypedDict):
    id: List[str]
    args: Dict
    default: Optional[Callable]


class EnvVarAction(argparse.Action):
    def __init__(self, envvar: str, required=True, default=None, **kwargs):
        if envvar is not None and envvar in os.environ:
            default = os.environ[envvar]
        if required and default is not None:
            required = False
        super(EnvVarAction, self).__init__(default=default, required=required, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)


class RegistryParameterAction(argparse.Action):
    def __init__(self, parameter: BaseParameter, required=True, default=None, **kwargs):
        if parameter is not None:
            default = ParameterRegistry.get_parameter_set("")(parameter.value)
        if required and default is not None:
            required = False
        super(RegistryParameterAction, self).__init__(default=default, required=required, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)