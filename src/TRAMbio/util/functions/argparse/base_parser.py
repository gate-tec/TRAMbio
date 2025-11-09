from typing import Dict, List
import argparse
import sys
import textwrap

from loguru import logger

from TRAMbio.services.parameter import BaseParameter
from TRAMbio.util.functions.argparse.config import handle_config_file
from TRAMbio.util.structure_library.argparse import OptionsDictionary, RegistryParameterAction


def parse_args_for(prog: str, description: str, usage: str, cli_options: Dict[str, OptionsDictionary], parameters: List[BaseParameter] = None):
    parser = argparse.ArgumentParser(prog=prog,
                                     usage="%(prog)s " + usage,
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    for name, option in cli_options.items():
        parser.add_argument(*option['id'], **option['args'], dest=name)

    parser.add_argument('--config', type=str, required=False, metavar='CONFIG_FILE', help=textwrap.dedent(
        """Set arguments from .ini style config file.
        """))

    if parameters is not None and len(parameters) > 0:
        env_var_group = parser.add_argument_group('Env-Vars', 'Relevant environment variables')
        for parameter in parameters:
            env_var_group.add_argument(
                parameter.value,
                action=RegistryParameterAction,
                parameter=parameter,
                help=parameter.description)

    args = vars(parser.parse_args())
    if args['config'] is not None:
        handle_config_file(args['config'], args, cli_options)

    for key, value in cli_options.items():
        if args[key] is None:
            if 'default' in value.keys() and value['default'] is not None:
                args[key] = value['default'](args)
            else:
                sys.exit(f"Option {value['id'][-1]} is required.")

    return args
