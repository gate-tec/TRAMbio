from typing import Dict
import argparse
import sys
import textwrap

from TRAMbio.util.functions.argparse.config import handle_config_file
from TRAMbio.util.structure_library.argparse import OptionsDictionary


def parse_args_for(prog: str, description: str, cli_options: Dict[str, OptionsDictionary]):
    parser = argparse.ArgumentParser(prog=prog,
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    for name, option in cli_options.items():
        parser.add_argument(*option['id'], **option['args'], dest=name)

    parser.add_argument('--config', type=str, required=False, metavar='CONFIG_FILE', help=textwrap.dedent(
        """Set arguments from .ini style config file.
        """))

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
