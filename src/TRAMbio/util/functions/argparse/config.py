import os
from typing import Dict

from loguru import logger
import configparser

from TRAMbio.util.structure_library.argparse import OptionsDictionary


def handle_config_file(config_file: str, arguments: dict, options_dict: Dict[str, OptionsDictionary]):
    if not os.path.exists(config_file):
        logger.warning(f'Config file "{config_file}" not found')

    config = configparser.ConfigParser()
    config.read(config_file)

    for section in config.sections():
        for key, value in config[section].items():
            if key in options_dict.keys() and arguments[key] is None:
                constraints = options_dict[key]['args']
                if 'choices' in constraints.keys():
                    type_func = constraints['type'] if 'type' in constraints.keys() else str
                    try:
                        choice = type_func(value)
                        if choice in constraints['choices']:
                            arguments[key] = choice
                    except ValueError:
                        pass
                elif 'type' in constraints.keys():
                    try:
                        arguments[key] = constraints['type'](value)
                    except ValueError:
                        pass
                else:
                    arguments[key] = value

    logger.debug(arguments)
