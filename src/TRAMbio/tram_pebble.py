#!/usr/bin/env python
import os
import multiprocessing as mp
import sys
import textwrap
import time
from typing import Dict, Literal

from loguru import logger
from pathvalidate import sanitize_filename

from TRAMbio import set_log_level

from TRAMbio.services import WorkflowServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import PebbleGameParameter, GeneralWorkflowParameter
from TRAMbio.util.functions.argparse.base_parser import parse_args_for

from TRAMbio.util.structure_library.argparse import OptionsDictionary


_CLI_OPTIONS: Dict[str, OptionsDictionary] = {
    'graph': OptionsDictionary(
        id=['-g', '--graph'], args=dict(type=str, metavar='GRAPH_FILE', help=textwrap.dedent(
            """Input graph in GRAPHML format.
            """)), default=None),
    'out_dir': OptionsDictionary(
        id=['-o', '--out-dir'], args=dict(type=str, metavar='OUTPUT_DIR', help=textwrap.dedent(
            """Directory for output files. (default: next to input file)
            """)),
        default=lambda argv: os.path.dirname(os.path.abspath(argv['graph']))),
    'name': OptionsDictionary(
        id=['-n', '--name'], args=dict(type=str, metavar='GRAPH_NAME', help=textwrap.dedent(
            """Alternate name for graph (used as output prefix). If not specified, name is derived from input file name.
            """)),
        default=lambda argv: os.path.basename(argv['graph'])
        if '.' not in os.path.basename(argv['graph']) else
        os.path.basename(argv['graph'])[:os.path.basename(argv['graph']).rindex('.')]),
    'k': OptionsDictionary(
        id=['-k', '--k-param'], args=dict(type=int, metavar='INTEGER', help=textwrap.dedent(
            """Parameter k for Pebble Game. (default: 2)
            """)),
        default=lambda argv: 2),
    'l': OptionsDictionary(
        id=['-l', '--l-param'], args=dict(type=int, metavar='INTEGER', help=textwrap.dedent(
            """Parameter l for Pebble Game. (default: 3)
            """)),
        default=lambda argv: 3),
    'random': OptionsDictionary(
        id=['-r', '--random'], args=dict(action='store_true', help=textwrap.dedent(
            """Specify for the graphs edges to be tested in random order.
            """)),
        default=lambda argv: False),
    'cores': OptionsDictionary(
        id=['-c', '--cores'], args=dict(type=int, metavar='THRESHOLD', help=textwrap.dedent(
            """Number of CPU cores for multiprocessing. (default: 1, no multiprocessing)
            """)),
        default=lambda argv: 1),
    'log_level': OptionsDictionary(
        id=['--log-level'], args=dict(type=str, choices=[
            'TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING', 'ERROR', 'CRITICAL', 'NONE'
        ], help=textwrap.dedent(
            """Set the logging level during execution. (default: INFO)
            """)),
        default=lambda argv: 'INFO')
}


def run_pipeline(
        input_path: str,
        out_dir: str,
        out_prefix: str,
        param_k: int,
        param_l: int,
        threads: int,
        verbose: bool = False,
        randomized_edges: bool = False
):
    parameter_id = f"TRAM_PEBBLE_{time.perf_counter()}"

    out_file = os.path.join(out_dir, f"{out_prefix}_results.xml")

    pebble_game_workflow_service = WorkflowServiceRegistry.PEBBLE.single_service()
    parameter_registry = ParameterRegistry.get_parameter_set(parameter_id=parameter_id)
    parameter_registry.set_parameter(PebbleGameParameter.K.value, param_k)
    parameter_registry.set_parameter(PebbleGameParameter.L.value, param_l)
    parameter_registry.set_parameter(PebbleGameParameter.THREADS.value, threads)
    parameter_registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, verbose)

    pebble_game_workflow_service.analyze_graph(
        graph_ml_path=input_path,
        out_path=out_file,
        parameter_id=parameter_id
    )


def main(default_log_level: Literal['TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING', 'ERROR', 'CRITICAL', 'NONE'] = "INFO"):
    set_log_level(default_log_level)

    try:
        mp.set_start_method('spawn', force=True)
    except RuntimeError:
        pass

    args = parse_args_for(
        'tram-pebble', "Run (k,l)-Pebble Game on arbitrary graphs.", _CLI_OPTIONS
    )

    ###################
    # Check input #####
    ###################

    input_file = args['graph']
    if not os.path.exists(input_file) or not str(input_file).endswith('.graphml'):
        sys.exit("Input file needs to be a valid .graphml file.")

    ####################
    # Check output #####
    ####################

    out_dir = args['out_dir']
    if os.path.isfile(out_dir):
        sys.exit(f'Given output directory "{out_dir}" is an existing filename.')
    out_prefix = sanitize_filename(args['name']).replace(' ', '_')

    ####################
    # Run pipeline #####
    ####################

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    set_log_level(args['log_level'])
    verbose = args['log_level'] in ['TRACE', 'DEBUG', 'INFO']  # verbosity argument for non-loguru messages

    try:
        run_pipeline(input_path=input_file,
                     out_dir=out_dir,
                     out_prefix=out_prefix,
                     param_k=args['k'],
                     param_l=args['l'],
                     threads=args['cores'],
                     verbose=verbose,
                     randomized_edges=args['random']
                     )
        logger.success('Done')
    except Exception as e:
        if args['log_level'] in ['TRACE']:
            raise e
        else:
            logger.critical(str(e))
