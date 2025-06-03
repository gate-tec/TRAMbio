#!/usr/bin/env python
import os
import sys
import textwrap
import time
from typing import Dict, Literal, Optional

from TRAMbio.services import IOServiceRegistry, WorkflowServiceRegistry, DefaultParameterRegistry, ParameterRegistry

from loguru import logger
from pathvalidate import sanitize_filename

import multiprocessing as mp

from TRAMbio import set_log_level
from TRAMbio.services.parameter import XtcParameter, GeneralWorkflowParameter, HydrogenBondParameter
from TRAMbio.util.functions.argparse.base_parser import parse_args_for
from TRAMbio.util.structure_library.argparse import OptionsDictionary


_tqdm_logger = logger.bind(task="tqdm")

_CLI_OPTIONS: Dict[str, OptionsDictionary] = {
    'xtc': OptionsDictionary(
        id=['-x', '--xtc'], args=dict(type=str, metavar='XTC_FILE', help=textwrap.dedent(
            """Trajectory file in XTC format.
            """)), default=None),
    'pdb': OptionsDictionary(
        id=['-p', '--pdb'], args=dict(type=str, metavar='PDB_FILE', help=textwrap.dedent(
            """Protein input file in PDB v3.3 format.
            """)), default=None),
    'out_dir': OptionsDictionary(
        id=['-o', '--out-dir'], args=dict(type=str, metavar='OUTPUT_DIR', help=textwrap.dedent(
            """Directory for output files. (default: next to input file)
            """)),
        default=lambda argv: os.path.dirname(os.path.abspath(argv['pdb']))),
    'name': OptionsDictionary(
        id=['-n', '--name'], args=dict(type=str, metavar='PDB_NAME', help=textwrap.dedent(
            """Alternate name for protein (used as output prefix). If not specified, name is derived from input file name.
            """)),
        default=lambda argv: os.path.basename(argv['pdb'])
        if '.' not in os.path.basename(argv['pdb']) else
        os.path.basename(argv['pdb'])[:os.path.basename(argv['pdb']).rindex('.')]),
    'edges': OptionsDictionary(
        id=['-e', '--edges'], args=dict(action='store_true', help=textwrap.dedent(
            """Specify to store graph edges as BND file."""
        )),
        default=lambda argv: False),
    'cores': OptionsDictionary(
        id=['-c', '--cores'], args=dict(type=int, metavar='CORES', help=textwrap.dedent(
            """Number of CPU cores for multiprocessing. (default: maximum available cores - 1)"""
        )),
        default=lambda argv: os.cpu_count() - 1),
    'threshold': OptionsDictionary(
        id=['-t', '--threshold'], args=dict(type=float, metavar='THRESHOLD', help=textwrap.dedent(
            f"""Energy threshold for inclusion of hydrogen bonds. All bonds with energy lower or equal to this threshold are included. (default: {DefaultParameterRegistry.get_parameter(HydrogenBondParameter.ENERGY_THRESHOLD.value):.3f})
            """)),
        default=lambda argv: None),
    'stride': OptionsDictionary(
        id=['-s', '--stride'], args=dict(type=int, metavar='STRIDE', help=textwrap.dedent(
            """Only processes every stride-th frame. (default: 50)
            Negative values result in only the first frame being processed.
            """)), default=lambda argv: 50),
    'module': OptionsDictionary(
        id=['-m', '--module'], args=dict(type=str, choices=IOServiceRegistry.XTC.list_service_names(), help=textwrap.dedent(
            f"""Base module for trajectory loading. (default: {DefaultParameterRegistry.get_parameter(XtcParameter.MODULE.value)})
            """)),
        default=lambda argv: DefaultParameterRegistry.get_parameter(XtcParameter.MODULE.value)),
    'log_level': OptionsDictionary(
        id=['-l', '--log-level'], args=dict(type=str, choices=[
            'TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING', 'ERROR', 'CRITICAL', 'NONE'
        ], help=textwrap.dedent(
            """Set the logging level during execution. (default: INFO)
            """)),
        default=lambda argv: 'INFO')
}


def run_pipeline(
        input_xtc: str,
        input_pdb: str,
        out_dir: str,
        out_prefix: str,
        stride: int,
        module: str,
        energy_threshold: Optional[float],
        cores: int = os.cpu_count() - 1,
        store_edges: bool = False,
        verbose: bool = False
):
    parameter_id = f"TRAM_XTC_{time.perf_counter()}"

    logger.debug(f"stride={stride}")

    out_file = os.path.join(out_dir, f"{out_prefix}_components.xml")
    temp_file = os.path.join(out_dir, f".temp_{out_prefix}_components.xml")
    edge_data_file = os.path.join(out_dir, f"{out_prefix}_edges.bnd") if store_edges else None

    xtc_workflow_service = WorkflowServiceRegistry.XTC.single_service()
    parameter_registry = ParameterRegistry.get_parameter_set(parameter_id=parameter_id)
    parameter_registry.set_parameter(XtcParameter.STRIDE.value, stride)
    parameter_registry.set_parameter(XtcParameter.MODULE.value, module)
    parameter_registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, verbose)

    if energy_threshold:
        parameter_registry.set_parameter(HydrogenBondParameter.ENERGY_THRESHOLD.value, energy_threshold)

    # Run Pebble Game on frames
    frame_generator = xtc_workflow_service.trajectory_to_components(
        xtc_path=input_xtc,
        pdb_path=input_pdb,
        edge_data_file=edge_data_file,
        cores=cores,
        parameter_id=parameter_id
    )

    # Process results per frame
    xtc_workflow_service.run_pipeline_on_generator(
        generator=frame_generator,
        out_file=out_file,
        temp_file=temp_file,
        parameter_id=parameter_id
    )


def main(default_log_level: Literal['TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING', 'ERROR', 'CRITICAL', 'NONE'] = "INFO"):
    set_log_level(default_log_level)
    try:
        mp.set_start_method('spawn', force=True)
    except RuntimeError:
        pass

    args = parse_args_for(
        'tram-xtc', "Calculate rigid components for MD trajectory.", _CLI_OPTIONS
    )

    ###################
    # Check input #####
    ###################

    input_xtc = args['xtc']
    if not os.path.exists(input_xtc) or not str(input_xtc).endswith('.xtc'):
        sys.exit("XTC-input file needs to be a valid .xtc file.")
    input_pdb = args['pdb']
    if not os.path.exists(input_pdb) or not str(input_pdb).endswith('.pdb'):
        sys.exit("PDB-input file needs to be a valid .pdb file.")

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

    num_cores = max(min(args['cores'], os.cpu_count()), 1)

    set_log_level(args['log_level'])
    verbose = args['log_level'] in ['TRACE', 'DEBUG', 'INFO']  # verbosity argument for non-loguru messages
    try:
        run_pipeline(
            input_xtc=input_xtc,
            input_pdb=input_pdb,
            out_dir=out_dir,
            out_prefix=out_prefix,
            energy_threshold=args['threshold'],
            cores=num_cores,
            stride=args['stride'],
            module=args['module'],
            store_edges=args['edges'],
            verbose=verbose
        )
        logger.success('Done')
    except Exception as e:
        if args['log_level'] in ['TRACE']:
            raise e
        else:
            logger.critical(str(e))
