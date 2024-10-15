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
from TRAMbio.services.parameter import HydrogenBondParameter, GeneralWorkflowParameter
from TRAMbio.util.functions.argparse.base_parser import parse_args_for
from TRAMbio.util.structure_library.argparse import OptionsDictionary


_CLI_OPTIONS: Dict[str, OptionsDictionary] = {
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
    'threshold': OptionsDictionary(
        id=['-t', '--threshold'], args=dict(type=float, metavar='THRESHOLD', help=textwrap.dedent(
            """Energy threshold for inclusion of hydrogen bonds. All bonds with energy lower or equal to this threshold are included. (default: #DEFAULT#)
            """).replace("#DEFAULT#", str(ParameterRegistry.get_parameter_set('')(HydrogenBondParameter.ENERGY_THRESHOLD.value)))),
        default=lambda argv: ParameterRegistry.get_parameter_set('')(HydrogenBondParameter.ENERGY_THRESHOLD.value)),
    'log_level': OptionsDictionary(
        id=['-l', '--log-level'], args=dict(type=str, choices=[
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
        energy_threshold: float,
        store_edges: bool = False,
        verbose: bool = False
):
    parameter_id = f"TRAM_PDB_{time.perf_counter()}"

    parameter_registry = ParameterRegistry.get_parameter_set(parameter_id=parameter_id)
    parameter_registry.set_parameter(HydrogenBondParameter.ENERGY_THRESHOLD.value, energy_threshold)
    parameter_registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, verbose)

    out_file = os.path.join(out_dir, f"{out_prefix}_components.xml")
    temp_file = os.path.join(out_dir, f".temp_{out_prefix}_components.xml")
    edge_data_file = os.path.join(out_dir, f"{out_prefix}_edges.bnd") if store_edges else None

    pdb_workflow_service = WorkflowServiceRegistry.PDB.single_service()

    structure_generator = pdb_workflow_service.pdb_to_components(
        pdb_path=input_path,
        edge_data_file=edge_data_file,
        parameter_id=parameter_id
    )

    pdb_workflow_service.run_pipeline_on_generator(
        generator=structure_generator,
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
        'tram-pdb', "Calculate rigid components for PDB files.", _CLI_OPTIONS
    )

    ###################
    # Check input #####
    ###################

    input_file = args['pdb']
    if len(input_file) != 4 and (not os.path.exists(input_file) or not str(input_file).endswith('.pdb')):
        sys.exit("Input file needs to be a valid .pdb file.")

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
                     energy_threshold=args['threshold'],
                     store_edges=args['edges'],
                     verbose=verbose
                     )
        logger.success('Done')
    except Exception as e:
        if args['log_level'] in ['TRACE']:
            raise e
        else:
            logger.critical(str(e))
