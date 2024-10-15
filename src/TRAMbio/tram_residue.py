#!/usr/bin/env python
import os
import sys
import textwrap
import time
from typing import Optional, Dict, Literal

from pathvalidate import sanitize_filename

from loguru import logger

from TRAMbio import set_log_level
from TRAMbio.util.functions.argparse.base_parser import parse_args_for
from TRAMbio.util.structure_library.argparse import OptionsDictionary

from TRAMbio.services import WorkflowServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import ResidueParameter, GeneralWorkflowParameter

_CLI_OPTIONS: Dict[str, OptionsDictionary] = {
    'xml_file': OptionsDictionary(
        id=['-x', '--xml'], args=dict(type=str, metavar='XML_FILE', help=textwrap.dedent(
            """XML-type component file from the component analysis.
            """)),
        default=None),
    'pdb_file': OptionsDictionary(
        id=['-p', '--pdb'], args=dict(type=str, metavar='PDB_FILE', help=textwrap.dedent(
            """Reference protein input file in PDB v3.3 format. (default: None)
            """)),
        default=lambda x: None),
    'name': OptionsDictionary(
        id=['-n', '--name'], args=dict(type=str, metavar='PDB_NAME', help=textwrap.dedent(
            """Alternate name for protein (used as output prefix). If not specified, name is derived from input file name.
            """)),
        default=lambda argv: os.path.basename(argv['xml_file'])
        if '.' not in os.path.basename(argv['xml_file']) else
        os.path.basename(argv['xml_file'])[:os.path.basename(argv['xml_file']).rindex('.')]),
    'out_dir': OptionsDictionary(
        id=['-o', '--out-dir'], args=dict(type=str, metavar='OUTPUT_DIR', help=textwrap.dedent(
            """Directory for output files. (default: next to input file)
            """)),
        default=lambda argv: os.path.dirname(os.path.abspath(argv['xml_file']))),
    'key': OptionsDictionary(
        id=['-k', '--key'], args=dict(type=str, metavar='MIN_KEY', help=textwrap.dedent(
            """Minimum value for state keys. Either a float for single frame PDBs, indicating the minimum strength for present hydrogen bonds, or a starting framenumber for trajectories (default: no limit)
            """)),
        default=lambda argv: None),
    'max_states': OptionsDictionary(
        id=['-m', '--max-states'], args=dict(type=int, metavar='MAX_STATES', help=textwrap.dedent(
            """Maximum number of states to visualize. Values lower than 1 indicate no limit. (default: 0)
            """)), default=lambda argv: 0),
    'log_level': OptionsDictionary(
        id=['-l', '--log-level'], args=dict(type=str, choices=[
            'TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING', 'ERROR', 'CRITICAL', 'NONE'
        ], help=textwrap.dedent(
            """Set the logging level during execution. (default: INFO)
            """)),
        default=lambda argv: 'INFO')
}


def pipeline(
        xml_path: str,
        pdb_path: Optional[str],
        out_dir: str,
        out_prefix: str,
        min_key: Optional[str],
        max_states: int,
        verbose: bool = True
):
    parameter_id = f"TRAM_RESIDUE_{time.perf_counter()}"
    out_file = os.path.join(out_dir, f"{out_prefix}_residue_components.json")

    parameter_registry = ParameterRegistry.get_parameter_set(parameter_id)
    residue_workflow_service = WorkflowServiceRegistry.RESIDUE.single_service()

    if min_key is not None:
        parameter_registry.set_parameter(ResidueParameter.MIN_KEY.value, min_key)
    parameter_registry.set_parameter(ResidueParameter.MAX_STATES.value, max_states)
    parameter_registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, verbose)

    residue_workflow_service.convert_to_residue_level(
        xml_path=xml_path, pdb_path=pdb_path, out_file=out_file, parameter_id=parameter_id
    )


def main(default_log_level: Literal['TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING', 'ERROR', 'CRITICAL', 'NONE'] = "INFO"):
    set_log_level(default_log_level)

    args = parse_args_for(
        'tram-residue', "Create residue-level components from atom-level component analysis results.", _CLI_OPTIONS
    )

    ######################
    # Check XML file #####
    ######################

    xml_file = args['xml_file']
    if not os.path.exists(xml_file):
        sys.exit(f'Components XML "{xml_file}" not found.')

    pdb_file = args['pdb_file']
    if pdb_file and not os.path.exists(pdb_file):
        sys.exit(f'PDB file "{pdb_file}" not found.')

    ####################
    # Check output #####
    ####################

    out_dir = args['out_dir']
    if os.path.isfile(out_dir):
        sys.exit(f'Given output directory "{out_dir}" is an existing filename.')
    out_prefix = sanitize_filename(args['name']).replace(' ', '_')

    if out_prefix.endswith('_components'):
        out_prefix = out_prefix[:-len('_components')]

    ####################
    # Run pipeline #####
    ####################

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    set_log_level(args['log_level'])
    verbose = args['log_level'] in ['TRACE', 'DEBUG', 'INFO']  # verbosity argument for non-loguru messages

    try:
        pipeline(
            xml_path=xml_file,
            pdb_path=pdb_file,
            out_dir=out_dir,
            out_prefix=out_prefix,
            min_key=args['key'],
            max_states=args['max_states'],
            verbose=verbose
        )
        logger.success('Done')
    except Exception as e:
        if args['log_level'] in ['TRACE']:
            raise e
        else:
            logger.critical(str(e))
