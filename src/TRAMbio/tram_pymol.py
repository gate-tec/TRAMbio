#!/usr/bin/env python
import os
import sys
import textwrap
import time
from typing import Optional, Dict, Literal

from pathvalidate import sanitize_filename

from TRAMbio import set_log_level
from TRAMbio.services import IOServiceRegistry, ParameterRegistry, WorkflowServiceRegistry
from TRAMbio.services.parameter import XtcParameter, GeneralWorkflowParameter
from TRAMbio.util.functions.argparse.base_parser import parse_args_for
from TRAMbio.util.structure_library.argparse import OptionsDictionary

from loguru import logger

_CLI_OPTIONS: Dict[str, OptionsDictionary] = {
    'pdb': OptionsDictionary(
        id=['-p', '--pdb'], args=dict(type=str, metavar='PDB_FILE', help=textwrap.dedent(
            """Protein input file in PDB v3.3 format.
            """)), default=None),
    'xtc': OptionsDictionary(
        id=['--xtc'], args=dict(type=str, metavar='XTC_FILE', help=textwrap.dedent(
            """Trajectory file in XTC format.
            """)), default=lambda argv: None),
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
    'xml_file': OptionsDictionary(
        id=['-x', '--xml'], args=dict(type=str, metavar='XML_FILE', help=textwrap.dedent(
            """XML-type component file from the component analysis.
            """)),
        default=lambda argv: os.path.join(argv['out_dir'], f"{argv['name']}_components.xml")),
    'bnd_file': OptionsDictionary(
        id=['-b', '--bnd-file'], args=dict(type=str, metavar='BOND_FILE', help=textwrap.dedent(
            """Bond list file (ending .bnd) from the component analysis.
            """)),
        default=lambda argv: None),
    'key': OptionsDictionary(
        id=['-k', '--key'], args=dict(type=str, metavar='MIN_KEY', help=textwrap.dedent(
            """Minimum value for state keys. Either a float for single frame PDBs, indicating the minimum strength for present hydrogen bonds, or a starting framenumber for trajectories (default: no limit)
            """)),
        default=lambda argv: None),
    'max_states': OptionsDictionary(
        id=['-s', '--max-states'], args=dict(type=int, metavar='MAX_STATES', help=textwrap.dedent(
            """Maximum number of states to visualize. Values lower than 1 indicate no limit. (default: 0)
            """)), default=lambda argv: 0),
    'module': OptionsDictionary(
        id=['-m', '--module'], args=dict(type=str, choices=IOServiceRegistry.XTC.list_service_names(), help=textwrap.dedent(
            """Base module for trajectory loading. Only important when visualizing Trajectories. (default: MDAnalysis)
            """)),
        default=lambda argv: 'MDAnalysis'),
    'log_level': OptionsDictionary(
        id=['-l', '--log-level'], args=dict(type=str, choices=[
            'TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING', 'ERROR', 'CRITICAL', 'NONE'
        ], help=textwrap.dedent(
            """Set the logging level during execution. (default: INFO)
            """)),
        default=lambda argv: 'INFO')
}


def pipeline(
        pdb_path: str,
        xtc_path: Optional[str],
        xml_path: str,
        bond_path: Optional[str],
        out_dir: str,
        out_prefix: str,
        module: str,
        verbose: bool = True
):
    parameter_id = f"TRAM_PYMOL_{time.perf_counter()}"

    out_pdb_file = f'{out_prefix}_components.pdb'
    out_pdb_path = os.path.join(out_dir, out_pdb_file)
    out_pml_file = f'{out_prefix}_components.pml'
    out_pml_path = os.path.join(out_dir, out_pml_file)

    pymol_workflow_service = WorkflowServiceRegistry.PYMOL.single_service()
    parameter_registry = ParameterRegistry.get_parameter_set(parameter_id=parameter_id)
    parameter_registry.set_parameter(XtcParameter.MODULE.value, module)
    parameter_registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, verbose)

    pymol_workflow_service.convert_to_pymol_files(
        in_pdb_path=pdb_path,
        in_xtc_path=xtc_path,
        in_xml_path=xml_path,
        in_bond_path=bond_path,
        out_pdb_path=out_pdb_path,
        out_pml_path=out_pml_path,
        out_prefix=out_prefix,
        rel_out_pdb_path=out_pdb_file,
        parameter_id=parameter_id
    )


def main(default_log_level: Literal['TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING', 'ERROR', 'CRITICAL', 'NONE'] = "INFO"):
    set_log_level(default_log_level)

    args = parse_args_for(
        'tram-pymol', "Create PyMol visualization scripts from component analysis results.", _CLI_OPTIONS
    )

    ###################
    # Check input #####
    ###################

    pdb_file = args['pdb']
    if len(pdb_file) != 4 and (not os.path.exists(pdb_file) or not str(pdb_file).endswith('.pdb')):
        sys.exit("Input file needs to be a valid .pdb file.")

    xtc_file = args['xtc']
    if xtc_file is not None and (not os.path.exists(xtc_file) or not str(xtc_file).endswith('.xtc')):
        sys.exit("Trajectory file needs to be a valid .xtc file.")

    ####################
    # Check output #####
    ####################

    out_dir = args['out_dir']
    if os.path.isfile(out_dir):
        sys.exit(f'Given output directory "{out_dir}" is an existing filename.')
    out_prefix = sanitize_filename(args['name']).replace(' ', '_')

    ######################
    # Check XML file #####
    ######################

    xml_file = args['xml_file']
    if not os.path.exists(xml_file):
        sys.exit(f'Components XML "{xml_file}" not found.')

    ######################
    # Check bnd file #####
    ######################

    bnd_file = args['bnd_file']
    if bnd_file is not None and not os.path.exists(xml_file):
        sys.exit(f'Bond-list file "{bnd_file}" not found.')

    ####################
    # Run pipeline #####
    ####################

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    set_log_level(args['log_level'])
    verbose = args['log_level'] in ['TRACE', 'DEBUG', 'INFO']  # verbosity argument for non-loguru messages

    if xtc_file is not None and bnd_file is not None:
        logger.warning("Additional hydrogen bond visualization is not advised for Trajectories.")
        bnd_file = None

    try:
        pipeline(
            pdb_path=pdb_file,
            xtc_path=xtc_file,
            xml_path=xml_file,
            bond_path=bnd_file,
            out_dir=out_dir,
            out_prefix=out_prefix,
            module=args['module'],
            verbose=verbose
        )
        logger.success('Done')
    except Exception as e:
        if args['log_level'] in ['TRACE']:
            raise e
        else:
            logger.critical(str(e))
