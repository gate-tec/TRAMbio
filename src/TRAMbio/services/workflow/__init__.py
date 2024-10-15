from TRAMbio.services.workflow.registry import WorkflowServiceRegistry, IBaseWorkflowService, IPdbWorkflowService, \
    IXtcWorkflowService, IResidueWorkflowService, IPyMolWorkflowService, IPebbleGameWorkflowService
from TRAMbio.util.errors import MissingDependencyError
from TRAMbio.services.workflow._base_workflow_service import BaseWorkflowService
import TRAMbio.services.workflow._pebble_game_workflow_service
import TRAMbio.services.workflow._pdb_workflow_service
import TRAMbio.services.workflow._xtc_workflow_service
import TRAMbio.services.workflow._residue_workflow_service
try:
    import TRAMbio.services.workflow._pymol_workflow_service
except MissingDependencyError as mDE:
    WorkflowServiceRegistry.PYMOL.register_exception("PyMolWorkflowService", mDE)
