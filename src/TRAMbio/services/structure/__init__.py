from TRAMbio.services.structure.registry import StructureServiceRegistry, IPdbStructureService, IXmlStructureService
from TRAMbio.util.errors import MissingDependencyError

try:
    import TRAMbio.services.structure._pdb_structure_service
except MissingDependencyError as mDE:
    StructureServiceRegistry.PDB.register_exception("PdbStructureService", mDE)

import TRAMbio.services.structure._xml_structure_service
