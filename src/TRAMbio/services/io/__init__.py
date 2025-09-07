from TRAMbio.services.io.registry import IOServiceRegistry, \
    IPdbIOService, IXtcIOService, IXmlIOService, IBondIOService, IPyMolIOService, \
    AbstractPdbIOContext
from TRAMbio.util.errors import MissingDependencyError

try:
    import TRAMbio.services.io._pdb_io_service
except MissingDependencyError as mDE1:
    IOServiceRegistry.PDB.register_exception('PdbIOService', mDE1)
except ModuleNotFoundError:
    pass

import TRAMbio.services.io._xml_io_service

try:
    import TRAMbio.services.io._xtc_mdanalysis_io_service
except MissingDependencyError as mDE2:
    IOServiceRegistry.XTC.register_exception('MDAnalysis', mDE2)
except ModuleNotFoundError:
    pass

try:
    import TRAMbio.services.io._xtc_mdtraj_io_service
except MissingDependencyError as mDE3:
    IOServiceRegistry.XTC.register_exception('mdtraj', mDE3)
except ModuleNotFoundError:
    pass

import TRAMbio.services.io._bond_io_service
import TRAMbio.services.io._pymol_io_service
