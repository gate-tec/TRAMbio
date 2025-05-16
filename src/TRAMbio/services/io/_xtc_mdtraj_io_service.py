from typing import Generator, Tuple

import pandas as pd
import numpy as np
from io import StringIO

from TRAMbio.services.io import IXtcIOService, IOServiceRegistry
from TRAMbio.services import StructureServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import HydrogenBondParameter

from TRAMbio.util.errors import MissingDependencyError

try:
    import mdtraj as md
    from mdtraj.formats.pdb.pdbfile import PDBTrajectoryFile
    from mdtraj.core.trajectory import Trajectory
    from mdtraj.utils import (in_units_of)
    # Known compatible versions are [1.9.9, 1.10.0]
    # TODO: Double check compatible version (mdtraj has no __version__) from:
    # sys.modules["mdtraj"].version.full_version
except ModuleNotFoundError:
    exc = MissingDependencyError(
        module="TRAMbio.services.io._xtc_mdtraj_io_service",
        dependency="mdtraj"
    )
    try:
        raise exc
    finally:
        exc.__context__ = None


__all__ = []


class PDBTrajectoryOutStream(PDBTrajectoryFile):
    """
    Adapter-Class for utilizing streams as output source
    """

    def __init__(self):
        try:
            # pass mode='x' to prevent file system access
            super().__init__(None, mode='x')
        except ValueError:
            pass
        # recreate write mode initialization
        self._mode = 'w'
        self._header_written = False
        self._footer_written = False
        self._file = StringIO()
        self._open = True

    def close(self) -> "StringIO":
        if self._mode == 'w' and not self._footer_written:
            self._write_footer()
        if self._open:
            pass  # removed automatic closing of stream
        self._open = False
        return self._file


def to_pdb_stream(self, bfactors=None) -> "StringIO":
    self._check_valid_unitcell()

    if bfactors is not None:
        if len(np.array(bfactors).shape) == 1:
            if len(bfactors) != self.n_atoms:
                raise ValueError(
                    "bfactors %s should be shaped as (n_frames, n_atoms) or (n_atoms,)" % str(np.array(bfactors).shape))

            bfactors = [bfactors] * self.n_frames

        else:
            if np.array(bfactors).shape != (self.n_frames, self.n_atoms):
                raise ValueError(
                    "bfactors %s should be shaped as (n_frames, n_atoms) or (n_atoms,)" % str(np.array(bfactors).shape))

    else:
        bfactors = [None] * self.n_frames

    f = PDBTrajectoryOutStream()
    for i in range(self.n_frames):

        if self._have_unitcell:
            f.write(in_units_of(self._xyz[i], self._distance_unit, f.distance_unit),
                    self.topology,
                    modelIndex=i,
                    bfactors=bfactors[i],
                    unitcell_lengths=in_units_of(self.unitcell_lengths[i], self._distance_unit,
                                                 f.distance_unit),
                    unitcell_angles=self.unitcell_angles[i])
        else:
            f.write(in_units_of(self._xyz[i], self._distance_unit, f.distance_unit),
                    self.topology,
                    modelIndex=i,
                    bfactors=bfactors[i])

    return f.close()


def trajectory_adapter(trajectory: Trajectory):
    trajectory.to_pdb_stream = to_pdb_stream.__get__(trajectory)
    return trajectory


class XtcMdtrajIOService(IXtcIOService):

    @property
    def name(self):
        return "mdtraj"

    def read(
            self,
            xtc_path: str,
            pdb_path: str,
            stride: int,
            parameter_id: str = ''
    ) -> Tuple[int, Generator[Tuple[int, pd.DataFrame], None, None]]:

        parameter_registry = ParameterRegistry.get_parameter_set(parameter_id)
        pdb_structure_service = StructureServiceRegistry.PDB.single_service()
        pdb_io_service = IOServiceRegistry.PDB.single_service()

        actual_stride = stride
        if stride < 1:
            actual_stride = 1

        traj = md.load(xtc_path, top=pdb_path, stride=actual_stride)

        frames = list(range(traj.n_frames))
        if stride < 1:
            frames = frames[:1]

        num_frames = len(frames)

        def generator():
            for frame in frames:
                raw_df = pdb_io_service.read(
                    trajectory_adapter(traj[frame]).to_pdb_stream()
                ).export_first_model()
                if parameter_registry(HydrogenBondParameter.INCLUDE.value) and not pdb_structure_service.has_hydrogen_atoms(raw_df, parameter_id=parameter_id):
                    raise ValueError(f"No hydrogen atoms in ATOM records of {pdb_path}")
                yield frame * actual_stride, raw_df

        return num_frames, generator()


IOServiceRegistry.XTC.register_service(XtcMdtrajIOService())
