from typing import Generator, Tuple

import pandas as pd

from TRAMbio.services import StructureServiceRegistry, ParameterRegistry
from TRAMbio.services.parameter import HydrogenBondParameter
from TRAMbio.util.wrapper.base.stream import ToggledStringIO
from TRAMbio.util.wrapper.pandas.warning import WarningWrapper
from TRAMbio.services.io import IXtcIOService, IOServiceRegistry

from TRAMbio.util.errors import MissingDependencyError

try:
    import MDAnalysis as mda
    from MDAnalysis.coordinates.PDB import PDBWriter
except ModuleNotFoundError:
    exc = MissingDependencyError(
        module="TRAMbio.services.io._xtc_mdanalysis_io_service",
        dependency="MDAnalysis"
    )
    try:
        raise exc
    finally:
        exc.__context__ = None

# Turn off logging messages from MDAnalysis.coordinates.PBD
import logging
logging.getLogger("MDAnalysis.coordinates.PBD").disabled = True


__all__ = []


class XtcMDAnalysisIOService(IXtcIOService):

    @property
    def name(self):
        return "MDAnalysis"

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

        u = mda.Universe(pdb_path, xtc_path)
        traj = u.trajectory

        start = traj.frame if traj.n_frames > 1 else 0

        frames = list(range(start, len(traj), actual_stride))
        if stride < 1:
            frames = frames[:1]

        num_frames = len(frames)

        def generator() -> Generator[Tuple[int, pd.DataFrame], None, None]:
            for frame in frames:

                # Block closing in PDBWriter
                out_stream = ToggledStringIO()
                out_stream.toggle = True

                with WarningWrapper(UserWarning):
                    #################
                    # IMPORTANT #####
                    #################
                    # Manually increment trajectory frame
                    # as PDBWriter.write always writes "current" timestep
                    traj[frame]
                    #################

                    # Pass "non-closable" Stream as an in-memory file-like
                    w = PDBWriter(out_stream, multiframe=False, start=frame)
                    w.write(obj=u)
                    w.close()

                out_stream.seek(0)

                raw_df = pdb_io_service.read(out_stream).export_first_model()
                if parameter_registry(HydrogenBondParameter.INCLUDE.value) and not pdb_structure_service.has_hydrogen_atoms(raw_df, parameter_id=parameter_id):
                    raise ValueError(f"No hydrogen atoms in ATOM records of {pdb_path}")
                yield frame, raw_df

                out_stream.toggle = False
                # Manually close stream after data is collected
                out_stream.close()

            traj.close()
            # TODO: teardown for offset files?

        return num_frames, generator()


IOServiceRegistry.XTC.register_service(XtcMDAnalysisIOService())
