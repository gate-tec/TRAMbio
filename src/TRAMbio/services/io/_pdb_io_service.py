import shutil
from contextlib import contextmanager
from io import StringIO
from typing import Union, TextIO, Generator
import os
from urllib.error import URLError, HTTPError

from TRAMbio.services.io import IPdbIOService, IOServiceRegistry, AbstractPdbIOContext
from loguru import logger

from TRAMbio.util.wrapper.biopandas.pandas_pdb import CustomPandasPdb


__all__ = []

from TRAMbio.util.wrapper.pandas.warning import WarningWrapper


class PdbIOContext(AbstractPdbIOContext):

    def __init__(self, pdb_stream: TextIO):
        self.__pdb_stream = pdb_stream

    def write_model(self, model: CustomPandasPdb, model_idx: int) -> None:
        self.__pdb_stream.write(f"MODEL     {model_idx:4d}".ljust(80, ' ') + '\n')
        out_stream = model.to_pdb_stream(records=['ATOM', 'HETATM', 'OTHERS'])
        shutil.copyfileobj(out_stream, self.__pdb_stream, -1)  # noqa
        self.__pdb_stream.write('ENDMDL'.ljust(80, ' ') + '\n')


class PdbIOService(IPdbIOService):

    @property
    def name(self):
        return "PdbIOService"

    def read(self, input_data: Union[str,  StringIO], verbose: bool = True) -> CustomPandasPdb:
        data: CustomPandasPdb
        if isinstance(input_data, StringIO):
            input_data.seek(0)
            with WarningWrapper(UserWarning):
                data =  CustomPandasPdb().read_pdb_from_list([
                    line + '\n' for line in input_data.getvalue().split('\n') if len(line) > 0
                ])
        elif os.path.exists(input_data) and os.path.isfile(input_data):
            if verbose:
                logger.info(f'Loading "{input_data}"...')
            with WarningWrapper(UserWarning):
                data =  CustomPandasPdb().read_pdb(path=input_data)
        elif len(input_data) == 4:
            if verbose:
                logger.info(f'Fetching "{input_data}"...')
            # PDB code
            try:
                data = CustomPandasPdb().fetch_pdb(pdb_code=input_data)
            except HTTPError as e:
                if str(e.code) == "404":
                    err = ValueError(f"Unknown PDB-code {input_data}")
                    try:
                        raise err
                    finally:
                        err.__context__ = None
                err = ValueError(f"HTTP Error {e.code} while fetching PDB.")
                try:
                    raise err
                finally:
                    err.__context__ = None
            except URLError as e:
                if "Name or service not known" in str(e.args):
                    err = ValueError("No network connection.")
                    try:
                        raise err
                    finally:
                        err.__context__ = None
                err = ValueError(f"Unknown PDB-code {input_data}")
                try:
                    raise err
                finally:
                    err.__context__ = None
            except AttributeError:
                err = ValueError(f"Unknown PDB-code {input_data}")
                try:
                    raise err
                finally:
                    err.__context__ = None
        else:
            raise ValueError(f"Unable to match PDB path input: {input_data}.")

        # TODO: quality check

        return data


    @contextmanager
    def pdb_file_context(
            self,
            pdb_path: str,
            header_stream:  StringIO
    ) -> Generator[AbstractPdbIOContext, None, None]:
        with open(pdb_path, 'w') as pdb_f:
            header_stream.seek(0)
            shutil.copyfileobj(header_stream, pdb_f, -1)  # noqa
            try:
                yield PdbIOContext(pdb_f)
            finally:
                pdb_f.write('END'.ljust(80, ' ') + '\n')


IOServiceRegistry.PDB.register_service(PdbIOService())
