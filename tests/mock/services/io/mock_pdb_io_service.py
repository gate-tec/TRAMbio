import operator
from io import StringIO
from typing import Generator, Union

import pytest

from TRAMbio.services.io import IPdbIOService, AbstractPdbIOContext
from TRAMbio.util.wrapper.biopandas.pandas_pdb import CustomPandasPdb
from tests.mock import MockServiceEvalException


class MockPdbIOService(IPdbIOService):

    def __init__(self, mock_data: CustomPandasPdb, expected_atom_count: int = None):
        self.mock_data = mock_data
        self.expected_atom_count = expected_atom_count
        self.calls = 0
        self.input_hashes = set()

    @property
    def name(self):
        return "PdbIOService"

    def read(self, input_data: Union[str, StringIO], verbose: bool = True) -> CustomPandasPdb:
        self.calls += 1
        __tracebackhide__ = operator.methodcaller("errisinstance", MockServiceEvalException)
        if not isinstance(input_data, StringIO):
            raise MockServiceEvalException(f"Wrong input type passed to MockService. Expected StringIO, got {type(input_data)} instead.")

        input_data.seek(0)
        atom_lines = [line for line in input_data.readlines() if line.startswith("ATOM")]

        if self.expected_atom_count is not None:
            num_atoms = len(atom_lines)

            if num_atoms != self.expected_atom_count:
                raise MockServiceEvalException(f"Wrong number of atoms loaded in trajectory. Expected {self.expected_atom_count}, got {num_atoms} instead.")

        self.input_hashes.add(hash("".join(atom_lines)) % (10 ** 8))

        return self.mock_data

    def get_result_counter(self):
        return self.calls, len(self.input_hashes)

    def pdb_file_context(self, pdb_path: str, header_stream: StringIO) -> Generator[AbstractPdbIOContext, None, None]:
        raise NotImplementedError
