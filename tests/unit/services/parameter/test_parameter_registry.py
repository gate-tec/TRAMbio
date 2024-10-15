import time
import pytest

from threading import Thread
from TRAMbio.services import lock_registry, ParameterRegistry
from TRAMbio.services.parameter import HydrogenBondParameter, GeneralWorkflowParameter

from numpy.ma.testutils import assert_equal

from TRAMbio.services.parameter.registry import verbosity_from_parameter
from tests.conftest import inject_pytest_logger

#######################
# Dummy Functions #####
#######################

class TestDummyFunctions:

    @pytest.fixture
    def locked_parameter(self):
        parameter_id = f"TRAM_PYTEST_PARAMETER_ID_{time.perf_counter()}"

        class DummyThread(Thread):

            def __init__(self, par_id):
                super().__init__()
                self._par_id = par_id
                self.local_lock = True

                @lock_registry(kwargs_name="par_id")
                def dummy_function(par_id: str):
                    # Theoretically use par_id here
                    while self.local_lock:
                        pass

                self._runnable = dummy_function

            def run(self):
                self._runnable(self._par_id)

        thread = DummyThread(parameter_id)
        thread.start()  # start parallel usage with same parameter_id

        yield parameter_id

        thread.local_lock = False


    @pytest.fixture
    def verbose_by_parameter(self):

        @verbosity_from_parameter()
        def dummy_function(verbose: bool = True, parameter_id: str = "TRAM_PYTEST_PARAMETER_ID"):
            return verbose

        yield dummy_function


#####################
# Mock Registry #####
#####################

class MockRegistry:

    @pytest.fixture
    def mock_registry_verbose(self):
        parameter_id = f"TRAM_PYTEST_PARAMETER_ID_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(parameter_id)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, True)

        yield parameter_id


    @pytest.fixture
    def mock_registry_non_verbose(self):
        parameter_id = f"TRAM_PYTEST_PARAMETER_ID_{time.perf_counter()}"
        registry = ParameterRegistry.get_parameter_set(parameter_id)
        registry.set_parameter(GeneralWorkflowParameter.VERBOSE.value, False)

        yield parameter_id


#############
# Tests #####
#############

class TestParameterLocked(TestDummyFunctions):

    def test_parameter_locked(self, locked_parameter, capsys):
        """Test that parameters can't be overwritten while locked function is running."""
        registry = ParameterRegistry.get_parameter_set(locked_parameter)
        initial_value = registry.get_parameter(HydrogenBondParameter.ENERGY_THRESHOLD.value)

        with inject_pytest_logger():
            registry.set_parameter(HydrogenBondParameter.ENERGY_THRESHOLD.value, initial_value+0.5)

        final_value = registry.get_parameter(HydrogenBondParameter.ENERGY_THRESHOLD.value)
        assert_equal(initial_value, final_value)

        captured = capsys.readouterr()[1]
        assert f"WARNING | Attempt to set \"{HydrogenBondParameter.ENERGY_THRESHOLD.value}\"" in captured

class TestVerbosityParameter(TestDummyFunctions, MockRegistry):

    # passed True + registry True => True
    def test_parameter_verbosity_by_value_true_1(self, verbose_by_parameter, mock_registry_verbose):
        assert verbose_by_parameter(verbose=True, parameter_id=mock_registry_verbose)

    # passed True + registry False => True
    def test_parameter_verbosity_by_value_true_2(self, verbose_by_parameter, mock_registry_non_verbose):
        assert verbose_by_parameter(verbose=True, parameter_id=mock_registry_non_verbose)

    # passed False + registry True => False
    def test_parameter_verbosity_by_value_false_1(self, verbose_by_parameter, mock_registry_verbose):
        assert not verbose_by_parameter(verbose=False, parameter_id=mock_registry_verbose)

    # passed False + registry False => False
    def test_parameter_verbosity_by_value_false_2(self, verbose_by_parameter, mock_registry_non_verbose):
        assert not verbose_by_parameter(verbose=False, parameter_id=mock_registry_non_verbose)

    # not passed + registry True => True
    def test_parameter_verbosity_by_registry_true(self, verbose_by_parameter, mock_registry_verbose):
        assert verbose_by_parameter(parameter_id=mock_registry_verbose)

    # not passed + registry False => False
    def test_parameter_verbosity_by_registry_false(self, verbose_by_parameter, mock_registry_non_verbose):
        assert not verbose_by_parameter(parameter_id=mock_registry_non_verbose)