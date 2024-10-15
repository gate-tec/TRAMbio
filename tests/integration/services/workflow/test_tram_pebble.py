import os

import pytest
import sys

from TRAMbio.tram_pebble import main


#################
# Arguments #####
#################

class TestArguments:

    @pytest.fixture
    def args_simple_2_3(self, tmp_path, path_graphml_simple):
        orig_argv = list(sys.argv)
        d = tmp_path / "pebble_sample"
        sys.argv = [
            "tram-pebble",
            "--graph", str(path_graphml_simple),
            "--out-dir", str(d),
            "--k-param", "2",
            "--l-param", "3",
            "--cores", "1",
            "--name", "sample",
            "--log-level", "TRACE"
        ]

        yield d

        sys.argv = orig_argv

        # clean-up files
        for file in os.listdir(d):
            os.remove(d / file)
        os.rmdir(d)

    @pytest.fixture
    def args_directed_2_3(self, tmp_path, path_graphml_directed):
        orig_argv = list(sys.argv)
        d = tmp_path / "pebble_sample"
        sys.argv = [
            "tram-pebble",
            "--graph", str(path_graphml_directed),
            "--out-dir", str(d),
            "--k-param", "2",
            "--l-param", "3",
            "--cores", "1",
            "--name", "sample",
            "--log-level", "TRACE"
        ]

        yield d

        sys.argv = orig_argv

        # clean-up files
        for file in os.listdir(d):
            os.remove(d / file)
        os.rmdir(d)

    @pytest.fixture
    def args_broken_trace(self, tmp_path, path_graphml_broken):
        orig_argv = list(sys.argv)
        d = tmp_path / "pebble_sample"
        sys.argv = [
            "tram-pebble",
            "--graph", str(path_graphml_broken),
            "--out-dir", str(d),
            "--k-param", "2",
            "--l-param", "3",
            "--cores", "1",
            "--name", "sample",
            "--log-level", "TRACE"
        ]

        yield d

        sys.argv = orig_argv

        # clean-up files
        for file in os.listdir(d):
            os.remove(d / file)
        os.rmdir(d)


#############
# Tests #####
#############

class TestTramPebble(TestArguments):

    def test_tram_pebble_simple(self, args_simple_2_3):
        main()

        assert len(list(args_simple_2_3.iterdir())) == 1  # result XML

    def test_tram_pebble_directed(self, args_directed_2_3):
        main()

        assert len(list(args_directed_2_3.iterdir())) == 1  # result XML


#####################################
# Test Error Handling & Logging #####
#####################################

class TestTramPebbleLogging(TestArguments):

    def test_error_with_trace(self, args_broken_trace):
        with pytest.raises(ValueError, match="Unable to load"):
            main()
