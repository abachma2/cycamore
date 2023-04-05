import os
import subprocess
import pytest


import run_inputs as ri
from helper import cyclus_has_coin


def coin_skipper(filename):
    raise pytest.skip(filename + " cannot be executed since Cyclus was not installed "
                   "with COIN support")


def test_inputs():
    files, _, _ = ri.get_files(ri.input_path)
    for f in files:
        absfile = os.path.join(ri.input_path, f)
        with open(absfile) as fh:
            src = fh.read()
        if cyclus_has_coin() or "GrowthRegion" not in src:
            testf = ri.TestFile(ri.cyclus_path, f, "-v0")
            testf.run()
            assert testf.passed, "Failed running {}".format(f)
        else:
            yield coin_skipper, absfile
