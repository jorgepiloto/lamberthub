""" Tests regarding package metadata """

import lamberthub


def test_version():
    assert lamberthub.__version__ == "0.2.dev0"
