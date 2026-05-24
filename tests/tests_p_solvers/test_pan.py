"""Tests for Pan and Ma's Bézier-function Lambert algorithm."""

import numpy as np
from numpy.testing import assert_allclose
import pytest

from lamberthub import pan2018


def test_pan2018_solves_zero_revolution_case():
    """Validate Pan and Ma's algorithm against a reference direct transfer."""
    mu_earth = 3.986004418e5
    r1 = np.array([15945.34, 0.0, 0.0])
    r2 = np.array([12214.83899, 10249.46731, 0.0])
    tof = 76.0 * 60

    v1, v2 = pan2018(mu_earth, r1, r2, tof)

    assert_allclose(v1, np.array([2.058913, 2.915965, 0.0]), atol=0.02, rtol=0.001)
    assert_allclose(v2, np.array([-3.451565, 0.910315, 0.0]), atol=0.02, rtol=0.02)


def test_pan2018_rejects_multirevolution_case():
    """Verify that Pan and Ma's solver reports unsupported revolutions."""
    r1 = np.array([1.0, 0.0, 0.0])
    r2 = np.array([0.0, 1.0, 0.0])

    with pytest.raises(ValueError) as excinfo:
        pan2018(1.0, r1, r2, 1.0, M=1)

    assert "Pan is not able to work within the multi-revolution scenario!" in str(
        excinfo.value
    )
