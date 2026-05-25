"""Tests for Jiang's semi-major-axis Lambert algorithm."""

import numpy as np
from numpy.testing import assert_allclose
import pytest

from lamberthub import jiang2016


def test_jiang2016_solves_zero_revolution_case():
    """Validate Jiang's algorithm against a reference direct transfer."""
    mu_earth = 3.986004418e5
    r1 = np.array([5000.0, 10000.0, 2100.0])
    r2 = np.array([-14600.0, 2500.0, 7000.0])
    tof = 3600

    v1, v2 = jiang2016(mu_earth, r1, r2, tof)

    assert_allclose(v1, np.array([-5.9925, 1.9254, 3.2456]), atol=0.02, rtol=0.001)
    assert_allclose(
        v2,
        np.array([-3.3125, -4.1966, -0.38529]),
        atol=0.02,
        rtol=0.001,
    )


def test_jiang2016_rejects_multirevolution_case():
    """Verify that Jiang's solver reports unsupported revolutions."""
    r1 = np.array([1.0, 0.0, 0.0])
    r2 = np.array([0.0, 1.0, 0.0])

    with pytest.raises(ValueError) as excinfo:
        jiang2016(1.0, r1, r2, 1.0, M=1)

    assert "Jiang is not able to work within the multi-revolution scenario!" in str(
        excinfo.value
    )
