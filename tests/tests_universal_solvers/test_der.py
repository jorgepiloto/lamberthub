"""Tests for Der's superior Lambert algorithm."""

import numpy as np
from numpy.testing import assert_allclose

from lamberthub import der2011


def test_der2011_solves_der_article_multirevolution_case():
    """Validate Der's algorithm against a reference multirevolution case."""
    mu_earth = 3.986004418e5
    r1 = np.array([22592.145603, -1599.915239, -19783.950506])
    r2 = np.array([1922.067697, 4054.157051, -8925.727465])
    tof = 36000

    v1, v2 = der2011(
        mu_earth,
        r1,
        r2,
        tof,
        M=1,
        is_prograde=True,
        is_low_path=False,
    )

    assert_allclose(v1, np.array([0.50335770, 0.61869408, -1.57176904]), rtol=5e-6)
    assert_allclose(v2, np.array([-4.18334626, -1.13262727, 6.13307091]), rtol=5e-6)
