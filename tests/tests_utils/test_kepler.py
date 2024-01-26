"""Unitary tests related to Kepler's equation"""

from numpy.testing import assert_allclose
import pytest

from lamberthub.utils.kepler import kepler_from_nu

TABLE_OF_SOLUTIONS = {
    "elliptic": [
        0.62749,
        0.80521,
        1.00220,
        0.24649,
    ],  # [M, E, nu, ecc] from Curtis 3.3
    "parabolic": [
        2 * 6.7737,
        3.1481,
        2.526364092261792,
        1,
    ],  # [M, B, nu, ecc] from Curtis 3.4
    "hyperbolic": [
        11.279,
        2.2927,
        1.7453292519943295,
        2.7696,
    ],  # [M, H, nu, ecc] from Curtis 3.5
}


@pytest.mark.parametrize("orbit_type", TABLE_OF_SOLUTIONS)
def test_kepler_from_nu(orbit_type):
    # Unpack expected data and initial one
    M_expected, _, nu, ecc = TABLE_OF_SOLUTIONS[orbit_type]

    # Compute the numerical predicted value
    M = kepler_from_nu(nu, ecc)

    # Check both quantities match
    assert_allclose(M, M_expected, atol=1e-2, rtol=1e-4)


def test_kepler_from_nu_raises_error_null_ecc():
    # Define a negative eccentricity at periapsis location
    ecc, nu = -0.5, 0
    with pytest.raises(ValueError) as excinfo:
        kepler_from_nu(nu, ecc)
    assert "Eccentricity cannot be negative!" in excinfo.exconly()
