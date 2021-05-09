import numpy as np
import pytest
from numpy.testing import assert_allclose

from lamberthub.utils.angles import B_to_nu, E_to_nu, H_to_nu, get_transfer_angle


@pytest.mark.parametrize("sense", [True, False])
def test_get_transfer_angle_collinear_vectors_same_sense(sense):

    # Build two vectors with the same direction and sense using scalar
    # proportion
    r1 = np.array([1, 0, 0])
    r2 = r1 if sense is True else -r1

    # Compute the transfer angle
    dtheta = get_transfer_angle(r1, r2, True)
    expected_dtheta = 0 if sense is True else np.pi

    # Check that transfer angle is zero
    assert_allclose(dtheta, expected_dtheta)


TABLE_OF_SOLUTIONS = {
    "elliptic": [0.80521, 1.00220, 0.24649],  # [E, nu, ecc] from Curtis 3.3
    "parabolic": [3.1481, 2.526364092261792, 1],  # [B, nu, ecc] from Curtis 3.4
    "hyperbolic": [2.2927, 1.7453292519943295, 2.7696],  # [H, nu, ecc] from Curtis 3.5
}


@pytest.mark.parametrize("orbit_type", TABLE_OF_SOLUTIONS)
def test_from_orbit_anomaly_to_true_anomaly(orbit_type):
    # Unpack expected values
    orbit_anomaly, expected_nu, ecc = TABLE_OF_SOLUTIONS[orbit_type]

    # Solve for the predicted numerical value
    if 0 <= ecc < 1:
        nu = E_to_nu(orbit_anomaly, ecc)
    elif ecc == 1:
        nu = B_to_nu(orbit_anomaly)
    else:
        nu = H_to_nu(orbit_anomaly, ecc)

    # Check values match
    assert_allclose(nu, expected_nu, rtol=1e-4)
