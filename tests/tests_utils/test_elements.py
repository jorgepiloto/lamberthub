import numpy as np
import pytest
from numpy.testing import assert_allclose

from lamberthub.utils.elements import rotation_matrix, rv2coe

DEG_TO_RAD = np.pi / 180


# This table was generated using NASA's software GMAT-2020a version
TABLE_OF_ELEMENTS_CONVERSIONS = {
    "non-equatorial-non-circular": [
        [
            np.array([10000, 0, 0]),  # [km]
            np.array([1.00, 5.171700990602514, 3.621264018986504]),  # [km / s]
        ],
        [
            10000.00002042279,
            0.1583912232069116,
            35 * DEG_TO_RAD,
            0,
            270.0000007387653 * DEG_TO_RAD,
            89.99999926123472 * DEG_TO_RAD,
        ],
    ],
    "equatorial-non-circular": [
        [
            np.array([10000, 0, 0]),  # [km]
            np.array([0, 7.35, 0]),  # [km / s]
        ],
        [
            13553.04570077853,
            0.355304570077853,
            0,
            0,
            0,
            0,
        ],
    ],
    "non-equatorial-circular": [
        [
            np.array([10000, 0, 0]),  # [km]
            np.array([0, 5.171700990602514, 3.621264018986504]),  # [km / s]
        ],
        [
            10000.00002042279,
            2.042278850497544e-09,
            35 * DEG_TO_RAD,
            0,
            0,
            0,
        ],
    ],
    "equatorial-circular": [
        [
            np.array([10000, 0, 0]),  # [km]
            np.array([0, 6.313481143553055, 0]),  # [km / s]
        ],
        [
            10000.00002042279,
            2.042278850497544e-09,
            0,
            0,
            0,
            0,
        ],
    ],
    "non-equatorial-hyperbolic": [
        [
            np.array([10000, 1250, 0]),  # [km]
            np.array([0, 10, 3.50]),  # [km / s]
        ],
        [
            28209.05210914074,
            1.829130032278414,
            19.4289577238344 * DEG_TO_RAD,
            7.125016348901807 * DEG_TO_RAD,
            349.6072413039819 * DEG_TO_RAD,
            10.39275869601813 * DEG_TO_RAD,
        ],
    ],
}


@pytest.mark.parametrize("orbit_type", TABLE_OF_ELEMENTS_CONVERSIONS)
def test_rv2coe_non_equatorial(orbit_type):
    # Unpack position and velocity vectors
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r, v = TABLE_OF_ELEMENTS_CONVERSIONS[orbit_type][0]

    # Compute predicted values
    p, ecc, inc, raan, argp, nu = rv2coe(mu_earth, r, v)

    # Expected values for COE
    (
        expected_p,
        expected_ecc,
        expected_inc,
        expected_raan,
        expected_argp,
        expected_nu,
    ) = TABLE_OF_ELEMENTS_CONVERSIONS[orbit_type][1]

    # Check all values
    assert_allclose(p, expected_p)
    assert_allclose(ecc, expected_ecc, atol=1e-8)
    assert_allclose(inc, expected_inc)
    assert_allclose(raan, expected_raan)
    assert_allclose(argp, expected_argp)
    assert_allclose(nu, expected_nu)


def test_rotation_matrix_x():
    result = rotation_matrix(0.218, 0)
    expected = np.array(
        [[1.0, 0.0, 0.0], [0.0, 0.97633196, -0.21627739], [0.0, 0.21627739, 0.97633196]]
    )
    assert np.allclose(expected, result)


def test_rotation_matrix_y():
    result = rotation_matrix(0.218, 1)
    expected = np.array(
        [[0.97633196, 0.0, 0.21627739], [0.0, 1.0, 0.0], [0.21627739, 0.0, 0.97633196]]
    )
    assert np.allclose(expected, result)


def test_rotation_matrix_z():
    result = rotation_matrix(0.218, 2)
    expected = np.array(
        [[0.97633196, -0.21627739, 0.0], [0.21627739, 0.97633196, 0.0], [0.0, 0.0, 1.0]]
    )
    assert np.allclose(expected, result)


def test_rotation_matrix_wrong_axis():
    with pytest.raises(ValueError) as excinfo:
        rotation_matrix(0.218, 3)
    assert "Invalid axis: must be one of 'x', 'y' or 'z'" in excinfo.exconly()
