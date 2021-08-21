import numpy as np
import pytest
from numpy.testing import assert_allclose

from lamberthub import gauss1809
from lamberthub.p_solvers.gauss import _X_at_x


def test_exceeded_maximum_number_of_iterations():
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([1.0, 0.0, 0.0])  # [km]
    r2 = np.array([0.0, 1.0, 0.0])  # [km]
    tof = 1000  # [s]

    # Solving the problem
    with pytest.raises(ValueError) as excinfo:
        v1, v2 = gauss1809(mu_earth, r1, r2, tof, prograde=True, low_path=True)
    assert "Exceeded maximum number of iterations." in excinfo.exconly()


def test_singular_for_180_degrees_transfer():
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([1.0, 0.0, 0.0])  # [km]
    r2 = -r1
    tof = 1000  # [s]

    # Solving the problem
    with pytest.raises(ValueError) as excinfo:
        v1, v2 = gauss1809(mu_earth, r1, r2, tof, prograde=True, low_path=True)
    assert "Transfer angle was found to be 180 degrees!" in excinfo.exconly()


@pytest.mark.parametrize("x", np.round(np.linspace(0, 1, 10 + 1), 1))
def test_X_at_x(x):

    # Use series expansion from lamberthub
    X = _X_at_x(x, order=5)

    # Use explicit formula for computing the value using three terms
    expected_X = (4 / 3) * (
        1
        + (6 / 5) * x
        + (6 * 8) / (5 * 7) * x ** 2
        + (6 * 8 * 10) / (5 * 7 * 9) * x ** 3
        + (6 * 8 * 10 * 12) / (5 * 7 * 9 * 11) * x ** 4
        + (6 * 8 * 10 * 12 * 14) / (5 * 7 * 9 * 11 * 13) * x ** 5
    )

    # Assert if both numbers raise the same value
    assert_allclose(X, expected_X)
