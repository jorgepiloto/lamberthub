import numpy as np
import pytest

from lamberthub import gauss1809


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
