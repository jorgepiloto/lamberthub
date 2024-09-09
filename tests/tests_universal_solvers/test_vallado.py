import numpy as np
import pytest

from lamberthub import vallado2013


def test_exception_for_180_transfer_angle():
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([1.0, 0.0, 0.0])  # [km]
    r2 = np.array([-1.0, 0.0, 0.0])  # [km]
    tof = 1000  # [s]

    # Solving the problem with only two iteration so an error is raised
    with pytest.raises(RuntimeError) as excinfo:
        v1, v2 = vallado2013(
            mu_earth, r1, r2, tof, maxiter=1, prograde=True, low_path=True
        )
    assert "Cannot compute orbit, phase angle is 180 degrees" in excinfo.exconly()


def test_raised_maximum_number_of_iterations():
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([1.0, 0.0, 0.0])  # [km]
    r2 = np.array([0.0, 1.0, 0.0])  # [km]
    tof = 1000  # [s]

    # Solving the problem with only two iteration so an error is raised
    with pytest.raises(ValueError) as excinfo:
        v1, v2 = vallado2013(
            mu_earth, r1, r2, tof, maxiter=1, prograde=True, low_path=True
        )
    assert "Exceeded maximum number of iterations!" in excinfo.exconly()
