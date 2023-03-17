import numpy as np
import pytest

from lamberthub import battin1984


def test_exceeded_maximum_number_of_iterations():
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([1.0, 0.0, 0.0])  # [km]
    r2 = np.array([0.0, 1.0, 0.0])  # [km]
    tof = 1000  # [s]

    # Solving the problem with only two iteration so an error is raised
    with pytest.raises(ValueError) as excinfo:
        v1, v2 = battin1984(
            mu_earth, r1, r2, tof, maxiter=1, prograde=True, low_path=True
        )
    assert "Exceeded maximum number of iterations!" in excinfo.exconly()
