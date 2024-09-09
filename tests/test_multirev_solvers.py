"""A collection of  tests only for multi-revolution solvers"""

import numpy as np
from numpy.testing import assert_allclose
import pytest

from lamberthub import MULTI_REV_SOLVERS

TABLE_OF_TRANSFERS = {
    "M1_prograde_high": [
        np.array([0.50335770, 0.61869408, -1.57176904]),  # [km / s]
        np.array([-4.18334626, -1.13262727, 6.13307091]),  # [km / s]
    ],
    "M1_prograde_low": [
        np.array([-2.45759553, 1.16945801, 0.43161258]),  # [km / s]
        np.array([-5.53841370, 0.01822220, 5.49641054]),  # [km / s]
    ],
    "M1_retrograde_high": [
        np.array([1.33645655, -0.94654565, 0.30211211]),  # [km / s]
        np.array([4.93628678, 0.39863416, -5.61593092]),  # [km / s]
    ],
    "M1_retrograde_low": [
        np.array([-1.38861608, -0.47836611, 2.21280154]),  # [km / s]
        np.array([3.92901545, 1.50871943, -6.52926969]),  # [km / s]
    ],
}
"""
Directly taken from example 1 from The Superior Lambert Algorithm (Der
Astrodynamics), by Gim J. Der, see https://amostech.com/TechnicalPapers/2011/Poster/DER.pdf
"""


@pytest.mark.parametrize("solver", MULTI_REV_SOLVERS)
@pytest.mark.parametrize("case", TABLE_OF_TRANSFERS)
def test_multirev_case(solver, case):
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([22592.145603, -1599.915239, -19783.950506])  # [km]
    r2 = np.array([1922.067697, 4054.157051, -8925.727465])  # [km]
    tof = 36000  # [s]

    # Unpack problem conditions
    M, sense, path = case.split("_")

    # Convert proper type
    M = int(M[-1])
    sense = True if sense == "prograde" else False
    path = True if path == "low" else False

    # Solve the problem
    v1, v2 = solver(mu_earth, r1, r2, tof, M=M, prograde=sense, low_path=path)

    # Expected final results
    expected_v1, expected_v2 = TABLE_OF_TRANSFERS[case]

    # Assert the results
    assert_allclose(v1, expected_v1, rtol=5e-6)
    assert_allclose(v2, expected_v2, rtol=5e-6)


@pytest.mark.parametrize("solver", MULTI_REV_SOLVERS)
def test_exception_try_lower_M(solver):
    """Test that solver does not find any solution for particular input"""
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([22592.145603, -1599.915239, -19783.950506])  # [km]
    r2 = np.array([1922.067697, 4054.157051, -8925.727465])  # [km]
    tof = 5 * 3600  # [s]

    with pytest.raises(ValueError) as excinfo:
        solver(mu_earth, r1, r2, tof, M=1)
    assert "ValueError: No feasible solution, try lower M!" in excinfo.exconly()
