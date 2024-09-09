"""A collection of hard tests only for robust solvers"""

import numpy as np
from numpy.testing import assert_allclose
import pytest

from lamberthub import ROBUST_SOLVERS

TABLE_OF_TRANSFERS = {
    "M0_prograde_high": [
        np.array([8.79257809, 0.27867677, 0.02581527]),  # [km / s]
        np.array([-8.68383320, -0.28592643, -0.03453010]),  # [km / s]
    ],
    "M1_prograde_high": [
        np.array([7.63353091, 0.24582764, 0.02569470]),  # [km / s]
        np.array([-7.50840227, -0.24335652, -0.02658981]),  # [km / s]
    ],
    "M1_prograde_low": [
        np.array([8.19519089, 2.30595215, 1.75229388]),  # [km / s]
        np.array([8.07984345, 2.30222567, 1.75189559]),  # [km / s]
    ],
    "M2_prograde_high": [
        np.array([6.51890385, 0.21496104, 0.02618989]),  # [km / s]
        np.array([-6.37230007, -0.20150975, -0.01832295]),  # [km / s]
    ],
    "M2_prograde_low": [
        np.array([7.00660748, 1.96687296, 1.49423471]),  # [km / s]
        np.array([6.87133644, 1.96250281, 1.49376762]),  # [km / s]
    ],
}
"""
Directly taken from example 2 from The Superior Lambert Algorithm (Der
Astrodynamics), by Gim J. Der, see https://amostech.com/TechnicalPapers/2011/Poster/DER.pdf
"""


@pytest.mark.parametrize("solver", ROBUST_SOLVERS)
@pytest.mark.parametrize("case", TABLE_OF_TRANSFERS)
def test_hard_case(solver, case):
    """
    Example 2a from The Superior Lambert Algorithm (Der Astrodynamics), by Gim
    J. Der, see: https://amostech.com/TechnicalPapers/2011/Poster/DER.pdf
    """
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([7231.58074563487, 218.02523761425, 11.79251215952])  # [km]
    r2 = np.array([7357.06485698842, 253.55724281562, 38.81222241557])  # [km]
    tof = 12300  # [s]

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
    assert_allclose(v1, expected_v1, rtol=1e-6)
    assert_allclose(v2, expected_v2, rtol=1e-6)
