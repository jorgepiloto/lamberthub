"""Basic unitary tests for checking all available solvers"""

import numpy as np
from numpy.testing import assert_allclose
import pytest

from lamberthub import ALL_SOLVERS as ALL_SOLVERS_LAMBERTHUB, NON_ROBUST_SOLVERS

ALL_SOLVERS = [
    solver for solver in ALL_SOLVERS_LAMBERTHUB if solver not in NON_ROBUST_SOLVERS
]
"""A list holding all solvers which present a decent accuracy. Most of the
solvers from the NON_ROBUST_SOLVERS macro were developed in the old days."""

(ATOL, RTOL) = (0.02, 0.001)
"""The absolute and relative tolerances for the tests."""


@pytest.mark.parametrize("solver", ALL_SOLVERS_LAMBERTHUB)
def test_case_from_vallado_book(solver):
    """
    Example 5.7 from Fundamentals of Astrodynamics and Applications (4th
    Edition), by David A. Vallado
    """
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([15945.34, 0.0, 0.0])  # [km]
    r2 = np.array([12214.83899, 10249.46731, 0.0])  # [km]
    tof = 76.0 * 60  # [s]

    # Solving the problem
    v1, v2 = solver(mu_earth, r1, r2, tof, prograde=True, low_path=True)

    # Expected final results
    expected_v1 = np.array([2.058913, 2.915965, 0.0])  # [km / s]
    expected_v2 = np.array([-3.451565, 0.910315, 0.0])  # [km / s]

    # Assert the results
    assert_allclose(v1, expected_v1, atol=ATOL, rtol=RTOL)
    assert_allclose(v2, expected_v2, atol=ATOL, rtol=ATOL)


@pytest.mark.parametrize("solver", ALL_SOLVERS)
def test_case_from_curtiss_book(solver):
    """
    Example 5.2 from Orbital Mechanics for Engineering Students (3rd
    Edition), by Howard D. Curtiss
    """
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([5000.0, 10000.0, 2100.0])  # [km]
    r2 = np.array([-14600.0, 2500.0, 7000.0])  # [km]
    tof = 3600  # [s]

    # Solve the problem
    v1, v2 = solver(mu_earth, r1, r2, tof, prograde=True)

    # Expected final results
    expected_v1 = np.array([-5.9925, 1.9254, 3.2456])  # [km / s]
    expected_v2 = np.array([-3.3125, -4.1966, -0.38529])  # [ km / s ]

    # Assert the results
    assert_allclose(v1, expected_v1, atol=ATOL, rtol=RTOL)
    assert_allclose(v2, expected_v2, atol=ATOL, rtol=RTOL)


@pytest.mark.parametrize("solver", ALL_SOLVERS_LAMBERTHUB)
def test_case_from_battin_book(solver):
    """
    Example 7.12 from An Introduction to the Mathematics and Methods of
    Astrodynamics (Revised Edition), by Richard H. Battin
    """
    # Initial conditions
    mu_sun = 39.47692641  # [AU ** 3 / year ** 2]
    r1 = np.array([0.159321004, 0.579266185, 0.052359607])  # [AU]
    r2 = np.array([0.057594337, 0.605750797, 0.068345246])  # [AU]
    tof = 0.010794065  # [year]

    # Solving the problem
    v1, v2 = solver(mu_sun, r1, r2, tof)

    # Expected final results
    expected_v1 = np.array([-9.303603251, 3.018641330, 1.536362143])  # [AU / year]

    # Assert the results
    assert_allclose(v1, expected_v1, atol=ATOL, rtol=RTOL)


@pytest.mark.parametrize("solver", ALL_SOLVERS)
def test_case_from_GMAT_hyperbolic_prograde(solver):
    """
    This test case was generated using GMAT2020a software from NASA by defining
    an initial position and velocity vectors around the Earth and propagating
    the resulting orbit in time.
    """
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([7100.0, 200.0, 1300.0])  # [km]
    # This vector was obtained after propagating the initial one with original
    # velocity a given amount of time.
    r2 = np.array([-38113.5870, 67274.1946, 29309.5799])  # [km]
    tof = 12000.0  # [s]

    # Solving the problem
    v1, v2 = solver(mu_earth, r1, r2, tof)

    # Expected final results
    expected_v1 = np.array([0.0, 10.35, 5.5])  # [km / s]
    expected_v2 = np.array([-3.6379, 4.4932, 1.7735])  # [km / s]

    # Assert the results
    assert_allclose(v1, expected_v1, atol=ATOL, rtol=RTOL)
    assert_allclose(v2, expected_v2, atol=ATOL, rtol=RTOL)


@pytest.mark.parametrize("solver", ALL_SOLVERS)
def test_case_from_GMAT_hyperbolic_retrograde(solver):
    """
    This test case was generated using GMAT2020a software from NASA by defining
    an initial position and velocity vectors around the Earth and propagating
    the resulting orbit in time.
    """
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([7100.0, 200.0, 1300.0])  # [km]
    # This vector was obtained after propagating the initial one with original
    # velocity a given amount of time.
    r2 = np.array([-47332.7499, -54840.2027, -37100.17067])  # [km]
    tof = 12000.0  # [s]

    # Solving the problem
    v1, v2 = solver(mu_earth, r1, r2, tof, prograde=False)

    # Expected final results
    expected_v1 = np.array([0.0, -10.35, -5.5])  # [km / s]
    expected_v2 = np.array([-4.3016, -3.4314, -2.5467])  # [km / s]

    # Assert the results
    assert_allclose(v1, expected_v1, atol=ATOL, rtol=RTOL)
    assert_allclose(v2, expected_v2, atol=ATOL, rtol=RTOL)


TABLE_OF_TRANSFERS_I = {
    "M0_prograde_low": [
        np.array([-2.09572809, 3.92602196, -4.94516810]),  # [km / s]
        np.array([2.46309613, 0.84490197, 6.10890863]),  # [km / s]
    ],
    "M0_retrograde_high": [
        np.array([1.94312182, -4.35300015, 4.54630439]),  # [km / s]
        np.array([-2.38885563, -1.42519647, -5.95772225]),  # [km / s]
    ],
}


@pytest.mark.parametrize("solver", ALL_SOLVERS)
@pytest.mark.parametrize("case", TABLE_OF_TRANSFERS_I)
def test_case_from_der_article_I(solver, case):
    """
    Example from Astrodynamics 102, by Gim J. Der. see:
    http://derastrodynamics.com/docs/astrodynamics_102_v2.pdf
    """
    # Initial conditions
    mu_earth = 3.986004418e5  # [km ** 3 / s ** 2]
    r1 = np.array([2249.171260, 1898.007100, 5639.599193])  # [km]
    r2 = np.array([1744.495443, -4601.556054, 4043.864391])  # [km]
    tof = 1618.50  # [s]

    # Unpack problem conditions
    M, sense, path = case.split("_")

    # Convert proper type
    M = int(M[-1])
    sense = True if sense == "prograde" else False
    path = True if path == "low" else False

    # Solve the problem
    v1, v2 = solver(mu_earth, r1, r2, tof, M=M, prograde=sense, low_path=path)

    # Expected final results
    expected_v1, expected_v2 = TABLE_OF_TRANSFERS_I[case]

    # Assert the results
    assert_allclose(v1, expected_v1, atol=ATOL, rtol=RTOL)
    assert_allclose(v2, expected_v2, atol=ATOL, rtol=RTOL)


TABLE_OF_TRANSFERS_II = {
    "M0_prograde_high": [
        np.array([2.000652697, 0.387688615, -2.666947760]),  # [km / s]
        np.array([-3.79246619, -1.77707641, 6.856814395]),  # [km / s]
    ],
    "M0_retrograde_high": [
        np.array([2.96616042, -1.27577231, -0.75545632]),  # [km / s]
        np.array([5.8437455, -0.20047673, -5.48615883]),  # [km / s]
    ],
}


@pytest.mark.parametrize("solver", ALL_SOLVERS)
@pytest.mark.parametrize("case", TABLE_OF_TRANSFERS_II)
def test_case_from_der_article_II(solver, case):
    """
    Example from Astrodynamics 102, by Gim J. Der. see:
    http://derastrodynamics.com/docs/astrodynamics_102_v2.pdf
    """
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
    expected_v1, expected_v2 = TABLE_OF_TRANSFERS_II[case]

    # Assert the results
    assert_allclose(v1, expected_v1, atol=ATOL, rtol=RTOL)
    assert_allclose(v2, expected_v2, atol=ATOL, rtol=RTOL)
