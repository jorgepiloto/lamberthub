"""This module holds all methods devised by Ruiye Jiang."""

import time

import numpy as np
from scipy.optimize import bisect, newton

from lamberthub.a_solvers.utils import _get_alpha0, _get_beta0, _get_orbit_parameter
from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import assert_parameters_are_valid


def jiang2016(
    mu,
    r1,
    r2,
    tof,
    M=0,
    prograde=True,
    low_path=True,
    maxiter=100,
    atol=1e-5,
    rtol=1e-7,
    full_output=False,
):
    r"""
    Battin's elegant algortihm for solving the Lambert's problem. This algorithm
    is known to improve Gauss original one by removing the singularity for 180
    transfer angles and increasing its performance.

    Parameters
    ----------
    mu: float
        Gravitational parameter, equivalent to :math:`GM` of attractor body.
    r1: numpy.array
        Initial position vector.
    r2: numpy.array
        Final position vector.
    M: int
        Number of revolutions. Must be equal or greater than 0 value.
    prograde: bool
        If `True`, specifies prograde motion. Otherwise, retrograde motion is imposed.
    low_path: bool
        If two solutions are available, it selects between high or low path.
    maxiter: int
        Maximum number of iterations.
    atol: float
        Absolute tolerance.
    rtol: float
        Relative tolerance.
    full_output: bool
        If True, the number of iterations and time per iteration are also returned.

    Returns
    -------
    v1: numpy.array
        Initial velocity vector.
    v2: numpy.array
        Final velocity vector.
    numiter: int
        Number of iterations.
    tpi: float
        Time per iteration in seconds.

    """

    # Check that input parameters are safe
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    # Retrieve the fundamental geometry of the problem
    r1_norm, r2_norm, c_norm = [np.linalg.norm(vec) for vec in [r1, r2, r2 - r1]]
    semiperimeter = (r1_norm + r2_norm + c_norm) / 2
    dtheta = get_transfer_angle(r1, r2, prograde)

    # Solve the minimum energy orbit
    a_min = semiperimeter / 2

    # Compute the time of flight for a parabolic transfer
    tof_parabolic = _lagrange_tof_parabolic(mu, r1_norm, r2_norm, c_norm)

    # Check if the current time of flight is greater (elliptic) or lower
    # (hyperbolic)
    if tof > tof_parabolic:
        a_lower, a_upper = _initial_guess_elliptic(
            a_min, mu, tof, dtheta, c_norm, semiperimeter, maxiter
        )
        _f_tof = _f_elliptic
    else:
        a_lower, a_upper = _initial_guess_hyperbolic(
            a_min, mu, tof, dtheta, c_norm, semiperimeter, maxiter
        )
        _f_tof = _f_hyperbolic

    # Apply the bisection method
    a = bisect(_f_tof, a_lower, a_upper, args=(mu, tof, dtheta, c_norm, semiperimeter))

    # Compute the value of alpha
    if tof > tof_parabolic:
        alpha = 2 * np.pi - _get_alpha0(a, semiperimeter)
    else:
        alpha = _get_alpha0(a, semiperimeter)

    # Compute the value of beta 
    if a > 0:
        beta = (
            _get_beta0(a, c_norm, semiperimeter)
            if dtheta < np.pi
            else -_get_beta0(a, c_norm, semiperimeter)
        )
    else:
        beta = _get_beta0

    # Get the orbit parameter
    p = _get_orbit_parameter(a, r1_norm, r2_norm, c_norm, semiperimeter, alpha, beta)

    # Assembly velocity vectors
    v1 = (
        np.sqrt(mu * p)
        / (r1_norm * r2_norm * np.sin(dtheta))
        * ((r2 - r1) + r2_norm / p * (1 - np.cos(dtheta)) * r1)
    )
    v2 = (
        np.sqrt(mu * p)
        / (r1_norm * r2_norm * np.sin(dtheta))
        * ((r2 - r1) - r1_norm / p * (1 - np.cos(dtheta)) * r2)
    )

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def _initial_guess_hyperbolic(a_min, mu, tof, dtheta, c, s, maxiter):
    a0, a1 = -(10 ** 6), -(10 ** 8)

    for _ in range(maxiter):
        f_a0 = _f_hyperbolic(a0, mu, tof, dtheta, c, s)
        f_a1 = _f_hyperbolic(a1, mu, tof, dtheta, c, s)

        if f_a0 * f_a1 < 0:
            break
        
        a0, a1 = a1, 2 * a1

    return a0, a1


def _initial_guess_elliptic(a_min, mu, tof, dtheta, c, s, maxiter):
    a0, a1 = a_min, 10 * a_min

    for _ in range(maxiter):
        f_a0 = _f_elliptic(a0, mu, tof, dtheta, c, s)
        f_a1 = _f_elliptic(a1, mu, tof, dtheta, c, s)

        if f_a0 * f_a1 < 0:
            break

        a0, a1 = a1, 2 * a1

    return a0, a1


def _f_elliptic(a, mu, tof, dtheta, c, s):
    alpha = 2 * np.pi - _get_alpha0(a, s)
    beta = _get_beta0(a, c, s) if dtheta <= np.pi else -_get_beta0(a, c, s)
    return tof - _lagrange_tof_elliptic(a, mu, alpha, beta)


def _f_hyperbolic(a, mu, tof, dtheta, c, s):
    alpha = _get_alpha0(a, s)
    beta = _get_beta0(a, c, s) if dtheta <= np.pi else -_get_beta0(a, c, s)
    return tof - _lagrange_tof_hyperbolic(a, mu, alpha, beta)


def _lagrange_tof_elliptic(a, mu, alpha, beta):
    tof = (a ** (3 / 2) * ((alpha - np.sin(alpha)) - (beta - np.sin(beta)))) / np.sqrt(
        mu
    )
    return tof


def _lagrange_tof_parabolic(mu, r1_norm, r2_norm, c_norm):
    tof = (
        (1 / 6) * (r1_norm + r2_norm + c_norm) ** (3 / 2)
        - (1 / 6) * (r1_norm + r2_norm - c_norm) ** (3 / 2)
    ) / np.sqrt(mu)
    return tof


def _lagrange_tof_hyperbolic(a, mu, alpha, beta):
    tof = (
        (-a) ** (3 / 2) * ((np.sinh(alpha) - alpha) - (np.sinh(beta) - beta))
    ) / np.sqrt(mu)
    return tof
