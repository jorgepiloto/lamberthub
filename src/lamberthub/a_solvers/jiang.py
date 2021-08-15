"""This module holds all methods devised by Ruiye Jiang."""

import time

import numpy as np
from scipy.optimize import bisect, newton

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
    maxiter=35,
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

    # Compute the parabolic transfer time using equation (3) from the report [1]
    mp = -1 if dtheta < np.pi else 1
    tof_p = (
        1 / 6 * (r1_norm + r2_norm + c_norm) ** (3 / 2)
        + mp * 1 / 6 * (r1_norm + r2_norm - c_norm) ** (3 / 2)
    ) / np.sqrt(mu)

    # Compute the semi-major axis of the minimum energy orbit. This is used in
    # the computation of the initial guess.
    a_m = semiperimeter / 2

    # The current time of flight is compared against the parabolic one. If it is
    # greater or lower it means that the transfer orbit is elliptic or
    # hyperbolic respectively.
    if tof > tof_p:
        print(f"Found elliptic orbit")

        a = newton(
            _kepler_elliptic,
            a_m,
            args=(mu, tof, tof_p, dtheta, c_norm, semiperimeter),
            tol=atol,
            maxiter=maxiter,
        )

        (alpha, beta) = _alpha_beta_elliptic(a, dtheta, tof, tof_p, c_norm, semiperimeter)
        p = (
            (4 * a * (semiperimeter - r1_norm) * (semiperimeter - r2_norm))
            / (c_norm)
            * np.sin((alpha + beta) / 2) ** 2
        )

    elif tof == tof_p:
        a = np.inf
    else:
        pass

    # Compute the velocity vectors
    v1 = (
        np.sqrt(mu * p)
        / (r1_norm * r2_norm * np.sin(dtheta))
        * ((r2 - r1) + r2_norm / p * (1 - np.cos(dtheta)) * r1)
    )
    v1 = (
        np.sqrt(mu * p)
        / (r1_norm * r2_norm * np.sin(dtheta))
        * ((r2 - r1) - r1_norm / p * (1 - np.cos(dtheta)) * r2)
    )
    return (v1, np.zeros(3))


def _alpha_beta_elliptic(a, tof, tof_p, dtheta, c, s):

    # Compute alpha
    if tof < tof_p:
        alpha = 2 * np.arcsin((s / (2 * a)) ** (1 / 2))
    else:
        alpha = 2 * (np.pi - np.arcsin((s / (2 * a)) ** (1 / 2)))

    # Compute beta
    if dtheta < np.pi:
        beta = 2 * np.arcsin(((s - c) / (2 * a)) ** (1 / 2))
    else:
        beta = -2 * np.arcsin(((s - c) / (2 * a)) ** (1 / 2))

    return (alpha, beta)


def _kepler_elliptic(a, mu, tof, tof_p, dtheta, c, s):
    (alpha, beta) = _alpha_beta_elliptic(a, tof, tof_p, dtheta, c, s)

    return tof - a ** (3 / 2) / np.sqrt(mu) * (
        (alpha - np.sin(alpha)) - (beta - np.sin(beta))
    )
