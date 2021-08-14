"""This module holds all methods devised by David A. Vallado."""

import time

import numpy as np

from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import assert_parameters_are_valid
from lamberthub.utils.stumpff import c2, c3


def vallado2013(
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
    Vallado's algorithm makes use of the universal formulation to solve for the
    Lambert's problem. By making use of a bisection method, it guarantees the
    convergence to the solution but the amount of iterations require
    dramatically increases.

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

    Notes
    -----
    This algorithm is presented as an alternative to the one developed by Bate
    in 1971. Bate did not impose a particular numerical solver for his algorithm
    but cited both bisection and Newton's one. However, for some values of the
    boundary problem, the initial guess might diverge if Newton's solver is
    used. That's why Vallado decided to employ a bisection method instead.
    Although detrimental from the point of view of performance, this algorithm
    properly reaches solution in the majority of the cases.

    All credits of the implementation go to Juan Luis Cano Rodríguez and the
    poliastro development team, from which this routine inherits. Some changes
    were made to adapt it to `lamberthub` API.

    Copyright (c) 2012-2021 Juan Luis Cano Rodríguez and the poliastro
    development team.

    References
    ----------
    [1] Vallado, D. A. (2001). Fundamentals of astrodynamics and applications
    (Vol. 12). Springer Science & Business Media.

    """

    # Check that input parameters are safe
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    # Retrieve the fundamental geometry of the problem
    r1_norm, r2_norm, c_norm = [np.linalg.norm(vec) for vec in [r1, r2, r2 - r1]]
    dtheta = get_transfer_angle(r1, r2, prograde)

    t_m = 1 if dtheta < np.pi else -1

    norm_r1_times_norm_r2 = r1_norm * r2_norm
    norm_r1_plus_norm_r2 = r1_norm + r2_norm

    cos_dnu = np.dot(r1, r2) / norm_r1_times_norm_r2

    A = t_m * (r1_norm * r2_norm * (1 + cos_dnu)) ** 0.5

    if A == 0.0:
        raise RuntimeError("Cannot compute orbit, phase angle is 180 degrees")

    psi = 0.0
    psi_low = -4 * np.pi ** 2
    psi_up = 4 * np.pi ** 2

    tic = time.perf_counter()
    for numiter in range(1, maxiter + 1):
        y = norm_r1_plus_norm_r2 + A * (psi * c3(psi) - 1) / c2(psi) ** 0.5
        if A > 0.0:
            # Readjust xi_low until y > 0.0
            # Translated directly from Vallado
            while y < 0.0:
                psi_low = psi
                psi = (
                    0.8
                    * (1.0 / c3(psi))
                    * (1.0 - norm_r1_times_norm_r2 * np.sqrt(c2(psi)) / A)
                )
                y = norm_r1_plus_norm_r2 + A * (psi * c3(psi) - 1) / c2(psi) ** 0.5

        xi = np.sqrt(y / c2(psi))
        tof_new = (xi ** 3 * c3(psi) + A * np.sqrt(y)) / np.sqrt(mu)

        # Convergence check
        if np.abs((tof_new - tof) / tof) < rtol:
            tac = time.perf_counter()
            tpi = (tac - tic) / numiter
            break

        # Bisection check
        condition = tof_new <= tof
        psi_low = psi_low + (psi - psi_low) * condition
        psi_up = psi_up + (psi - psi_up) * (not condition)

        psi = (psi_up + psi_low) / 2
    else:
        raise ValueError("Exceeded maximum number of iterations!")

    f = 1 - y / r1_norm
    g = A * np.sqrt(y / mu)

    gdot = 1 - y / r2_norm

    v1 = (r2 - f * r1) / g
    v2 = (gdot * r2 - r1) / g

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)
