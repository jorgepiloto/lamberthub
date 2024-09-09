"""
Lambert's problem solver using the method proposed by Giulio Avanzini in 2008.

This module holds the implementation of the algorithm devised by Avanzini, which
appeared for the first time in his article [1]. The solver takes advantage of
the conservation of the eccentricity projection along the chord direction, so it
iterates over this orbital element till the time of transfer from Kepler's
equation matches the desired transfer one.

However, the original routines was not developed to work under multi-revolutions
conditions. In addition, Avanzini did not provide the derivative of Kepler's
equation with respect to the transverse eccentricity component, forcing a
bounded numerical method (bisection or regulafalsi). This was solved by Quan He,
which proposed new changes to Avanzini's original algorithm in an article [2]
proposed in 2010.

Finally, a new method [3] was proposed by Changxuan Wen in 2014, simplifying all
previous routines.

[1] Avanzini, G. (2008). A simple Lambert algorithm. Journal of guidance,
    control, and dynamics, 31(6), 1587-1594.

[2] He, Q., Li, J., & Han, C. (2010). Multiple-revolution solutions of the
    transverse-eccentricity-based lambert problem. Journal of guidance,
    control, and dynamics, 33(1), 265-269.

[3] Wen, C., Zhao, Y., & Shi, P. (2014). Derivative analysis and algorithm
    modification of transverse-eccentricity-based Lambert problem. Journal of
    Guidance, Control, and Dynamics, 37(4), 1195-1201.

"""

import time

import numpy as np
from scipy.optimize import newton

from lamberthub.ecc_solvers.utils import (
    _f,
    coe_at_eccT,
    get_fundamental_ellipse_properties,
    get_geometry,
)
from lamberthub.utils.assertions import assert_parameters_are_valid
from lamberthub.utils.elements import coe2rv


def avanzini2008(
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
    Solves Lambert problem using Gooding's devised algorithm.

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
    The following routine might be simplified making use of private functions.
    However, we decided to expose all the auxiliary routines to properly
    reproduce original report figures.

    """
    # Check proper value of multi-revolution. Although we could not accept M at
    # all, this ensures all routines within the library work under the same
    # number and type of parameters.
    if M > 0:
        raise ValueError(
            "Avanzini is not able to work within the multi-revolution scenario!"
        )

    # Check that input parameters are safe
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    # Start by computing an auxiliary set of geometric parameters
    geometry = get_geometry(r1, r2, prograde)

    # Filter out the evolution of ecc_T w.r.t. the independent variable based on
    # the transfer angle criteria, see original report flowchart.
    eccT_at_x = _get_eccT_at_x(geometry)

    # Iterative process based on Newton-Raphson secant method begins. Avanzini
    # did not provided a derivative for Kepler's equation of time w.r.t. the ecc_T
    # variable, so this root solver is the one required rather than the pure N-R
    # one.
    tic = time.perf_counter()
    x_sol, r = newton(
        _f,
        0,
        args=(eccT_at_x, mu, geometry, tof),
        maxiter=maxiter,
        tol=atol,
        rtol=rtol,
        full_output=True,
    )
    tac = time.perf_counter()

    # Extract the number of iterations
    numiter = r.iterations

    # Compute the time per iteration
    tpi = (tac - tic) / numiter

    # Solve the actual value of ecc_T at solved x and retrieve COE elements
    ecc_T = eccT_at_x(x_sol)
    p, ecc, inc, raan, argp, nu_1, nu_2 = coe_at_eccT(ecc_T, r1, r2, prograde)

    # Compute the velocity vectors from the classic orbital elements
    (_, v1), (_, v2) = [coe2rv(mu, p, ecc, inc, raan, argp, nu) for nu in [nu_1, nu_2]]

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def _get_eccT_limits(geometry):
    """
    Computes transverse eccentricity value limits as of problem geometry.
    """
    # Solve for the fundamental ellipse properties
    r1_norm, r2_norm, c_norm, dtheta, _ = geometry
    ecc_F, a_F, p_F = get_fundamental_ellipse_properties(r1_norm, r2_norm, c_norm)

    # Compute the limits
    ecc_max = -1 / np.cos(dtheta / 2)
    ecc_H = np.sqrt(ecc_max**2 - ecc_F**2) if dtheta > np.pi else np.inf
    ecc_P = np.sqrt(1 - ecc_F**2)

    return ecc_H, ecc_P


def _get_eccT_at_x(geometry):
    """
    Returns proper transverse eccentricity function depending on the value
    of the transfer angle.

    Parameters
    ----------
    ecc_H: float
        Lower value limit.
    ecc_P: float
        Upper value limit.
    dtheta: float
        Transfer angle.

    Returns
    -------
    eccT_at_x: function
        A function to be evaluated at independent variable values.

    Notes
    -----
    These are equations (16) and (18) from the official report [1].

    """
    # Compute the limits of the ecc_T value
    ecc_H, ecc_P = _get_eccT_limits(geometry)
    _, _, _, dtheta, _ = geometry

    if dtheta > np.pi:
        # Equation (16) is applied
        def eccT_at_x(x):
            X = np.exp((1 / ecc_H + 1 / ecc_P) * x)
            return ecc_P * ecc_H * (X - 1) / (ecc_P + ecc_H * X)

    else:
        # Equation (18) is applied
        def eccT_at_x(x):
            return ecc_P * (1 - np.exp(-x / ecc_P))

    return eccT_at_x
