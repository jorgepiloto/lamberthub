"""This module holds all methods devised by Carl Fredrich Gauss."""

import time

import numpy as np
from numpy.linalg import norm

from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import (assert_parameters_are_valid,
                                         assert_transfer_angle_not_pi,
                                         assert_transfer_angle_not_zero)


def gauss1809(
    mu,
    r1,
    r2,
    tof,
    M=0,
    prograde=True,
    low_path=True,
    maxiter=250,
    atol=1e-5,
    rtol=1e-7,
    full_output=False,
):
    r"""
    Lambert's problem solver devised by Carl Friedrich Gauss in 1809. The method
    has been implemented according to Bate's book (see [2]) and extended to the
    hyperbolic case. This method shows poor accuracy, being only suitable for
    low transfer angles.

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
    The algorithm originally devised by Gauss exploits the so-called ratio of
    sector to triangle area, which is a numerical value related with the orbital
    parameter. This Algorithm was used to the discovery of the orbit of Ceres by
    the genius and adopted by many other authors of his time due to its
    simplicity. However, the Algorithm is found to be singular for transfer
    angles of 180 degrees and shows a low performance for really small angles.

    References
    ----------
    [1] Gauss, C. F. (1809). Theoria motus corporum coelestium in sectionibus
    conicis solem ambientium auctore Carolo Friderico Gauss. sumtibus Frid.
    Perthes et IH Besser.

    [2] Bate, R. R., Mueller, D. D., White, J. E., & Saylor, W. W. (2020).
    Fundamentals of astrodynamics. Courier Dover Publications.

    """

    # Check that input parameters are safe
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    # Norm of the initial and final position vectors
    r1_norm, r2_norm = [norm(r) for r in [r1, r2]]

    # Compute the cosine of the transfer angle and check
    dtheta = get_transfer_angle(r1, r2, prograde)
    [
        check_angle(dtheta)
        for check_angle in [
            assert_transfer_angle_not_zero,
            assert_transfer_angle_not_pi,
        ]
    ]

    # Obtain the constants l and m
    ll = (r1_norm + r2_norm) / (
        4 * np.sqrt(r1_norm * r2_norm) * np.cos(dtheta / 2)
    ) - 1 / 2
    m = (mu * tof ** 2) / (2 * np.sqrt(r1_norm * r2_norm) * np.cos(dtheta / 2)) ** 3

    # Initial guess formulation is of the arbitrary type
    y0 = 1.00

    # The iterative procedure can start now
    tic = time.perf_counter()
    for numiter in range(1, maxiter + 1):

        # Find for the value of capital X
        x = m / y0 ** 2 - ll
        X = _X_at_x(x)

        # Use previously computed value of X for a better approximation of y
        y = 1 + X * (ll + x)

        # Check the convergence of the method
        if np.abs(y - y0) <= atol:
            tac = time.perf_counter()
            tpi = (tac - tic) / numiter
            break
        else:
            # The new initial guess is the previously computed y value
            y0 = y
    else:
        raise ValueError("Exceeded maximum number of iterations.")

    # Once the value of y has been found, the shape and nature of the orbit can
    # be determined from the sign of x variable
    if x > 0:
        # The orbit is an ellipse, being deltaAnomaly = deltaE
        deltaAnomaly = 2 * np.arccos(1 - 2 * x)
    elif x == 0:
        # The orbit is a parabola
        pass
    else:
        # The orbit is an hyperbola, being deltaAnomaly = deltaF
        deltaAnomaly = 2 * np.arccosh(1 - 2 * x)

    # Compute the orbital parameter
    p = (r1_norm * r2_norm * (1 - np.cos(dtheta))) / (
        r1_norm
        + r2_norm
        - 2 * np.sqrt(r1_norm * r2_norm) * np.cos(dtheta / 2) * np.cos(deltaAnomaly / 2)
    )

    # Evaluate f, g, f_dot and g_dot functions for future solving v1 and v2
    f = 1 - r2_norm * (1 - np.cos(dtheta)) / p
    g = (r1_norm * r2_norm * np.sin(dtheta)) / np.sqrt(mu * p)
    f_dot = (
        np.sqrt(mu / p)
        * np.tan(dtheta / 2)
        * ((1 - np.cos(dtheta)) / p - 1 / r1_norm - 1 / r1_norm)
    )
    g_dot = 1 - r1_norm * (1 - np.cos(dtheta)) / p

    # Compute the initial and final velocity vectors
    v1 = (r2 - f * r1) / g
    v2 = f_dot * r1 + g_dot * v1

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def _X_at_x(x):
    """Computes capital X as function of lower x.

    Parameters
    ----------
    x: float

    Returns
    -------
    X: float

    """
    # The first two coefficients are stored so the list can be updated
    # dynamically in the following statement
    coeff = [1, 6 / 5 * x]
    # Compute only up to the tenth power. Gauss' method is a non-robust one and
    # computing additional series elements does not change its nature. In fact,
    # it is seen to be detrimental from the point of view of memory.
    [
        coeff.append(coeff[i - 1] * (2 * i + 4) / (2 * i + 3) * x ** i)
        for i in range(2, 10)
    ]
    X = 4 / 3 * sum(coeff)
    return X
