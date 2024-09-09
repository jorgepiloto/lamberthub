"""This module holds all methods devised by Carl Fredrich Gauss."""

import time

import numpy as np
from numpy.linalg import norm

from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import (
    assert_parameters_are_valid,
    assert_transfer_angle_not_pi,
    assert_transfer_angle_not_zero,
)


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

    # Compute the s and w constants
    s = _get_s(r1_norm, r2_norm, dtheta)
    w = _get_w(mu, tof, r1_norm, r2_norm, dtheta)

    # Initial guess formulation is of the arbitrary type
    y0 = 1.00

    # The iterative procedure can start now
    tic = time.perf_counter()
    for numiter in range(1, maxiter + 1):
        # Compute the value of the free-parameter
        x = _gauss_first_equation(y0, s, w)

        # Evaluate the new value of the dependent variable
        y = _gauss_second_equation(x, s)

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


def _get_s(r1_norm, r2_norm, dtheta):
    """Returns the s auxiliary constant.

    Parameters
    ----------
    r1_norm: float
        The norm of the initial position vector.
    r2_norm: float
        The norm of the final position vector.
    dtheta: float
        The transfer angle.

    Returns
    -------
    s: float
        An auxiliary constant.

    Notes
    -----
    This is equation number (5.6-2) from Bate's book [2].

    """
    s = (r1_norm + r2_norm) / (
        4 * np.sqrt(r1_norm * r2_norm) * np.cos(dtheta / 2)
    ) - 1 / 2
    return s


def _get_w(mu, tof, r1_norm, r2_norm, dtheta):
    """Returns the w auxiliary constant.

    Parameters
    ----------
    mu: float
        The gravitational constant.
    tof: float
        The time of flight.
    r1_norm: float
        The norm of the initial position vector.
    r2_norm: float
        The norm of the final position vector.
    dtheta: float
        The transfer angle.

    Returns
    -------
    w: float
        An auxiliary constant.

    Notes
    -----
    This is equation number (5.6-3) from Bate's book [2].

    """
    w = (mu * tof**2) / (2 * np.sqrt(r1_norm * r2_norm) * np.cos(dtheta / 2)) ** 3
    return w


def _gauss_first_equation(y, s, w):
    """Evaluates Gauss' first equation.

    Parameters
    ----------
    y: float
        The dependent variable.
    s: float
        First auxiliary variable.
    w: float
        Second auxiliary variable.

    Returns
    -------
    x: float
        The independent variable or free-parameter.

    Notes
    -----
    This is equation (5.6-13) from Bate's book [2].

    """
    x = w / y**2 - s
    return x


def _gauss_second_equation(x, s):
    """Evaluates Gauss' second equation.

    Parameters
    ----------
    x: float
        The independent variable.
    s: float
        First auxiliary variable

    Returns
    -------
    y: float
        Dependent variable.

    Notes
    -----
    This is equation (5.6-14) from Bate's book, reference [2].

    """
    y = 1 + _X_at_x(x) * (s + x)
    return y


def _X_at_x(x, order=50):
    """Computes capital X as function of lower x.

    Parameters
    ----------
    x: float
        The independent variable.
    n_items: int
        The number of terms to be considered in the series.

    Returns
    -------
    X: float
        The series summation.

    Notes
    -----
    This is equation (5.6-15) from Bate's book, in reference [2].

    """
    coefficients = [1]
    for n in range(3, (3 + order)):
        coeff = (2 * n) / (2 * n - 1)
        coefficients.append(np.prod(coefficients[-1]) * coeff * x)
    X = (4 / 3) * np.sum(coefficients)
    return X
