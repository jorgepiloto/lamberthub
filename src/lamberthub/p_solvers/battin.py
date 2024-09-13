"""This module holds all methods devised by R. H. Battin."""

import time

import numpy as np

from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import assert_parameters_are_valid


def battin1984(
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
    Battin's elegant algorithm for solving the Lambert's problem. This algorithm
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
    [1] Battin, R. H., & Vaughan, R. M. (1984). An elegant Lambert algorithm.
    Journal of Guidance, Control, and Dynamics, 7(6), 662-670.

    [2] Battin, R. H. (1999). An introduction to the mathematics and methods of
    astrodynamics. Aiaa.

    [3] Vaughan, R. M. (1983). An improvement of Gauss' method for solving
    Lambert's problem (Doctoral dissertation, Massachusetts Institute of
    Technology).

    """
    # Check that input parameters are safe
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    # Retrieve the fundamental geometry of the problem
    r1_norm, r2_norm, c_norm = [np.linalg.norm(vec) for vec in [r1, r2, r2 - r1]]
    semiperimeter = (r1_norm + r2_norm + c_norm) / 2
    dtheta = get_transfer_angle(r1, r2, prograde)

    # Compute the auxiliary lambda parameter. The variable _lambda is used here
    # as lambda is a reserved keyword in Python for declaring inline functions.
    # This expression is equivalent to equation (1) from [1]. Equation (7.122a)
    # applies here. The sign is selected according to the value of the transfer
    # angle, being positive if lower than 180 degrees or negative if greater.
    # Therefore, Battin's lambda variable is the so-called transfer angle
    # parameter.
    _lambda = _get_lambda(c_norm, semiperimeter, dtheta)

    # Solve for the new auxiliary variables developed by Battin. Equations (30)
    # and (31) from the report [1] apply here.
    ll = _get_ll(_lambda)
    m = _get_m(mu, tof, semiperimeter, _lambda)

    # The non-dimensional transfer time is computed using equation (54) from
    # Battin's report [1]. The non-dimensional time for the parabolic case is
    # also computed here using equation appearing in page 340 of Battin's book
    # [2]. The author did not provide this in the original paper [1], which is a
    # critical errata, as it is used in the initial guess computation.
    T = np.sqrt(8 * mu / semiperimeter**3) * tof
    T_p = (4 / 3) * (1 - _lambda**3)

    # The initial guess procedure is set according to piece-wise no-numbered
    # equation appearing in the last page of the report [1]. This relation is
    # also shown in page 340 from book [2], where it is claimed that if T is
    # greater than T_p, the transfer orbit is expected to be elliptic type,
    # otherwise parabolic or hyperbolic are assumed.
    x0 = ll if T > T_p else 0

    # The iterative procedure starts
    tic = time.perf_counter()
    for numiter in range(1, maxiter + 1):
        # Evaluate the h coefficients given by equations (47) and (48) from the
        # report [1] or the relations (7.111) and (7.112) from the book [2]. The
        # denominator is common for both h1 and h2 coefficients, so it is
        # to avoid computing it twice
        h1, h2 = _get_h_coefficients(x0, ll, m)

        # Compute the auxiliary variable u, which is used in the computation of
        # the real positive root in Battin's second equation.
        u = _u_at_h(h1, h2)

        y = _battin_second_equation(u, h1, h2)

        # Compute the new value for the free-parameter
        x = _battin_first_equation(y, ll, m)

        # Check if the new computed value lies within desired tolerance. If so,
        # stop the iteration procedure, otherwise update the initial guess and
        # keep performing the computation.
        if np.abs(x - x0) <= atol:
            tac = time.perf_counter()
            tpi = (tac - tic) / numiter
            break
        else:
            x0 = x
    else:
        raise ValueError("Exceeded maximum number of iterations!")

    r11 = (1 + _lambda) ** 2 / (4 * tof * _lambda)
    s11 = y * (1 + x)
    t11 = (m * semiperimeter * (1 + _lambda) ** 2) / s11

    # Compute the initial and final velocity vectors
    v1 = -r11 * (s11 * (r1 - r2) - t11 * r1 / r1_norm)
    v2 = -r11 * (s11 * (r1 - r2) + t11 * r2 / r2_norm)

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def _battin_first_equation(y, ll, m):
    """Battin's first equation.

    Parameters
    ----------
    y: float
        The dependent variable.
    ll: float
        First auxiliary variable.
    m: float
        Second auxiliary variable.

    Returns
    -------
    x: float
        The independent variable.

    Notes
    -----
    This is equaiton (49) from original report [1], which is alternative from
    the book [2] is expression (7.113)

    """
    x = np.sqrt(((1 - ll) / 2) ** 2 + m / y**2) - (1 + ll) / 2
    return x


def _battin_second_equation(u, h1, h2):
    """Battin's second equation.

    Parameters
    ----------
    u: float
        The dependent variable.
    h1: float
        The first of the h coefficients.
    h2: float
        The second of the h coefficients.

    Returns
    -------
    y: float
        The dependent variable.

    Notes
    -----
    Computes the desired positive root via an explicit equation defined in
    [1] under number (57).
    """
    # Evaluate the auxiliary parameter B and the K function
    B, K = _B_at_h(h1, h2), _K_at_u(u)
    y = ((1 + h1) / 3) * (2 + np.sqrt(B + 1) / (1 - 2 * u * K))
    return y


def _get_lambda(c, s, dtheta):
    """Compute the transfer angle parameter.

    Parameters
    ----------
    c: float
        The norm of the chord vector.
    s: float
        The semiperimeter.
    dtheta: float
        The transfer angle in radians.

    Returns
    -------
    _lambda: float
        The transfer angle parameter.

    Notes
    -----
    This is inline equation in the very first page of report [1].

    """
    _lambda = np.sqrt(s * (s - c)) / s
    _lambda = np.abs(_lambda) if dtheta < np.pi else -np.abs(_lambda)
    return _lambda


def _get_ll(_lambda):
    """Computes the l variable.

    Parameters
    ----------
    _lambda: float
        The transfer angle parameter.

    Returns
    -------
    ll: float
        Auxiliary variable.

    Notes
    -----
    This is equation (30) from original report.

    """
    ll = ((1 - _lambda) / (1 + _lambda)) ** 2
    return ll


def _get_m(mu, tof, s, _lambda):
    """Computes the m auxiliary variable.

    Parameters
    ----------
    mu: float
        The gravitational parameter.
    tof: float
        The time of flight.
    s: float
        The semiperimeter of the orbit.
    _lambda: float
        The transfer angle parameter.

    Returns
    -------
    m: float
        Auxiliary variable.

    Notes
    -----
    This is equation (31) from official report [1].

    """
    m = (8 * mu * tof**2) / (s**3 * (1 + _lambda) ** 6)
    return m


def _get_h_coefficients(x, ll, m):
    """Evaluates the h1 and h2 coefficients.

    Parameters
    ----------
    x: float
        The free-parameter.
    ll: float
        The first auxiliary variable.
    m: float
        The second auxiliary variable.

    Returns
    -------
    h1: float
        The first of the h coefficients.
    h2: float
        The second of the h coefficients.

    Notes
    -----
    These are equations (47) and (48) from report [1].

    """
    # Evaluate the xi function at particular value of x. This is equation
    # (53) from original report [1] or expression (7.121) from the book [2].
    # The function follows the series technique for an easy implementation
    # from the computational point of view.
    xi = _xi_at_x(x)

    # Compute the value of the h coefficients
    h_denominator = (1 + 2 * x + ll) * (4 * x + xi * (3 + x))
    h1 = ((ll + x) ** 2 * (1 + 3 * x + xi)) / h_denominator
    h2 = (m * (x - ll + xi)) / (h_denominator)
    return (h1, h2)


def _u_at_h(h1, h2):
    """Evaluates u at h coefficients.

    Parameters
    ----------
    h1: float
        The first of the h coefficients.
    h2: float
        The second of the h coefficients.

    Returns
    -------
    u: float
        Auxiliary variable.

    """
    u = _u_at_B(_B_at_h(h1, h2))
    return u


def _u_at_B(B):
    """Evaluates u auxiliary variable at given B.

    Parameters
    ----------
    B: float
        Auxiliary variable.

    Returns
    -------
    u: float
        Auxiliary variable.

    Notes
    -----
    This is equation (55) from the original report [1].

    """
    u = -B / (2 * np.sqrt(1 + B) + 1)
    return u


def _B_at_h(h1, h2):
    """Evaluates B auxiliary variable at given h coefficients.

    Parameters
    ----------
    h1: float
        The first of the h coefficients.
    h2: float
        The second of the h coefficients.

    Returns
    -------
    B: float
        Auxiliary variable.

    Notes
    -----
    This is equation (56) from the original report.

    """
    B = (27 * h2) / (4 * (1 + h1) ** 3)
    return B


def _xi_at_x(x, levels=125):
    """Evaluates the xi function at a particular value of x.

    Parameters
    ----------
    x: float
        The independent variable.

    Returns
    -------
    xi: float
        The value of the xi function.

    Notes
    -----
    This is equation (53) from original report [1]. However, the method
    presented in [3] which makes use of a series is the one used here as it is
    much simpler to implement.

    """
    # Compute the value of eta, given by equation (52) in [1].
    eta = x / (np.sqrt(1 + x) + 1) ** 2

    # Initial values
    (
        delta,
        u,
        sigma,
        m1,
    ) = (
        0,
        1,
        1,
        1,
    )

    while np.abs(u) > 1e-18 and m1 <= levels:
        m1 += 1
        gamma = (m1 + 3) ** 2 / (4 * (m1 + 3) ** 2 - 1)
        delta = 1 / (1 + gamma * eta * delta)
        u = u * (delta - 1)
        sigma = sigma + u

    xi = 8 * (np.sqrt(1 + x) + 1) / (3 + 1 / (5 + eta + (9 * eta / 7) * sigma))
    return xi


def _K_at_u(u, levels=1000):
    """Evaluates the K function at a particular value of u.

    Parameters
    ----------
    u: float
        Battin's auxiliary variable.

    Returns
    -------
    K: float
        The value of the K function.

    Notes
    -----
    This is equation (58) from original report [1]. However, the method
    presented in [3] which makes use of a series is the one used here as it is
    much simpler to implement.

    """
    # Initial values
    (
        delta,
        u0,
        sigma,
        n1,
    ) = (
        1,
        1,
        1,
        0,
    )

    while np.abs(u0) > 1e-18 and n1 <= levels:
        if n1 == 0:
            gamma = 4 / 27
            delta = 1 / (1 - gamma * u * delta)
            u0 = u0 * (delta - 1)
            sigma = sigma + u0
        else:
            for val in [1, 2]:
                if val == 1:
                    gamma = (
                        2
                        * (3 * n1 + 1)
                        * (6 * n1 - 1)
                        / (9 * (4 * n1 - 1) * (4 * n1 + 1))
                    )
                else:
                    gamma = (
                        2
                        * (3 * n1 + 2)
                        * (6 * n1 + 1)
                        / (9 * (4 * n1 + 1) * (4 * n1 + 3))
                    )
                delta = 1 / (1 - gamma * u * delta)
                u0 = u0 * (delta - 1)
                sigma = sigma + u0

        n1 = n1 + 1

    K = (sigma / 3) ** 2
    return K
