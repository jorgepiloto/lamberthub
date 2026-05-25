"""Lambert's problem solver using the method proposed by Jiang et al. in 2016."""

import time

import numpy as np
from numpy.linalg import norm

from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import (
    assert_parameters_are_valid,
    assert_transfer_angle_not_pi,
    assert_transfer_angle_not_zero,
)


def jiang2016(
    mu,
    r1,
    r2,
    tof,
    M=0,
    is_prograde=True,
    is_low_path=True,
    maxiter=35,
    atol=1e-5,
    rtol=1e-7,
    full_output=False,
):
    r"""
    Solve Lambert's problem using Jiang's semi-major-axis iteration.

    Parameters
    ----------
    mu: float
        Gravitational parameter, equivalent to :math:`GM` of attractor body.
    r1: numpy.array
        Initial position vector.
    r2: numpy.array
        Final position vector.
    tof: float
        Time of flight between the initial and final position vectors.
    M: int
        Number of revolutions. Must be equal or greater than 0.
    is_prograde: bool
        If `True`, specifies prograde motion. Otherwise, retrograde motion is imposed.
    is_low_path: bool
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
    Jiang, Chao, Wang, and Yang classify the transfer with the parabolic time of
    flight, solve the elliptic or hyperbolic Lagrange time equation by iterating
    on the semi-major axis, and reconstruct the boundary velocities from the
    resulting conic parameter.

    References
    ----------
    Jiang, R., Chao, T., Wang, S., & Yang, M. (2016). Improved semi-major Axis
    iterated method for Lambert's problem. 2016 IEEE Chinese Guidance,
    Navigation and Control Conference, 1423-1428.

    """
    if M > 0:
        raise ValueError(
            "Jiang is not able to work within the multi-revolution scenario!"
        )

    assert_parameters_are_valid(mu, r1, r2, tof, M)

    theta = get_transfer_angle(r1, r2, is_prograde)
    assert_transfer_angle_not_zero(theta)
    assert_transfer_angle_not_pi(theta)

    r1_norm, r2_norm = norm(r1), norm(r2)
    c = norm(r2 - r1)
    s = (r1_norm + r2_norm + c) / 2

    tof_parabolic = _tof_parabolic(mu, s, c, theta)

    tic = time.perf_counter() if full_output else 0.0
    if np.isclose(tof, tof_parabolic, rtol=rtol, atol=atol):
        p = _parabolic_p(r1_norm, r2_norm, c, s, theta)
        numiter = 0
    elif tof > tof_parabolic:
        a_min = s / 2
        tof_min_energy = _tof_elliptic(mu, a_min, s, c, theta, long_period=False)
        long_period = tof > tof_min_energy
        lower, upper = _get_elliptic_bracket(mu, tof, a_min, s, c, theta, long_period)
        a, numiter = _solve_semimajor_axis(
            lambda axis: _tof_elliptic(mu, axis, s, c, theta, long_period) - tof,
            lower,
            upper,
            maxiter,
            atol,
            rtol,
        )
        p = _elliptic_p(a, r1_norm, r2_norm, c, s, theta, long_period)
    else:
        lower, upper = _get_hyperbolic_bracket(mu, tof, s, c, theta)
        a, numiter = _solve_semimajor_axis(
            lambda axis: _tof_hyperbolic(mu, axis, s, c, theta) - tof,
            lower,
            upper,
            maxiter,
            atol,
            rtol,
        )
        p = _hyperbolic_p(a, r1_norm, r2_norm, c, s, theta)

    if full_output:
        tac = time.perf_counter()
        tpi = 0.0 if numiter == 0 else (tac - tic) / numiter
    else:
        tpi = 0.0

    v1, v2 = _reconstruct_velocities(mu, r1, r2, r1_norm, r2_norm, theta, p)
    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def _tof_parabolic(mu, s, c, theta):
    """Compute Jiang's parabolic branch time from equation (3).

    Parameters
    ----------
    mu: float
        Gravitational parameter.
    s: float
        Semiperimeter of the triangle formed by the two position vectors and
        the chord.
    c: float
        Chord length between the initial and final position vectors.
    theta: float
        Transfer angle in radians.

    Returns
    -------
    tof: float
        Time of flight for the limiting parabolic transfer.

    Notes
    -----
    The parabolic time is used by Jiang et al. to classify the transfer as
    elliptic, hyperbolic, or parabolic before solving for the semi-major axis.

    """
    beta_sign = _beta_sign(theta)
    return np.sqrt(2 / mu) * (s**1.5 - beta_sign * (s - c) ** 1.5) / 3


def _tof_elliptic(mu, a, s, c, theta, long_period):
    """Evaluate the elliptic Lagrange time equation.

    Parameters
    ----------
    mu: float
        Gravitational parameter.
    a: float
        Positive semi-major axis of the elliptic transfer orbit.
    s: float
        Semiperimeter of the transfer triangle.
    c: float
        Chord length between the terminal position vectors.
    theta: float
        Transfer angle in radians.
    long_period: bool
        If `True`, use the long-period alpha branch.

    Returns
    -------
    tof: float
        Elliptic transfer time associated with ``a``.

    Notes
    -----
    This is equation (2) of Jiang et al., with alpha and beta computed from
    equation (6).

    """
    alpha, beta = _elliptic_alpha_beta(a, s, c, theta, long_period)
    return np.sqrt(a**3 / mu) * ((alpha - np.sin(alpha)) - (beta - np.sin(beta)))


def _tof_hyperbolic(mu, a, s, c, theta):
    """Evaluate the hyperbolic Lagrange time equation.

    Parameters
    ----------
    mu: float
        Gravitational parameter.
    a: float
        Negative semi-major axis of the hyperbolic transfer orbit.
    s: float
        Semiperimeter of the transfer triangle.
    c: float
        Chord length between the terminal position vectors.
    theta: float
        Transfer angle in radians.

    Returns
    -------
    tof: float
        Hyperbolic transfer time associated with ``a``.

    Notes
    -----
    This is equation (4) of Jiang et al., with alpha and beta computed from
    equation (11).

    """
    alpha, beta = _hyperbolic_alpha_beta(a, s, c, theta)
    return np.sqrt((-a) ** 3 / mu) * ((np.sinh(alpha) - alpha) - (np.sinh(beta) - beta))


def _elliptic_alpha_beta(a, s, c, theta, long_period):
    """Compute elliptic alpha and beta from Jiang's equation (6).

    Parameters
    ----------
    a: float
        Positive semi-major axis of the elliptic transfer orbit.
    s: float
        Semiperimeter of the transfer triangle.
    c: float
        Chord length between the terminal position vectors.
    theta: float
        Transfer angle in radians.
    long_period: bool
        If `True`, use the ``2 * pi - alpha`` branch needed after the minimum
        energy time.

    Returns
    -------
    alpha: float
        Elliptic alpha angle in radians.
    beta: float
        Signed elliptic beta angle in radians.

    """
    alpha = 2 * np.arcsin(np.sqrt(s / (2 * a)))
    if long_period:
        alpha = 2 * np.pi - alpha

    beta = 2 * np.arcsin(np.sqrt((s - c) / (2 * a)))
    beta *= _beta_sign(theta)
    return alpha, beta


def _hyperbolic_alpha_beta(a, s, c, theta):
    """Compute hyperbolic alpha and beta from Jiang's equation (11).

    Parameters
    ----------
    a: float
        Negative semi-major axis of the hyperbolic transfer orbit.
    s: float
        Semiperimeter of the transfer triangle.
    c: float
        Chord length between the terminal position vectors.
    theta: float
        Transfer angle in radians.

    Returns
    -------
    alpha: float
        Hyperbolic alpha parameter.
    beta: float
        Signed hyperbolic beta parameter.

    """
    alpha = 2 * np.arcsinh(np.sqrt(s / (-2 * a)))
    beta = 2 * np.arcsinh(np.sqrt((s - c) / (-2 * a)))
    beta *= _beta_sign(theta)
    return alpha, beta


def _beta_sign(theta):
    """Return the beta sign selected by the transfer angle.

    Parameters
    ----------
    theta: float
        Transfer angle in radians.

    Returns
    -------
    beta_sign: float
        Positive for transfers up to pi radians and negative for transfers
        greater than pi radians.

    Notes
    -----
    Jiang et al. define beta as nonnegative for short-way transfers and
    nonpositive for transfers whose angle exceeds pi.

    """
    return 1.0 if theta <= np.pi else -1.0


def _get_elliptic_bracket(mu, tof, a_min, s, c, theta, long_period):
    """Find a positive semi-major-axis interval containing the elliptic root.

    Parameters
    ----------
    mu: float
        Gravitational parameter.
    tof: float
        Desired transfer time.
    a_min: float
        Minimum-energy elliptic semi-major axis.
    s: float
        Semiperimeter of the transfer triangle.
    c: float
        Chord length between the terminal position vectors.
    theta: float
        Transfer angle in radians.
    long_period: bool
        If `True`, bracket the long-period elliptic time equation.

    Returns
    -------
    lower: float
        Lower semi-major-axis bound.
    upper: float
        Upper semi-major-axis bound.

    Notes
    -----
    Jiang et al. describe expanding a fixed upper bound before applying a
    scalar iteration. This implementation keeps that bracketing idea but uses a
    dynamic upper bound so the solver is scale independent.

    """
    eps = np.sqrt(np.finfo(float).eps)
    lower = a_min * (1 + eps)
    upper = max(2 * lower, 1.0)

    def residual(axis):
        """Evaluate the elliptic time-of-flight residual at ``axis``."""
        return _tof_elliptic(mu, axis, s, c, theta, long_period) - tof

    while residual(lower) * residual(upper) > 0:
        upper *= 2

    return lower, upper


def _get_hyperbolic_bracket(mu, tof, s, c, theta):
    """Find a negative semi-major-axis interval containing the hyperbolic root.

    Parameters
    ----------
    mu: float
        Gravitational parameter.
    tof: float
        Desired transfer time.
    s: float
        Semiperimeter of the transfer triangle.
    c: float
        Chord length between the terminal position vectors.
    theta: float
        Transfer angle in radians.

    Returns
    -------
    lower: float
        More negative semi-major-axis bound.
    upper: float
        Negative bound close to the parabolic limit.

    Notes
    -----
    The paper gives fixed trial values for the hyperbolic interval. The bounds
    here follow the same sign convention while scaling with the problem
    geometry.

    """
    eps = np.sqrt(np.finfo(float).eps)
    lower = -max(s, c, 1.0)
    upper = -eps * max(s, c, 1.0)

    def residual(axis):
        """Evaluate the hyperbolic time-of-flight residual at ``axis``."""
        return _tof_hyperbolic(mu, axis, s, c, theta) - tof

    while residual(lower) < 0:
        lower *= 2

    return lower, upper


def _solve_semimajor_axis(function, lower, upper, maxiter, atol, rtol):
    """Solve the scalar semi-major-axis equation with safeguarded Newton steps.

    Parameters
    ----------
    function: callable
        Residual of the Lagrange time equation, expressed as ``tof(a) - tof``.
    lower: float
        Lower bracket bound for the semi-major axis.
    upper: float
        Upper bracket bound for the semi-major axis.
    maxiter: int
        Maximum number of iterations.
    atol: float
        Absolute tolerance used for residual and step convergence.
    rtol: float
        Relative tolerance used for step convergence.

    Returns
    -------
    a: float
        Semi-major axis that satisfies the selected time equation.
    numiter: int
        Number of iterations used.

    Notes
    -----
    Jiang et al. replace bisection with Newton iteration. The implementation is
    safeguarded by the current bracket, so a Newton step that leaves the
    feasible interval falls back to bisection.

    """
    x = (lower + upper) / 2
    f_lower = function(lower)

    for numiter in range(1, maxiter + 1):
        f_x = function(x)
        if np.abs(f_x) <= atol:
            return x, numiter

        if f_lower * f_x <= 0:
            upper = x
        else:
            lower = x
            f_lower = f_x

        derivative = _numerical_derivative(function, x, lower, upper)
        x_new = x - f_x / derivative if derivative != 0 else np.nan

        if not np.isfinite(x_new) or not lower < x_new < upper:
            x_new = (lower + upper) / 2

        f_new = function(x_new)
        if np.abs(f_new) <= atol or (
            np.abs(x_new - x) <= atol + rtol * np.abs(x_new)
            and np.abs(f_new) <= np.abs(f_x)
        ):
            return x_new, numiter

        x = x_new

    raise ValueError("Exceeded maximum number of iterations!")


def _numerical_derivative(function, x, lower, upper):
    """Differentiate the time equation without stepping outside the bracket.

    Parameters
    ----------
    function: callable
        Scalar residual function to differentiate.
    x: float
        Current semi-major-axis iterate.
    lower: float
        Current lower bracket bound.
    upper: float
        Current upper bracket bound.

    Returns
    -------
    derivative: float
        Centered finite-difference derivative clipped to the active bracket.

    """
    step = np.sqrt(np.finfo(float).eps) * max(1.0, np.abs(x))
    x_minus = max(lower, x - step)
    x_plus = min(upper, x + step)

    if x_plus == x_minus:
        return 0.0

    return (function(x_plus) - function(x_minus)) / (x_plus - x_minus)


def _elliptic_p(a, r1_norm, r2_norm, c, s, theta, long_period):
    """Compute the elliptic semi-latus rectum from Jiang's equation (14).

    Parameters
    ----------
    a: float
        Positive semi-major axis of the elliptic transfer orbit.
    r1_norm: float
        Norm of the initial position vector.
    r2_norm: float
        Norm of the final position vector.
    c: float
        Chord length between the terminal position vectors.
    s: float
        Semiperimeter of the transfer triangle.
    theta: float
        Transfer angle in radians.
    long_period: bool
        If `True`, use the long-period elliptic alpha branch.

    Returns
    -------
    p: float
        Elliptic semi-latus rectum.

    """
    alpha, beta = _elliptic_alpha_beta(a, s, c, theta, long_period)
    factor = 4 * a * (s - r1_norm) * (s - r2_norm) / c**2
    return factor * np.sin((alpha + beta) / 2) ** 2


def _hyperbolic_p(a, r1_norm, r2_norm, c, s, theta):
    """Compute the hyperbolic semi-latus rectum from Jiang's equation (15).

    Parameters
    ----------
    a: float
        Negative semi-major axis of the hyperbolic transfer orbit.
    r1_norm: float
        Norm of the initial position vector.
    r2_norm: float
        Norm of the final position vector.
    c: float
        Chord length between the terminal position vectors.
    s: float
        Semiperimeter of the transfer triangle.
    theta: float
        Transfer angle in radians.

    Returns
    -------
    p: float
        Hyperbolic semi-latus rectum.

    """
    alpha, beta = _hyperbolic_alpha_beta(a, s, c, theta)
    factor = -4 * a * (s - r1_norm) * (s - r2_norm) / c**2
    return factor * np.sinh((alpha + beta) / 2) ** 2


def _parabolic_p(r1_norm, r2_norm, c, s, theta):
    """Compute the parabolic semi-latus rectum as the elliptic limiting value.

    Parameters
    ----------
    r1_norm: float
        Norm of the initial position vector.
    r2_norm: float
        Norm of the final position vector.
    c: float
        Chord length between the terminal position vectors.
    s: float
        Semiperimeter of the transfer triangle.
    theta: float
        Transfer angle in radians.

    Returns
    -------
    p: float
        Parabolic semi-latus rectum.

    Notes
    -----
    Equation (16) in the PDF is typeset as the root of a quadratic in ``p``.
    This expression is the equivalent parabolic limit of equation (14), using
    Jiang's beta sign convention.

    """
    factor = 2 * (s - r1_norm) * (s - r2_norm) / c**2
    return factor * (np.sqrt(s) + _beta_sign(theta) * np.sqrt(s - c)) ** 2


def _reconstruct_velocities(mu, r1, r2, r1_norm, r2_norm, theta, p):
    """Reconstruct endpoint velocities from Jiang's equations (17) and (18).

    Parameters
    ----------
    mu: float
        Gravitational parameter.
    r1: numpy.array
        Initial position vector.
    r2: numpy.array
        Final position vector.
    r1_norm: float
        Norm of the initial position vector.
    r2_norm: float
        Norm of the final position vector.
    theta: float
        Transfer angle in radians.
    p: float
        Semi-latus rectum of the solved conic.

    Returns
    -------
    v1: numpy.array
        Initial velocity vector.
    v2: numpy.array
        Final velocity vector.

    """
    scale = np.sqrt(mu * p) / (r1_norm * r2_norm * np.sin(theta))
    one_minus_cos = 1 - np.cos(theta)

    v1 = scale * ((r2 - r1) + (r2_norm / p) * one_minus_cos * r1)
    v2 = scale * ((r2 - r1) - (r1_norm / p) * one_minus_cos * r2)
    return v1, v2
