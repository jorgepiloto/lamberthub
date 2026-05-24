"""Lambert's problem solver using the method proposed by Gim J. Der in 2011."""

import time

import numpy as np
from numpy.linalg import norm
from scipy.optimize import brentq

from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import (
    assert_parameters_are_valid,
    assert_transfer_angle_not_pi,
    assert_transfer_angle_not_zero,
)


def der2011(
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
    Solve Lambert's problem using Der's superior Lambert algorithm.

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
    This implementation follows Der's ``lambert2`` computational procedure. It
    solves Sun's normalized Lambert equation in the path parameter ``x`` with a
    Laguerre iteration and reconstructs the terminal velocities from radial and
    chord-wise components.

    References
    ----------
    Der, G. J. (2011). The superior Lambert algorithm. AMOS, Maui, Hawaii.

    """
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    r1_norm, r2_norm, c_norm = [norm(vec) for vec in (r1, r2, r2 - r1)]
    psi = get_transfer_angle(r1, r2, is_prograde)
    assert_transfer_angle_not_zero(psi)
    assert_transfer_angle_not_pi(psi)

    m = r1_norm + r2_norm + c_norm
    n = r1_norm + r2_norm - c_norm
    sigma = _compute_sigma(r1_norm, r2_norm, m, psi)
    tau = _compute_normalized_time(mu, m, tof)
    tau_p = _compute_parabolic_time(sigma)

    if M > 0 and tau <= tau_p:
        raise ValueError("No feasible solution, try lower M!")

    if M > 0:
        tau_min = _compute_minimum_time(sigma, M, maxiter, atol, rtol)
        if tau < tau_min:
            raise ValueError("No feasible solution, try lower M!")

    tic = time.perf_counter()
    x, y, numiter = _find_path_parameter(
        sigma,
        tau,
        tau_p,
        M,
        is_low_path,
        maxiter,
        atol,
        rtol,
    )
    tac = time.perf_counter()
    tpi = (tac - tic) / numiter

    v1, v2 = _reconstruct_velocities(mu, r1, r2, r1_norm, r2_norm, c_norm, m, n, x, y)

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def _compute_sigma(r1_norm, r2_norm, m, psi):
    """Compute Der's signed angle parameter."""
    return 2 * np.sqrt(r1_norm * r2_norm) * np.cos(psi / 2) / m


def _compute_normalized_time(mu, m, tof):
    """Compute Sun's normalized time of flight."""
    return 4 * tof * np.sqrt(mu / m**3)


def _compute_parabolic_time(sigma):
    """Compute the normalized parabolic time of flight."""
    return 2 * (1 - sigma**3) / 3


def _find_path_parameter(
    sigma,
    tau,
    tau_p,
    M,
    is_low_path,
    maxiter,
    atol,
    rtol,
):
    """Solve for Der's path parameter with the Laguerre method."""
    if M == 0 and np.isclose(tau, tau_p, atol=atol, rtol=rtol):
        return 1.0, sigma, 1

    is_elliptic = M > 0 or tau > tau_p
    tau_me = _compute_minimum_energy_time(sigma, M) if is_elliptic else None

    if M == 0:
        x0 = -0.5 if is_elliptic and tau > tau_me else 0.5
    else:
        x0 = 0.5 if is_low_path else -0.5

    if is_elliptic:
        return _solve_elliptic(sigma, tau, M, x0, maxiter, atol, rtol)

    return _solve_hyperbolic(sigma, tau, x0, maxiter, atol, rtol)


def _solve_elliptic(sigma, tau, M, x0, maxiter, atol, rtol):
    """Solve Sun's elliptic Lambert equation."""
    x = x0
    for numiter in range(1, maxiter + 1):
        y = _compute_y(sigma, x)
        tof = _tof_elliptic(sigma, x, y, M)
        f = tof - tau

        if np.abs(f) <= atol + rtol * np.abs(tau):
            return x, y, numiter

        dt = _tof_elliptic_prime(sigma, x, y, tof)
        d2t = _tof_elliptic_second(sigma, x, y, dt)
        try:
            x = _laguerre_step(x, f, dt, d2t)
        except ValueError:
            return _solve_elliptic_bracketed(sigma, tau, M, x0, maxiter, atol, rtol)

        if np.abs(x) >= 1:
            return _solve_elliptic_bracketed(sigma, tau, M, x0, maxiter, atol, rtol)

    raise ValueError("Exceeded maximum number of iterations.")


def _solve_elliptic_bracketed(sigma, tau, M, x0, maxiter, atol, rtol):
    """Solve Sun's elliptic Lambert equation with Der's path-region bounds."""
    eps = np.finfo(float).eps
    lower, upper = (-1 + eps, 0.0) if x0 < 0 else (0.0, 1 - eps)

    def elliptic_equation(x):
        """Evaluate Sun's elliptic Lambert equation."""
        y = _compute_y(sigma, x)
        return _tof_elliptic(sigma, x, y, M) - tau

    try:
        x = brentq(
            elliptic_equation,
            lower,
            upper,
            xtol=atol,
            rtol=rtol,
            maxiter=maxiter,
        )
    except ValueError as exc:
        raise ValueError("Exceeded maximum number of iterations.") from exc

    y = _compute_y(sigma, x)
    return x, y, maxiter


def _solve_hyperbolic(sigma, tau, x0, maxiter, atol, rtol):
    """Solve Sun's hyperbolic Lambert equation."""
    x = max(1.1, x0 + 1.0)
    for numiter in range(1, maxiter + 1):
        y = _compute_y(sigma, x)
        tof = _tof_hyperbolic(sigma, x, y)
        f = tof - tau

        if np.abs(f) <= atol + rtol * np.abs(tau):
            return x, y, numiter

        dt, d2t = _numerical_derivatives(_tof_hyperbolic, sigma, x, y)
        try:
            x = max(1.0 + np.finfo(float).eps, _laguerre_step(x, f, dt, d2t))
        except ValueError:
            return _solve_hyperbolic_bracketed(sigma, tau, maxiter, atol, rtol)

    raise ValueError("Exceeded maximum number of iterations.")


def _solve_hyperbolic_bracketed(sigma, tau, maxiter, atol, rtol):
    """Solve Sun's hyperbolic Lambert equation on the valid path interval."""
    lower = 1.0 + np.sqrt(np.finfo(float).eps)
    upper = 2.0

    def hyperbolic_equation(x):
        """Evaluate Sun's hyperbolic Lambert equation."""
        y = _compute_y(sigma, x)
        return _tof_hyperbolic(sigma, x, y) - tau

    while hyperbolic_equation(upper) > 0:
        upper *= 2

    try:
        x = brentq(
            hyperbolic_equation,
            lower,
            upper,
            xtol=atol,
            rtol=rtol,
            maxiter=maxiter,
        )
    except ValueError as exc:
        raise ValueError("Exceeded maximum number of iterations.") from exc

    y = _compute_y(sigma, x)
    return x, y, maxiter


def _laguerre_step(x, f, df, d2f):
    """Compute one Laguerre update for the current scalar equation value."""
    for degree in range(2, 13):
        discriminant = (degree - 1) ** 2 * df**2 - degree * (degree - 1) * f * d2f
        if discriminant < 0:
            continue

        root = np.sqrt(discriminant)
        denominator = df + np.copysign(root, df)
        if denominator == 0:
            continue

        candidate = x - degree * f / denominator
        if np.isfinite(candidate):
            return candidate

    raise ValueError("Laguerre iteration failed to compute a finite update.")


def _compute_y(sigma, x):
    """Compute Sun's dependent variable from the path and angle parameters."""
    radicand = 1 - sigma**2 * (1 - x**2)
    y_abs = np.sqrt(max(0.0, radicand))
    if sigma == 0:
        return y_abs
    return np.copysign(y_abs, sigma)


def _tof_elliptic(sigma, x, y, M):
    """Evaluate Sun's normalized elliptic time of flight."""
    one_minus_x2 = 1 - x**2
    one_minus_y2 = 1 - y**2
    numerator = (
        _acot_elliptic_x(x)
        - _acot_elliptic_y(y)
        - x * np.sqrt(one_minus_x2)
        + y * np.sqrt(max(0.0, one_minus_y2))
        + M * np.pi
    )
    return numerator / one_minus_x2 ** (3 / 2)


def _tof_hyperbolic(sigma, x, y):
    """Evaluate Sun's normalized hyperbolic time of flight."""
    x2_minus_one = x**2 - 1
    y2_minus_one = y**2 - 1
    numerator = (
        -_acoth_sun(x)
        + _acoth_sun(y)
        + x * np.sqrt(x2_minus_one)
        - y * np.sqrt(max(0.0, y2_minus_one))
    )
    return numerator / x2_minus_one ** (3 / 2)


def _tof_elliptic_prime(sigma, x, y, tof):
    """Evaluate the first derivative of Sun's elliptic time equation."""
    return (3 * x * tof - 2 * (1 - sigma**3 * x / y)) / (1 - x**2)


def _tof_elliptic_second(sigma, x, y, dtof):
    """Evaluate the second derivative of Sun's elliptic time equation."""
    if x == 0:
        x = np.finfo(float).eps
    return ((1 + 4 * x**2) * dtof + 2 * (1 - sigma**5 * x**3 / y**3)) / (x * (1 - x**2))


def _numerical_derivatives(tof_function, sigma, x, y):
    """Approximate first and second derivatives of a time equation."""
    step = np.sqrt(np.finfo(float).eps) * max(1.0, np.abs(x))
    x_minus = max(1.0 + np.finfo(float).eps, x - step)
    x_plus = x + step

    y_minus = _compute_y(sigma, x_minus)
    y_plus = _compute_y(sigma, x_plus)
    f_minus = tof_function(sigma, x_minus, y_minus)
    f = tof_function(sigma, x, y)
    f_plus = tof_function(sigma, x_plus, y_plus)

    first = (f_plus - f_minus) / (x_plus - x_minus)
    second = (
        2
        * ((f_plus - f) / (x_plus - x) - (f - f_minus) / (x - x_minus))
        / (x_plus - x_minus)
    )
    return first, second


def _acot_elliptic_x(x):
    """Evaluate the elliptic arc-cotangent branch for the path parameter."""
    return np.arccos(np.clip(x, -1.0, 1.0))


def _acot_elliptic_y(y):
    """Evaluate the elliptic arc-cotangent branch for the angle parameter."""
    if y >= 0:
        return np.arccos(np.clip(y, -1.0, 1.0))
    return -np.arccos(np.clip(-y, -1.0, 1.0))


def _acoth_sun(u):
    """Evaluate the inverse hyperbolic cotangent branch used by Sun."""
    if u >= 0:
        return np.arccosh(u)
    return -np.arccosh(-u)


def _phi(u):
    """Evaluate Der's minimum-time auxiliary function."""
    sqrt_term = np.sqrt(max(0.0, 1 - u**2))
    if u == 0:
        u = np.finfo(float).eps
    return _acot_elliptic_y(u) - (2 + u**2) * sqrt_term / (3 * u)


def _compute_minimum_time(sigma, M, maxiter, atol, rtol):
    """Compute the normalized minimum time for a multirevolution transfer."""

    def minimum_time_equation(x):
        """Evaluate Der's minimum-time equation."""
        y = _compute_y(sigma, x)
        return _phi(y) - _phi(x) - M * np.pi

    try:
        x = brentq(
            minimum_time_equation,
            np.finfo(float).eps,
            1 - np.finfo(float).eps,
            xtol=atol,
            rtol=rtol,
            maxiter=maxiter,
        )
    except ValueError:
        x = 0.1
    else:
        y = _compute_y(sigma, x)
        return 2 * (1 / x - sigma**3 / y) / 3

    for _ in range(maxiter):
        y = _compute_y(sigma, x)
        phi = minimum_time_equation(x)
        dphi = 2 * (1 - x**2) * (1 - sigma**5 * x**3 / y**3) / (3 * x**2)

        if np.abs(phi) <= atol + rtol * M * np.pi:
            return 2 * (1 / x - sigma**3 / y) / 3

        x -= phi / dphi
        x = np.clip(x, np.finfo(float).eps, 0.999999)

    raise ValueError("Exceeded maximum number of iterations.")


def _compute_minimum_energy_time(sigma, M):
    """Compute the normalized minimum-energy time."""
    return (
        M * np.pi
        + np.arccos(np.clip(sigma, -1.0, 1.0))
        + sigma * np.sqrt(max(0.0, 1 - sigma**2))
    )


def _reconstruct_velocities(mu, r1, r2, r1_norm, r2_norm, c_norm, m, n, x, y):
    """Compute the terminal velocity vectors from radial and chord components."""
    vc = np.sqrt(mu) * (y / np.sqrt(n) + x / np.sqrt(m))
    vr = np.sqrt(mu) * (y / np.sqrt(n) - x / np.sqrt(m))

    i_c = (r2 - r1) / c_norm
    i_r1 = r1 / r1_norm
    i_r2 = r2 / r2_norm

    v1 = vc * i_c + vr * i_r1
    v2 = vc * i_c - vr * i_r2
    return v1, v2
