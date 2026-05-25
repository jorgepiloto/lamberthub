"""A module hosting the universal solver by McElreath, Down, and Majji."""

import time

from numba import njit as jit
import numpy as np

from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import (
    assert_parameters_are_valid,
    assert_transfer_angle_not_zero,
)


def mcelreath2025(
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
    Solve Lambert's problem using McElreath, Down, and Majji's method.

    Parameters
    ----------
    mu: float
        Gravitational parameter, equivalent to :math:`GM` of attractor body.
    r1: numpy.array
        Initial position vector.
    r2: numpy.array
        Final position vector.
    tof: float
        Time of flight between initial and final position vectors.
    M: int
        Number of revolutions. Must be equal or greater than 0 value.
    is_prograde: bool
        If `True`, specifies prograde motion. Otherwise, retrograde motion is imposed.
    is_low_path: bool
        If two solutions are available, it selects between high or low path.
    maxiter: int
        Maximum number of iterations for the scalar root solvers.
    atol: float
        Absolute tolerance used by the scalar root solvers.
    rtol: float
        Relative tolerance used by the scalar root solvers.
    full_output: bool
        If True, the number of iterations and time per iteration are also returned.

    Returns
    -------
    v1: numpy.array
        Initial velocity vector.
    v2: numpy.array
        Final velocity vector.
    numiter: int
        Number of scalar root-solver iterations.
    tpi: float
        Time per iteration in seconds.

    Notes
    -----
    This implementation uses the paper's universal transfer-time equation,
    :math:`T(x) = 2 U_1^* m + U_3`, together with the Stark velocity
    reconstruction. The root is found with bracketing methods over the bounded
    universal-variable intervals described in the paper.

    References
    ----------
    [1] McElreath, J., Down, I. M., & Majji, M. (2025). A universal approach
        for solving the multi-revolution Lambert's problem. Celestial Mechanics
        and Dynamical Astronomy, 137, 22.
    """
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)
    dtheta = get_transfer_angle(r1, r2, is_prograde)
    assert_transfer_angle_not_zero(dtheta)

    if np.isclose(np.sin(dtheta), 0.0, atol=1e-14):
        raise RuntimeError("Cannot compute orbit, phase angle is 180 degrees")

    m = (-1) ** M * np.sqrt(r1_norm * r2_norm) * np.cos(dtheta / 2)
    tof = np.sqrt(mu) * tof

    tic = time.perf_counter() if full_output else 0.0
    x, numiter = _find_x(
        tof, r1_norm, r2_norm, m, dtheta, M, is_low_path, maxiter, atol, rtol
    )
    tpi = (time.perf_counter() - tic) / numiter if full_output else 0.0

    _, _, u0_star, _ = _alpha_and_universal_values(
        np.complex128(x), r1_norm, r2_norm, m
    )
    p = (
        r1_norm
        * r2_norm
        * (1 - np.cos(dtheta))
        / (r1_norm + r2_norm - 2 * m * u0_star.real)
    )

    c_norm = np.linalg.norm(r2 - r1)
    v_c = c_norm * np.sqrt(mu * p) / (r1_norm * r2_norm * np.sin(dtheta))
    v_r = np.sqrt(mu / p) * (1 - np.cos(dtheta)) / np.sin(dtheta)

    i_r1 = r1 / r1_norm
    i_r2 = r2 / r2_norm
    i_c = (r2 - r1) / c_norm

    v1 = v_c * i_c + v_r * i_r1
    v2 = v_c * i_c - v_r * i_r2

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


@jit(cache=True, fastmath=True)
def _find_x(tof, r1_norm, r2_norm, m, dtheta, M, is_low_path, maxiter, atol, rtol):
    """Find McElreath's universal variable for the requested transfer."""
    r_norm_sum = r1_norm + r2_norm
    if M == 0:
        return _find_single_revolution_x(
            tof, r_norm_sum, m, dtheta, maxiter, atol, rtol
        )

    return _find_multi_revolution_x(
        tof, r_norm_sum, m, dtheta, M, is_low_path, maxiter, atol, rtol
    )


@jit(cache=True, fastmath=True)
def _find_single_revolution_x(tof, r_norm_sum, m, dtheta, maxiter, atol, rtol):
    """Find the single-revolution elliptic or hyperbolic root."""
    tof_p = _get_parabolic_tof(r_norm_sum, m)
    if abs(tof - tof_p) <= atol + rtol * abs(tof_p):
        return np.complex128(0.0), 1

    if tof > tof_p:
        tof_theta = _tof_equation(np.complex128(dtheta), r_norm_sum, m)
        beta = np.sqrt((tof - tof_p) / (tof_theta.real - tof_p))
        beta *= dtheta / (2 * np.pi - dtheta)
        x_0 = np.complex128(2 * np.pi * beta / (1 + beta))
        lo, hi = 0.0, 2 * np.pi
    elif dtheta < np.pi:
        tau = r_norm_sum**2 / (2 * m**2) - 1
        x_tau = 1j * np.log(tau + np.sqrt(tau**2 - 1))
        x_1 = x_tau * np.exp(-2.5)
        tof_x1 = _tof_equation(np.complex128(x_1), r_norm_sum, m)
        tof_p_16 = np.complex128(tof_p) ** 1.6
        tof_x1_16 = np.complex128(tof_x1) ** 1.6
        exponent = 2.5 / (np.log(-tof_p_16) - np.log(tof_x1_16 - tof_p_16))
        x_0 = x_tau * (1 - (tof / tof_p) ** 1.6) ** exponent
        x_0 = np.complex128(1j * x_0.imag)
        lo, hi = 0.0, x_tau.imag
    else:
        tof_theta = _tof_equation(np.complex128(1j * dtheta), r_norm_sum, m)
        x_0 = np.complex128(
            1j * dtheta * np.sqrt(np.log(tof / tof_p) / np.log(tof_theta / tof_p))
        )
        lo, hi = 0.0, np.inf

    return _solve_desired_x(x_0, tof, r_norm_sum, m, 0, lo, hi, maxiter, rtol)


@jit(cache=True, fastmath=True)
def _find_multi_revolution_x(
    tof, r_norm_sum, m, dtheta, M, is_low_path, maxiter, atol, rtol
):
    """Find one branch of the multi-revolution root pair."""
    two_pi = 2 * np.pi
    two_pi_revs = two_pi * M
    sqrt_revs = np.sqrt(M)
    lo, hi = two_pi_revs, two_pi_revs + two_pi

    if dtheta < np.pi:
        x_min = (
            2 * (dtheta + (1 / 6) ** (2.5 * sqrt_revs)) ** 0.4
            - (1 / 3) ** sqrt_revs
            + two_pi_revs
        )
    else:
        x_min = (
            -2 * (-dtheta + two_pi + 0.05**sqrt_revs) ** 0.4
            + (two_pi - 0.25**sqrt_revs)
            + two_pi_revs
        )

    x_min, minimum_iterations = _solve_minimum_x(
        np.complex128(x_min), tof, r_norm_sum, m, M, lo, hi, maxiter
    )
    tof_min = _tof_equation(np.complex128(x_min), r_norm_sum, m).real

    if tof < tof_min - atol:
        raise ValueError("No feasible solution, try lower M!")

    if abs(tof - tof_min) <= atol + rtol * abs(tof_min):
        if minimum_iterations < 1:
            minimum_iterations = 1
        return np.complex128(x_min), minimum_iterations

    g_short, g_long = _multi_revolution_initial_guesses(
        tof, tof_min, x_min, r_norm_sum, m, dtheta, M
    )
    x_0 = np.complex128(g_short if is_low_path else g_long)
    x, desired_iterations = _solve_desired_x(
        x_0, tof, r_norm_sum, m, M, lo, hi, maxiter, rtol
    )
    total = minimum_iterations + desired_iterations
    if total < 1:
        total = 1
    return x, total


@jit(cache=True, fastmath=True)
def _solve_minimum_x(x, tof, r_norm_sum, m, M, lo, hi, maxiter):
    """Solve for the multi-revolution minimum transfer-time location."""
    two_pi_revs = 2 * np.pi * M
    tolerance = 1e-2
    for numiter in range(1, maxiter + 1):
        tof_current, tof_prime, tof_2prime, tof_3prime = _tof_equation_full(
            x, r_norm_sum, m
        )
        step = (
            -2 * tof_prime * tof_2prime / (2 * tof_2prime**2 - tof_prime * tof_3prime)
        )
        x_old = x
        x = x + step

        if _outside_bounds(x, lo, hi):
            x = x_old - tof_prime / tof_2prime

        err = abs((x - x_old) / (x - two_pi_revs + 2 * np.pi))
        relative_gap = (tof - tof_current.real) / tof
        if relative_gap < 1e-4:
            tolerance = 1e-8
        if err < tolerance:
            return x.real, numiter

    raise RuntimeError("Failed to converge")


@jit(cache=True, fastmath=True)
def _multi_revolution_initial_guesses(tof, tof_min, x_min, r_norm_sum, m, dtheta, M):
    """Compute the original short- and long-period initial guesses."""
    two_pi = 2 * np.pi
    two_pi_revs = two_pi * M
    x_theta = dtheta + two_pi_revs
    tof_theta = _tof_equation(np.complex128(x_theta), r_norm_sum, m).real
    scale = (x_min - x_theta) ** 2 / (
        abs(tof_theta - tof_min)
        * (two_pi_revs - x_theta) ** 2
        * (two_pi_revs + two_pi - x_theta) ** 2
    )
    gamma = -np.sqrt((tof - tof_min) * scale)
    next_revs = M + 1

    b_short = -two_pi * (M + next_revs) * gamma + 1
    c_short = two_pi**2 * M * next_revs * gamma - x_min
    x_short = (-b_short + np.sqrt(b_short**2 - 4 * gamma * c_short)) / (2 * gamma)

    b_long = two_pi * (M + next_revs) * gamma + 1
    c_long = -(two_pi**2) * M * next_revs * gamma - x_min
    x_long = (-b_long + np.sqrt(b_long**2 + 4 * gamma * c_long)) / (-2 * gamma)

    return x_short, x_long


@jit(cache=True, fastmath=True)
def _solve_desired_x(x, tof, r_norm_sum, m, M, lo, hi, maxiter, rtol):
    """Solve for the requested transfer time with Householder updates."""
    two_pi_revs = 2 * np.pi * M
    tolerance = rtol if rtol < 1e-10 else 1e-10

    for numiter in range(1, maxiter + 1):
        tof_current, tof_prime, tof_2prime, tof_3prime = _tof_equation_full(
            x, r_norm_sum, m
        )
        residual = tof_current - tof
        step = -(
            (6 * residual * tof_prime**2 - 3 * residual**2 * tof_2prime)
            / (
                6 * tof_prime**3
                - 6 * residual * tof_prime * tof_2prime
                + residual**2 * tof_3prime
            )
        )
        x_old = x
        x = x + step

        if _outside_bounds(x, lo, hi):
            x = x_old - 2 * residual * tof_prime / (
                2 * tof_prime**2 - residual * tof_2prime
            )
            if _outside_bounds(x, lo, hi):
                x = x_old - residual / tof_prime
                if _outside_bounds(x, lo, hi):
                    rt = -(residual / tof_prime)
                    sign = np.sign(rt.real + rt.imag)
                    boundary = hi if sign > 0 else lo
                    x = x_old + 0.5 * (boundary - x_old)

        err = abs((x - x_old) / (x - two_pi_revs + 2 * np.pi))
        if err < tolerance:
            return _clean_root(x), numiter

    raise RuntimeError("Failed to converge")


@jit(cache=True, fastmath=True)
def _outside_bounds(x, lo, hi):
    """Check if an update stepped outside the solver bounds."""
    coordinate = x.real + x.imag
    return coordinate > hi or coordinate < lo


@jit(cache=True, fastmath=True)
def _clean_root(x):
    """Return a real root (as complex) when the imaginary component is noise."""
    if abs(x.imag) < 1e-14:
        return np.complex128(x.real)
    return x


@jit(cache=True, fastmath=True)
def _get_parabolic_tof(r_norm_sum, m):
    """Compute the scaled parabolic time of flight."""
    chi = np.sqrt(2 * (r_norm_sum - 2 * m))
    return m * chi + chi**3 / 6


@jit(cache=True, fastmath=True)
def _tof_at_x(x, r1_norm, r2_norm, m):
    """Evaluate McElreath's scaled transfer-time equation at x."""
    if abs(x) < 1e-8:
        return _get_parabolic_tof(r1_norm + r2_norm, m)

    _, u1_star, _, u3 = _alpha_and_universal_values(
        np.complex128(x), r1_norm, r2_norm, m
    )
    return (2 * u1_star * m + u3).real


@jit(cache=True, fastmath=True)
def _tof_equation(x, r_norm_sum, m):
    """Evaluate the reference scaled transfer-time function (complex result)."""
    u0 = np.cos(x)
    u0_star = np.cos(0.5 * x)
    alpha = (1 - u0) / (r_norm_sum - 2 * m * u0_star)
    alpha_inv = 1 / alpha
    alpha_sqrt = np.sqrt(alpha)
    alpha_sqrt_inv = 1 / alpha_sqrt

    u1_star = np.sin(0.5 * x) * alpha_sqrt_inv
    u3 = (x * alpha_sqrt_inv - np.sin(x) * alpha_sqrt_inv) * alpha_inv
    return 2 * u1_star * m + u3


@jit(cache=True, fastmath=True)
def _tof_equation_full(x, r_norm_sum, m):
    """Evaluate the reference scaled transfer-time function and derivatives."""
    u0 = np.cos(x)
    u0_star = np.cos(0.5 * x)
    alpha = (1 - u0) / (r_norm_sum - 2 * m * u0_star)
    alpha_inv = 1 / alpha
    alpha_sqrt = np.sqrt(alpha)
    alpha_sqrt_inv = 1 / alpha_sqrt

    u1 = np.sin(x) * alpha_sqrt_inv
    u1_star = np.sin(0.5 * x) * alpha_sqrt_inv
    u2 = (1 - u0) * alpha_inv
    u3 = (x * alpha_sqrt_inv - u1) * alpha_inv
    tof_current = 2 * u1_star * m + u3

    u1_star_inv = 1 / u1_star
    alpha_prime = 0.5 * u1_star_inv * alpha_sqrt * (2 * u0_star - m * alpha)
    alpha_2prime = (
        -0.5 * alpha
        + (0.5 * u0_star * alpha_sqrt_inv - m * alpha_sqrt) * alpha_prime * u1_star_inv
    )
    alpha_3prime = (
        0.25 * m * u1_star_inv**2 * (m * alpha - 0.5 * u0_star)
        - 0.5
        * alpha_sqrt_inv
        * u1_star_inv
        * (0.5 * u0_star * alpha_inv + m)
        * alpha_prime
        - 0.75
    ) * alpha_prime + (0.5 * alpha_sqrt_inv * u1_star_inv) * (
        u0_star - 2 * m * alpha
    ) * alpha_2prime

    tof_prime = alpha_inv * (
        (m * u0_star + u2) * alpha_sqrt - 0.5 * (tof_current + 2 * u3) * alpha_prime
    )
    tof_2prime = alpha_inv * (
        (-3 * tof_prime + 2 * m * u0_star * alpha_sqrt_inv) * alpha_prime
        - 0.75 * tof_current * alpha_inv * alpha_prime**2
        + u1
        - 0.5 * m * alpha * u1_star
        - (0.5 * tof_current + u3) * alpha_2prime
    )
    tof_3prime = alpha_inv * (
        u0 * alpha_sqrt_inv
        - 0.25 * m * u0_star * alpha_sqrt
        - alpha_prime
        * (
            4.5 * tof_2prime
            + 1.5 * m * u1_star
            + (
                2.25 * tof_prime * alpha_inv
                - 0.375 * alpha_inv**2 * tof_current * alpha_prime
            )
            * alpha_prime
        )
        + (
            2 * m * u0_star * alpha_sqrt_inv
            - u2 * alpha_sqrt_inv
            - 3.5 * tof_prime
            + (u3 * alpha_inv - 1.75 * tof_current * alpha_inv) * alpha_prime
        )
        * alpha_2prime
        - (0.5 * tof_current + u3) * alpha_3prime
    )
    return tof_current, tof_prime, tof_2prime, tof_3prime


@jit(cache=True, fastmath=True)
def _alpha_and_universal_values(x, r1_norm, r2_norm, m):
    """Evaluate alpha and the universal functions needed by the solver."""
    u0 = np.cos(x)
    u0_star = np.cos(x / 2)
    alpha = (1 - u0) / (r1_norm + r2_norm - 2 * m * u0_star)
    alpha_sqrt = np.sqrt(alpha)
    u1_star = np.sin(x / 2) / alpha_sqrt
    u3 = (x - np.sin(x)) / (alpha * alpha_sqrt)
    return alpha, u1_star, u0_star, u3
