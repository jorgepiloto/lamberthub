"""A module hosting the regularized Lambert solver by De La Torre et al."""

from math import factorial
import time

import numpy as np

from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import assert_parameters_are_valid


def delatorre2018(
    mu,
    r1,
    r2,
    tof,
    M=0,
    is_prograde=True,
    is_low_path=True,
    maxiter=50,
    atol=1e-12,
    rtol=1e-12,
    full_output=False,
):
    r"""
    Solve Lambert's problem using De La Torre's regularized formulation.

    Parameters
    ----------
    mu: float
        Gravitational parameter, equivalent to :math:`GM` of attractor body.
    r1: numpy.array
        Initial position vector.
    r2: numpy.array
        Final position vector.
    tof: float
        Time of flight.
    M: int
        Number of revolutions. Must be equal or greater than 0.
    is_prograde: bool
        If ``True``, specifies prograde motion. Otherwise, retrograde motion is
        imposed.
    is_low_path: bool
        If two multi-revolution solutions are available, selects the low-path
        branch.
    maxiter: int
        Maximum number of iterations.
    atol: float
        Absolute tolerance.
    rtol: float
        Relative tolerance.
    full_output: bool
        If True, the number of iterations and time per iteration are returned.

    Returns
    -------
    v1: numpy.array
        Initial velocity vector.
    v2: numpy.array
        Final velocity vector.

    Notes
    -----
    The free parameter ``z`` is the square of half the variation in eccentric
    anomaly on the conic solving the transfer.  It is obtained from Eq. (68) of
    [1] with a safeguarded Newton iteration on a bracketed interval.

    References
    ----------
    [1] De La Torre, D., Flores, R., & Fantino, E. (2018). On the solution of
        Lambert's problem by regularization. Acta Astronautica, 153, 26-38.

    """
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    tic = time.perf_counter()

    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)
    dtheta = get_transfer_angle(r1, r2, is_prograde)

    P = r1_norm + r2_norm
    Q = 2 * np.sqrt(r1_norm * r2_norm) * np.cos(dtheta / 2)

    z, numiter = _find_z(P, Q, mu, tof, M, is_low_path, maxiter, atol, rtol)

    V = P - Q * _stumpff0(z)
    p = r1_norm * r2_norm * (1 - np.cos(dtheta)) / V

    f = 1 - r2_norm * (1 - np.cos(dtheta)) / p
    g = r1_norm * r2_norm * np.sin(dtheta) / np.sqrt(mu * p)
    g_dot = 1 - r1_norm * (1 - np.cos(dtheta)) / p

    v1 = (r2 - f * r1) / g
    v2 = (g_dot * r2 - r1) / g

    tac = time.perf_counter()
    tpi = (tac - tic) / numiter

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def _find_z(P, Q, mu, tof, M, is_low_path, maxiter, atol, rtol):
    """Find the regularized free parameter ``z``."""
    if M == 0:
        z_low, z_high = _initial_bracket(P, Q, mu, tof)
    else:
        z_min, tof_min = _find_minimum_time(P, Q, mu, M, atol, rtol)
        if tof < tof_min:
            raise ValueError("No feasible solution, try lower M!")

        eps = 1e-12
        if is_low_path:
            z_low, z_high = eps, z_min
        else:
            z_low, z_high = z_min, np.pi**2 - eps

    return _solve_bracketed(P, Q, mu, tof, M, z_low, z_high, maxiter, atol, rtol)


def _initial_bracket(P, Q, mu, tof):
    """Return a monotonic bracket for the single-revolution root."""
    eps = 1e-12
    z_high = np.pi**2 - eps
    z_low = -1.0

    if Q > 0:
        z_floor = -(np.arccosh(P / Q) ** 2)
        z_low = z_floor + max(1e-12, abs(z_floor) * 1e-12)

    while _tof_equation_y(P, Q, mu, z_low, 0) > tof:
        z_low *= 2

    if _tof_equation_y(P, Q, mu, z_high, 0) < tof:
        raise ValueError("No feasible solution, try lower M!")

    return z_low, z_high


def _solve_bracketed(P, Q, mu, tof, M, z_low, z_high, maxiter, atol, rtol):
    """Safeguarded Newton iteration on a bracketed time-of-flight interval."""
    f_low = _tof_equation_y(P, Q, mu, z_low, M) - tof
    f_high = _tof_equation_y(P, Q, mu, z_high, M) - tof

    if f_low == 0:
        return z_low, 1
    if f_high == 0:
        return z_high, 1
    if f_low * f_high > 0:
        raise ValueError("No feasible solution, try lower M!")

    z = 0.5 * (z_low + z_high)
    for numiter in range(1, maxiter + 1):
        f = _tof_equation_y(P, Q, mu, z, M) - tof
        if abs(f) <= atol + rtol * tof:
            return z, numiter

        if f_low * f <= 0:
            z_high = z
            f_high = f
        else:
            z_low = z
            f_low = f

        df = _tof_derivative(P, Q, mu, z, M, z_low, z_high)
        z_new = z - f / df if df != 0 and np.isfinite(df) else np.nan
        if not np.isfinite(z_new) or z_new <= z_low or z_new >= z_high:
            z_new = 0.5 * (z_low + z_high)

        z = z_new

    raise ValueError("Exceeded maximum number of iterations!")


def _find_minimum_time(P, Q, mu, M, atol, rtol):
    """Locate the multi-revolution minimum of the transfer-time curve."""
    z_low = 1e-12
    z_high = np.pi**2 - 1e-12
    gr = (np.sqrt(5) - 1) / 2

    z1 = z_high - gr * (z_high - z_low)
    z2 = z_low + gr * (z_high - z_low)
    t1 = _tof_equation_y(P, Q, mu, z1, M)
    t2 = _tof_equation_y(P, Q, mu, z2, M)

    for _ in range(100):
        if t1 > t2:
            z_low = z1
            z1 = z2
            t1 = t2
            z2 = z_low + gr * (z_high - z_low)
            t2 = _tof_equation_y(P, Q, mu, z2, M)
        else:
            z_high = z2
            z2 = z1
            t2 = t1
            z1 = z_high - gr * (z_high - z_low)
            t1 = _tof_equation_y(P, Q, mu, z1, M)

        if abs(z_high - z_low) <= atol + rtol * abs(z1):
            break

    z_min = 0.5 * (z_low + z_high)
    return z_min, _tof_equation_y(P, Q, mu, z_min, M)


def _tof_derivative(P, Q, mu, z, M, z_low, z_high):
    """Compute a local numerical derivative of the transfer-time equation."""
    dz = max(1e-7, np.sqrt(np.finfo(float).eps) * max(1.0, abs(z)))
    dz = min(dz, 0.25 * (z - z_low), 0.25 * (z_high - z))
    if dz <= 0:
        return 0.0
    return (
        _tof_equation_y(P, Q, mu, z + dz, M) - _tof_equation_y(P, Q, mu, z - dz, M)
    ) / (2 * dz)


def _tof_equation_y(P, Q, mu, z, M):
    """Evaluate Eq. (68), the regularized transfer-time equation."""
    c0 = _stumpff0(z)
    c1 = _stumpff1(z)
    c2_4z = _stumpff2(4 * z)
    c3_4z = _stumpff3(4 * z)
    V = P - Q * c0

    first = (
        (2 * P * c3_4z + Q * (c1 * c2_4z - 2 * c0 * c3_4z))
        * np.sqrt(2 * V / mu)
        / c1**3
    )

    if M == 0:
        return first

    second = M * np.pi * np.sqrt(V**3 / (2 * mu * z**3)) / c1**3
    return first + second


def _stumpff0(z):
    """Evaluate the zeroth Stumpff function."""
    if z > 1e-8:
        return np.cos(np.sqrt(z))
    if z < -1e-8:
        return np.cosh(np.sqrt(-z))
    return _stumpff(0, z)


def _stumpff1(z):
    """Evaluate the first Stumpff function."""
    if z > 1e-8:
        sqrt_z = np.sqrt(z)
        return np.sin(sqrt_z) / sqrt_z
    if z < -1e-8:
        sqrt_neg_z = np.sqrt(-z)
        return np.sinh(sqrt_neg_z) / sqrt_neg_z
    return _stumpff(1, z)


def _stumpff2(z):
    """Evaluate the second Stumpff function."""
    if z > 1e-8:
        return (1 - np.cos(np.sqrt(z))) / z
    if z < -1e-8:
        return (np.cosh(np.sqrt(-z)) - 1) / (-z)
    return _stumpff(2, z)


def _stumpff3(z):
    """Evaluate the third Stumpff function."""
    if z > 1e-8:
        sqrt_z = np.sqrt(z)
        return (sqrt_z - np.sin(sqrt_z)) / (z * sqrt_z)
    if z < -1e-8:
        sqrt_neg_z = np.sqrt(-z)
        return (np.sinh(sqrt_neg_z) - sqrt_neg_z) / (-z * sqrt_neg_z)
    return _stumpff(3, z)


def _stumpff(k, z):
    """Evaluate a Stumpff function by its series expansion."""
    term = 1.0 / factorial(k)
    value = term
    i = 1
    while True:
        term *= -z / ((k + 2 * i - 1) * (k + 2 * i))
        new_value = value + term
        if new_value == value:
            return new_value
        value = new_value
        i += 1
