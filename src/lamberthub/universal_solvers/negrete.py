"""A module hosting the contour-integral solver by Negrete and Abdelkhalik."""

from math import factorial
import time

import numpy as np

from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import (
    assert_parameters_are_valid,
    assert_transfer_angle_not_zero,
)


def negrete2024(
    mu,
    r1,
    r2,
    tof,
    M=0,
    is_prograde=True,
    is_low_path=True,
    maxiter=64,
    atol=1e-5,
    rtol=1e-7,
    full_output=False,
):
    r"""
    Solve Lambert's problem using the contour-integral method.

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
        Number of trapezoidal quadrature points used in each contour integral.
    atol: float
        Absolute tolerance used when checking if a multi-revolution root is real.
    rtol: float
        Relative tolerance used when checking if a multi-revolution root is real.
    full_output: bool
        If True, the number of quadrature points and time per point are also returned.

    Returns
    -------
    v1: numpy.array
        Initial velocity vector.
    v2: numpy.array
        Final velocity vector.
    numiter: int
        Number of quadrature points.
    tpi: float
        Time per quadrature point in seconds.

    Notes
    -----
    This is the exact solution proposed by Negrete and Abdelkhalik, where the
    root of the universal-variable time-of-flight equation is recovered from
    contour integrals in the complex plane. The integrals are approximated with
    the composite trapezoidal rule.

    References
    ----------
    [1] Negrete, A., & Abdelkhalik, O. An Exact Solution to Lambert's Problem
        Using Contour Integrals.
    """
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)
    dtheta = get_transfer_angle(r1, r2, is_prograde)
    assert_transfer_angle_not_zero(dtheta)

    A = _get_A(r1_norm, r2_norm, dtheta)
    if A == 0.0:
        raise RuntimeError("Cannot compute orbit, phase angle is 180 degrees")

    numiter = max(4, int(maxiter))
    tic = time.perf_counter()
    z = _find_z(mu, r1_norm, r2_norm, A, tof, M, is_low_path, numiter, atol, rtol)
    tac = time.perf_counter()

    y = _y_at_z(z, r1_norm, r2_norm, A).real
    f = 1.0 - y / r1_norm
    g = A * np.sqrt(y / mu)
    gdot = 1.0 - y / r2_norm

    v1 = (r2 - f * r1) / g
    v2 = (gdot * r2 - r1) / g
    tpi = (tac - tic) / (3 * numiter)

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def _find_z(mu, r1_norm, r2_norm, A, tof, M, is_low_path, numiter, atol, rtol):
    """Find the universal anomaly root from contour-integral moments."""
    if M == 0:
        tp = _tof_at_z(0.0, mu, r1_norm, r2_norm, A)
        if tof < tp:
            center, radius = _hyperbolic_contour(mu, r1_norm, r2_norm, A, tof)
            z_min, z_max = center - radius, 0.0
        else:
            center = 2.0 * np.pi**2
            radius = 2.0 * np.pi**2 - 0.5
            z_min, z_max = 0.0, 4.0 * np.pi**2

        I0, I1 = _contour_integrals(
            center, radius, mu, r1_norm, r2_norm, A, tof, 1, numiter
        )
        z = _real_solution(I1 / I0, atol, rtol)
        return _refine_single_root(
            z, z_min, z_max, mu, r1_norm, r2_norm, A, tof, numiter, atol, rtol
        )

    z_min = (2.0 * M * np.pi) ** 2
    z_max = (2.0 * (M + 1) * np.pi) ** 2
    center = 2.0 * np.pi**2 * (2 * M**2 + 2 * M + 1)
    radius = 2.0 * np.pi**2 * (2 * M + 1) - 0.5
    I0, I1, I2, I3 = _contour_integrals(
        center, radius, mu, r1_norm, r2_norm, A, tof, 3, numiter
    )

    a = I1**2 - I0 * I2
    b = I3 * I0 - I1 * I2
    c = I2**2 - I1 * I3
    roots = np.roots(np.array([a, b, c], dtype=np.complex128))

    if np.any(np.abs(roots.imag) > atol + rtol * np.abs(roots.real)):
        raise ValueError("No feasible solution, try lower M!")

    roots = np.sort(roots.real)
    roots = np.array(
        [
            _refine_single_root(
                root, z_min, z_max, mu, r1_norm, r2_norm, A, tof, numiter, atol, rtol
            )
            for root in roots
        ]
    )
    roots.sort()
    return roots[0] if is_low_path is True else roots[1]


def _refine_single_root(
    z, z_min, z_max, mu, r1_norm, r2_norm, A, tof, numiter, atol, rtol
):
    """Refine a candidate root using smaller local contours."""
    for _ in range(2):
        center = _refinement_center(z, z_min, z_max)
        radius = _refinement_radius(center, z_min, z_max)
        I0, I1 = _contour_integrals(
            center, radius, mu, r1_norm, r2_norm, A, tof, 1, numiter
        )
        z = _real_solution(I1 / I0, atol, rtol)
    return z


def _refinement_center(z, z_min, z_max):
    """Choose the center for a root-refinement contour."""
    if z_min == 0.0 and z < z_min:
        return min(0.5, 0.25 * z_max)
    return z


def _refinement_radius(z, z_min, z_max):
    """Choose a bounded contour radius around a candidate root."""
    radius = 2.0
    if z_min == 0.0 and 0.0 < z < 1.0:
        radius = min(radius, z, 0.5 * (z_max - z))
    elif z_min < z < z_max:
        radius = min(radius, 0.5 * (z - z_min), 0.5 * (z_max - z))
    elif z_min < z < 0.0 < z_max:
        radius = min(radius, 0.9 * abs(z))
    return max(radius, 1e-8)


def _contour_integrals(center, radius, mu, r1_norm, r2_norm, A, tof, degree, numiter):
    """Compute contour-integral moments around a circular path."""
    x = np.arange(numiter) / numiter
    exp = np.exp(2j * np.pi * x)
    z = center + radius * exp
    dzdx = 2j * np.pi * radius * exp
    fz = _pole_function(z, mu, r1_norm, r2_norm, A, tof)

    integrals = []
    for k in range(degree + 1):
        integrals.append(np.mean(z**k * fz * dzdx))
    return integrals


def _pole_function(z, mu, r1_norm, r2_norm, A, tof):
    """Evaluate the reciprocal time-of-flight residual in the complex plane."""
    y = _y_at_z(z, r1_norm, r2_norm, A)
    chi = np.sqrt(y / _c2(z))
    return 1.0 / (chi**3 * _c3(z) + A * np.sqrt(y) - np.sqrt(mu) * tof)


def _tof_at_z(z, mu, r1_norm, r2_norm, A):
    """Evaluate the universal-variable time of flight at a z value."""
    y = _y_at_z(z, r1_norm, r2_norm, A)
    chi = np.sqrt(y / _c2(z))
    return ((chi**3 * _c3(z) + A * np.sqrt(y)) / np.sqrt(mu)).real


def _hyperbolic_contour(mu, r1_norm, r2_norm, A, tof):
    """Build a negative-z contour for hyperbolic single-revolution transfers."""
    z_low = -1.0
    for _ in range(64):
        y = _y_at_z(z_low, r1_norm, r2_norm, A).real
        if y > 0.0 and _tof_at_z(z_low, mu, r1_norm, r2_norm, A) < tof:
            break
        z_low *= 2.0
    else:
        raise ValueError("Could not bracket hyperbolic contour!")

    center = 0.5 * z_low
    radius = abs(center) * (1.0 - 1e-6)
    return center, radius


def _real_solution(z, atol, rtol):
    """Return a contour root if its imaginary component is negligible."""
    if abs(z.imag) > atol + rtol * abs(z.real):
        raise ValueError("Contour integral did not return a real solution!")
    return z.real


def _get_A(r1_norm, r2_norm, dtheta):
    """Compute Vallado's geometry parameter for the transfer angle."""
    t_m = 1 if dtheta < np.pi else -1
    return t_m * np.sqrt(r1_norm * r2_norm * (1 + np.cos(dtheta)))


def _y_at_z(z, r1_norm, r2_norm, A):
    """Compute the universal-variable y parameter at a z value."""
    return r1_norm + r2_norm + A * (z * _c3(z) - 1.0) / np.sqrt(_c2(z))


def _c2(z):
    """Evaluate the second Stumpff function for real or complex z values."""
    z = np.asarray(z, dtype=np.complex128)
    out = np.empty_like(z)
    mask = np.abs(z) < 1e-8
    out[~mask] = (1.0 - np.cos(np.sqrt(z[~mask]))) / z[~mask]
    out[mask] = _stumpff_series(z[mask], 2)
    return out.item() if out.ndim == 0 else out


def _c3(z):
    """Evaluate the third Stumpff function for real or complex z values."""
    z = np.asarray(z, dtype=np.complex128)
    out = np.empty_like(z)
    mask = np.abs(z) < 1e-8
    sqrt_z = np.sqrt(z[~mask])
    out[~mask] = (sqrt_z - np.sin(sqrt_z)) / (z[~mask] * sqrt_z)
    out[mask] = _stumpff_series(z[mask], 3)
    return out.item() if out.ndim == 0 else out


def _stumpff_series(z, order):
    """Evaluate a Stumpff function series for small z values."""
    term = np.ones_like(z, dtype=np.complex128) / factorial(order)
    total = term.copy()
    for k in range(1, 20):
        term = (-z) ** k / factorial(2 * k + order)
        total += term
    return total
