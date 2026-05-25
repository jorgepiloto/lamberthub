"""A module hosting the contour-integral solver by Negrete and Abdelkhalik."""

from math import factorial
import time

from numba import njit as jit
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
    tic = time.perf_counter() if full_output else 0.0
    z = _find_z(mu, r1_norm, r2_norm, A, tof, M, is_low_path, numiter, atol, rtol)
    tpi = (time.perf_counter() - tic) / (3 * numiter) if full_output else 0.0

    y = _y_at_z(np.complex128(z), r1_norm, r2_norm, A).real
    f = 1.0 - y / r1_norm
    g = A * np.sqrt(y / mu)
    gdot = 1.0 - y / r2_norm

    v1 = (r2 - f * r1) / g
    v2 = (gdot * r2 - r1) / g

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


@jit(cache=True, fastmath=True)
def _find_z(mu, r1_norm, r2_norm, A, tof, M, is_low_path, numiter, atol, rtol):
    """Find the universal anomaly root from contour-integral moments."""
    if M == 0:
        tp = _tof_at_z(0.0, mu, r1_norm, r2_norm, A)
        if tof < tp:
            center, radius, ok = _hyperbolic_contour(mu, r1_norm, r2_norm, A, tof)
            if not ok:
                raise ValueError("Could not bracket hyperbolic contour!")
            z_min, z_max = center - radius, 0.0
        else:
            center = 2.0 * np.pi**2
            radius = 2.0 * np.pi**2 - 0.5
            z_min, z_max = 0.0, 4.0 * np.pi**2

        moments = _contour_integrals(
            center, radius, mu, r1_norm, r2_norm, A, tof, 1, numiter
        )
        I0 = moments[0]
        I1 = moments[1]
        z = _real_solution(I1 / I0, atol, rtol)
        return _refine_single_root(
            z, z_min, z_max, mu, r1_norm, r2_norm, A, tof, numiter, atol, rtol
        )

    z_min = (2.0 * M * np.pi) ** 2
    z_max = (2.0 * (M + 1) * np.pi) ** 2
    center = 2.0 * np.pi**2 * (2 * M**2 + 2 * M + 1)
    radius = 2.0 * np.pi**2 * (2 * M + 1) - 0.5
    moments = _contour_integrals(
        center, radius, mu, r1_norm, r2_norm, A, tof, 3, numiter
    )
    I0 = moments[0]
    I1 = moments[1]
    I2 = moments[2]
    I3 = moments[3]

    a = I1**2 - I0 * I2
    b = I3 * I0 - I1 * I2
    c = I2**2 - I1 * I3

    # Closed-form quadratic; coefficients are complex.
    disc = b * b - 4.0 * a * c
    sqrt_disc = np.sqrt(disc)
    root1 = (-b + sqrt_disc) / (2.0 * a)
    root2 = (-b - sqrt_disc) / (2.0 * a)

    if abs(root1.imag) > atol + rtol * abs(root1.real) or abs(
        root2.imag
    ) > atol + rtol * abs(root2.real):
        raise ValueError("No feasible solution, try lower M!")

    r1r = root1.real
    r2r = root2.real
    if r1r > r2r:
        r1r, r2r = r2r, r1r

    refined_lo = _refine_single_root(
        r1r, z_min, z_max, mu, r1_norm, r2_norm, A, tof, numiter, atol, rtol
    )
    refined_hi = _refine_single_root(
        r2r, z_min, z_max, mu, r1_norm, r2_norm, A, tof, numiter, atol, rtol
    )
    if refined_lo > refined_hi:
        refined_lo, refined_hi = refined_hi, refined_lo
    return refined_lo if is_low_path else refined_hi


@jit(cache=True, fastmath=True)
def _refine_single_root(
    z, z_min, z_max, mu, r1_norm, r2_norm, A, tof, numiter, atol, rtol
):
    """Refine a candidate root using smaller local contours."""
    for _ in range(2):
        center = _refinement_center(z, z_min, z_max)
        radius = _refinement_radius(center, z_min, z_max)
        moments = _contour_integrals(
            center, radius, mu, r1_norm, r2_norm, A, tof, 1, numiter
        )
        I0 = moments[0]
        I1 = moments[1]
        z = _real_solution(I1 / I0, atol, rtol)
    return z


@jit(cache=True, fastmath=True)
def _refinement_center(z, z_min, z_max):
    """Choose the center for a root-refinement contour."""
    if z_min == 0.0 and z < z_min:
        v = 0.25 * z_max
        return 0.5 if 0.5 < v else v
    return z


@jit(cache=True, fastmath=True)
def _refinement_radius(z, z_min, z_max):
    """Choose a bounded contour radius around a candidate root."""
    radius = 2.0
    if z_min == 0.0 and 0.0 < z < 1.0:
        a = z
        b = 0.5 * (z_max - z)
        if a < radius:
            radius = a
        if b < radius:
            radius = b
    elif z_min < z < z_max:
        a = 0.5 * (z - z_min)
        b = 0.5 * (z_max - z)
        if a < radius:
            radius = a
        if b < radius:
            radius = b
    elif z_min < z < 0.0 < z_max:
        a = 0.9 * abs(z)
        if a < radius:
            radius = a
    return radius if radius > 1e-8 else 1e-8


@jit(cache=True, fastmath=True)
def _contour_integrals(center, radius, mu, r1_norm, r2_norm, A, tof, degree, numiter):
    """Compute contour-integral moments around a circular path."""
    out = np.zeros(degree + 1, dtype=np.complex128)
    two_pi_j = 2j * np.pi
    inv_n = 1.0 / numiter
    for i in range(numiter):
        t = i * inv_n
        exp_t = np.exp(two_pi_j * t)
        z = center + radius * exp_t
        dzdx = two_pi_j * radius * exp_t
        fz = _pole_function(z, mu, r1_norm, r2_norm, A, tof)
        w = fz * dzdx
        zk = 1.0 + 0j
        for k in range(degree + 1):
            out[k] += zk * w
            zk = zk * z
    for k in range(degree + 1):
        out[k] *= inv_n
    return out


@jit(cache=True, fastmath=True)
def _pole_function(z, mu, r1_norm, r2_norm, A, tof):
    """Evaluate the reciprocal time-of-flight residual in the complex plane."""
    y = _y_at_z(z, r1_norm, r2_norm, A)
    chi = np.sqrt(y / _c2(z))
    return 1.0 / (chi**3 * _c3(z) + A * np.sqrt(y) - np.sqrt(mu) * tof)


@jit(cache=True, fastmath=True)
def _tof_at_z(z, mu, r1_norm, r2_norm, A):
    """Evaluate the universal-variable time of flight at a z value."""
    zc = np.complex128(z)
    y = _y_at_z(zc, r1_norm, r2_norm, A)
    chi = np.sqrt(y / _c2(zc))
    return ((chi**3 * _c3(zc) + A * np.sqrt(y)) / np.sqrt(mu)).real


@jit(cache=True, fastmath=True)
def _hyperbolic_contour(mu, r1_norm, r2_norm, A, tof):
    """Build a negative-z contour for hyperbolic single-revolution transfers."""
    z_low = -1.0
    success = False
    for _ in range(64):
        y = _y_at_z(np.complex128(z_low), r1_norm, r2_norm, A).real
        if y > 0.0 and _tof_at_z(z_low, mu, r1_norm, r2_norm, A) < tof:
            success = True
            break
        z_low *= 2.0

    center = 0.5 * z_low
    radius = abs(center) * (1.0 - 1e-6)
    return center, radius, success


@jit(cache=True, fastmath=True)
def _real_solution(z, atol, rtol):
    """Return a contour root if its imaginary component is negligible."""
    if abs(z.imag) > atol + rtol * abs(z.real):
        raise ValueError("Contour integral did not return a real solution!")
    return z.real


@jit(cache=True, fastmath=True)
def _get_A(r1_norm, r2_norm, dtheta):
    """Compute Vallado's geometry parameter for the transfer angle."""
    t_m = 1.0 if dtheta < np.pi else -1.0
    return t_m * np.sqrt(r1_norm * r2_norm * (1 + np.cos(dtheta)))


@jit(cache=True, fastmath=True)
def _y_at_z(z, r1_norm, r2_norm, A):
    """Compute the universal-variable y parameter at a z value."""
    return r1_norm + r2_norm + A * (z * _c3(z) - 1.0) / np.sqrt(_c2(z))


@jit(cache=True, fastmath=True)
def _c2(z):
    """Evaluate the second Stumpff function for a complex scalar z."""
    if abs(z) < 1e-8:
        return _stumpff_series(z, 2)
    return (1.0 - np.cos(np.sqrt(z))) / z


@jit(cache=True, fastmath=True)
def _c3(z):
    """Evaluate the third Stumpff function for a complex scalar z."""
    if abs(z) < 1e-8:
        return _stumpff_series(z, 3)
    sqrt_z = np.sqrt(z)
    return (sqrt_z - np.sin(sqrt_z)) / (z * sqrt_z)


# Precompute factorial reciprocals at import time (order=2 and order=3 series).
_FACT_INV_C2 = np.array(
    [1.0 / float(factorial(2 * k + 2)) for k in range(20)], dtype=np.float64
)
_FACT_INV_C3 = np.array(
    [1.0 / float(factorial(2 * k + 3)) for k in range(20)], dtype=np.float64
)


@jit(cache=True, fastmath=True)
def _stumpff_series(z, order):
    """Evaluate a Stumpff function series for small z values (complex scalar)."""
    if order == 2:
        coeffs = _FACT_INV_C2
    else:
        coeffs = _FACT_INV_C3
    total = np.complex128(coeffs[0])
    neg_z = -z
    term = np.complex128(1.0)
    for k in range(1, 20):
        term = term * neg_z
        total = total + term * coeffs[k]
    return total
