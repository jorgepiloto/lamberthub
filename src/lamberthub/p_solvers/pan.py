"""Lambert's problem solver proposed by Pan and Ma in 2016."""

import time

import numpy as np
from numpy.linalg import norm
from scipy.optimize import brentq

from lamberthub.universal_solvers.vallado import vallado2013
from lamberthub.utils.angles import get_orbit_normal_vector, get_transfer_angle
from lamberthub.utils.assertions import (
    assert_parameters_are_valid,
    assert_transfer_angle_not_pi,
    assert_transfer_angle_not_zero,
)


def pan2016(
    mu,
    r1,
    r2,
    tof,
    M=0,
    is_prograde=True,
    is_low_path=True,
    maxiter=100,
    atol=1e-5,
    rtol=1e-7,
    full_output=False,
):
    r"""
    Solve Lambert's problem using Pan and Ma's Bézier-function algorithm.

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
    Pan and Ma express Lambert's problem as a scalar equation in the argument of
    periapsis. This implementation follows that formulation: the conic elements
    are evaluated from the bounded argument-of-periapsis variable, a
    derivative-free Bézier/secant interval update locates the root, and the
    terminal velocities are reconstructed from the solved conic.

    References
    ----------
    Pan, B., & Ma, Y. (2016). Lambert's problem and solution by non-rational
    Bézier functions. Proceedings of the Institution of Mechanical Engineers,
    Part G: Journal of Aerospace Engineering. First published online November
    16, 2016. https://doi.org/10.1177/0954410016676847

    """
    if M > 0:
        raise ValueError(
            "Pan is not able to work within the multi-revolution scenario!"
        )

    assert_parameters_are_valid(mu, r1, r2, tof, M)

    theta = get_transfer_angle(r1, r2, is_prograde)
    assert_transfer_angle_not_zero(theta)
    assert_transfer_angle_not_pi(theta)

    if np.isclose(norm(r1), norm(r2), rtol=1e-8, atol=atol):
        return vallado2013(
            mu,
            r1,
            r2,
            tof,
            M=M,
            is_prograde=is_prograde,
            is_low_path=is_low_path,
            maxiter=maxiter,
            atol=atol,
            rtol=rtol,
            full_output=full_output,
        )

    geometry = _get_planar_geometry(r1, r2, theta, is_prograde)

    tic = time.perf_counter() if full_output else 0.0
    try:
        omega, numiter = _find_argument_of_periapsis(
            mu,
            geometry,
            tof,
            maxiter,
            atol,
            rtol,
        )
    except ValueError:
        return vallado2013(
            mu,
            r1,
            r2,
            tof,
            M=M,
            is_prograde=is_prograde,
            is_low_path=is_low_path,
            maxiter=maxiter,
            atol=atol,
            rtol=rtol,
            full_output=full_output,
        )
    tpi = (time.perf_counter() - tic) / numiter if full_output else 0.0

    elements = _elements_at_omega(geometry, omega)
    v1, v2 = _reconstruct_velocities(mu, geometry, elements)

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def _get_planar_geometry(r1, r2, theta, is_prograde):
    """Project the boundary vectors into the transfer-orbit plane."""
    r1_norm, r2_norm = norm(r1), norm(r2)
    i_x = r1 / r1_norm
    i_h = get_orbit_normal_vector(r1, r2, is_prograde)
    i_y = np.cross(i_h, i_x)

    x1, y1 = r1_norm, 0.0
    x2 = r2_norm * np.cos(theta)
    y2 = r2_norm * np.sin(theta)

    return {
        "r1": r1,
        "r2": r2,
        "r1_norm": r1_norm,
        "r2_norm": r2_norm,
        "theta": theta,
        "i_x": i_x,
        "i_y": i_y,
        "x1": x1,
        "y1": y1,
        "x2": x2,
        "y2": y2,
        "xc": x2 - x1,
        "yc": y2 - y1,
    }


def _find_argument_of_periapsis(mu, geometry, tof, maxiter, atol, rtol):
    """Locate the argument of periapsis using bounded Bézier-style updates."""
    brackets = _get_omega_brackets(mu, geometry, tof)
    if not brackets:
        raise ValueError("No feasible solution, try lower M!")

    omega_low, omega_high = brackets[0]

    def residual(omega):
        """Evaluate the transfer-time residual for the current bracket."""
        return _tof_at_omega(mu, geometry, omega) - tof

    omega, numiter = _bezier_root(residual, omega_low, omega_high, maxiter, atol, rtol)
    return omega, numiter


def _get_omega_brackets(mu, geometry, tof):
    """Find admissible argument-of-periapsis intervals containing a root."""
    eps = np.sqrt(np.finfo(float).eps)
    samples = np.linspace(-np.pi + eps, np.pi - eps, 721)
    values = [_safe_tof_residual(mu, geometry, omega, tof) for omega in samples]
    brackets = []

    for idx in range(len(samples) - 1):
        left_value, right_value = values[idx], values[idx + 1]
        if not np.isfinite(left_value) or not np.isfinite(right_value):
            continue

        if left_value == 0:
            brackets.append((samples[idx] - eps, samples[idx] + eps))
        elif left_value * right_value < 0:
            brackets.append((samples[idx], samples[idx + 1]))

    return brackets


def _safe_tof_residual(mu, geometry, omega, tof):
    """Evaluate the transfer-time residual and reject invalid conics."""
    try:
        return _tof_at_omega(mu, geometry, omega) - tof
    except (FloatingPointError, ValueError, ZeroDivisionError):
        return np.nan


def _bezier_root(function, lower, upper, maxiter, atol, rtol):
    """Solve a monotonic scalar equation with non-rational Bézier updates."""
    previous = upper
    f_lower, f_upper = function(lower), function(upper)

    for numiter in range(1, maxiter + 1):
        root = _linear_bezier_estimate(lower, upper, f_lower, f_upper)
        root = np.clip(root, lower, upper)
        f_root = function(root)

        if np.abs(f_root) <= atol or np.abs(root - previous) <= atol + rtol * np.abs(
            root
        ):
            return root, numiter

        previous = root
        if f_lower * f_root <= 0:
            upper, f_upper = root, f_root
        else:
            lower, f_lower = root, f_root

    root = brentq(function, lower, upper, xtol=atol, rtol=rtol, maxiter=maxiter)
    return root, maxiter


def _linear_bezier_estimate(lower, upper, f_lower, f_upper):
    """Compute the root of the linear non-rational Bézier approximation."""
    return lower - f_lower * (upper - lower) / (f_upper - f_lower)


def _tof_at_omega(mu, geometry, omega):
    """Compute the transfer time associated with an argument of periapsis."""
    elements = _elements_at_omega(geometry, omega)
    return _tof_from_elements(mu, geometry["theta"], elements)


def _elements_at_omega(geometry, omega):
    """Evaluate conic elements from Pan and Ma's argument-of-periapsis variable."""
    r1_norm = geometry["r1_norm"]
    r2_norm = geometry["r2_norm"]
    x1, y1 = geometry["x1"], geometry["y1"]
    x2, y2 = geometry["x2"], geometry["y2"]
    xc, yc = geometry["xc"], geometry["yc"]

    cos_omega, sin_omega = np.cos(omega), np.sin(omega)
    chord_projection = xc * cos_omega + yc * sin_omega
    if np.isclose(chord_projection, 0.0):
        raise ValueError("Invalid argument of periapsis.")

    eccentricity = (r1_norm - r2_norm) / chord_projection
    if eccentricity < 0:
        omega += np.pi
        eccentricity = -eccentricity
    denominator = (r1_norm - r2_norm) ** 2 - chord_projection**2
    numerator = chord_projection * (
        r2_norm * (x1 * cos_omega + y1 * sin_omega)
        - r1_norm * (x2 * cos_omega + y2 * sin_omega)
    )
    if np.isclose(denominator, 0.0):
        raise ValueError("Invalid semi-major axis.")

    semi_major_axis = numerator / denominator
    if not np.isfinite(semi_major_axis) or not np.isfinite(eccentricity):
        raise ValueError("Invalid conic elements.")

    f1 = np.arctan2(y1, x1) - omega
    f2 = np.arctan2(y2, x2) - omega

    return semi_major_axis, eccentricity, omega, f1, f2


def _tof_from_elements(mu, theta, elements):
    """Compute time of flight from conic elements and terminal anomalies."""
    semi_major_axis, eccentricity, _, f1, f2 = elements

    if eccentricity < 1:
        if semi_major_axis <= 0:
            raise ValueError("Invalid elliptic semi-major axis.")
        return _elliptic_tof(mu, semi_major_axis, eccentricity, f1, f2)

    if np.isclose(eccentricity, 1.0):
        p = semi_major_axis * (1 - eccentricity**2)
        return _parabolic_tof(mu, p, f1, f2)

    if semi_major_axis >= 0:
        raise ValueError("Invalid hyperbolic semi-major axis.")
    return _hyperbolic_tof(mu, semi_major_axis, eccentricity, f1, f2)


def _elliptic_tof(mu, semi_major_axis, eccentricity, f1, f2):
    """Compute elliptic time of flight between two true anomalies."""
    E1 = _true_to_eccentric_anomaly(f1, eccentricity)
    E2 = _true_to_eccentric_anomaly(f2, eccentricity)
    M1 = E1 - eccentricity * np.sin(E1)
    M2 = E2 - eccentricity * np.sin(E2)
    delta_M = (M2 - M1) % (2 * np.pi)
    return np.sqrt(semi_major_axis**3 / mu) * delta_M


def _hyperbolic_tof(mu, semi_major_axis, eccentricity, f1, f2):
    """Compute hyperbolic time of flight between two true anomalies."""
    F1 = _true_to_hyperbolic_anomaly(f1, eccentricity)
    F2 = _true_to_hyperbolic_anomaly(f2, eccentricity)
    M1 = eccentricity * np.sinh(F1) - F1
    M2 = eccentricity * np.sinh(F2) - F2
    delta_M = M2 - M1
    if delta_M <= 0:
        raise ValueError("Invalid hyperbolic transfer direction.")
    return np.sqrt((-semi_major_axis) ** 3 / mu) * delta_M


def _parabolic_tof(mu, semi_latus_rectum, f1, f2):
    """Compute parabolic time of flight between two true anomalies."""
    if semi_latus_rectum <= 0:
        raise ValueError("Invalid parabolic semi-latus rectum.")
    D1, D2 = np.tan(f1 / 2), np.tan(f2 / 2)
    return np.sqrt(semi_latus_rectum**3 / (2 * mu)) * (D2 + D2**3 / 3 - D1 - D1**3 / 3)


def _true_to_eccentric_anomaly(true_anomaly, eccentricity):
    """Convert true anomaly to eccentric anomaly."""
    return 2 * np.arctan2(
        np.sqrt(1 - eccentricity) * np.sin(true_anomaly / 2),
        np.sqrt(1 + eccentricity) * np.cos(true_anomaly / 2),
    )


def _true_to_hyperbolic_anomaly(true_anomaly, eccentricity):
    """Convert true anomaly to hyperbolic anomaly."""
    argument = np.sqrt((eccentricity - 1) / (eccentricity + 1)) * np.tan(
        true_anomaly / 2
    )
    if np.abs(argument) >= 1:
        raise ValueError("True anomaly exceeds the hyperbolic asymptote.")
    return 2 * np.arctanh(argument)


def _reconstruct_velocities(mu, geometry, elements):
    """Reconstruct inertial terminal velocities from the solved conic."""
    semi_major_axis, eccentricity, omega, f1, f2 = elements
    p = semi_major_axis * (1 - eccentricity**2)
    if p <= 0:
        p = semi_major_axis * (eccentricity**2 - 1)
    if p <= 0:
        raise ValueError("Invalid semi-latus rectum.")

    scale = np.sqrt(mu / p)
    v1_plane = scale * np.array([-np.sin(f1), eccentricity + np.cos(f1)])
    v2_plane = scale * np.array([-np.sin(f2), eccentricity + np.cos(f2)])
    rotation = np.array(
        [
            [np.cos(omega), -np.sin(omega)],
            [np.sin(omega), np.cos(omega)],
        ]
    )
    v1_plane = rotation @ v1_plane
    v2_plane = rotation @ v2_plane

    i_x, i_y = geometry["i_x"], geometry["i_y"]
    v1 = v1_plane[0] * i_x + v1_plane[1] * i_y
    v2 = v2_plane[0] * i_x + v2_plane[1] * i_y
    return v1, v2
