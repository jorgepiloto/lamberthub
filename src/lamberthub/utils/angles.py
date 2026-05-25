"""Utilities related to angles computations"""

from numba import njit as jit
import numpy as np

from lamberthub.linalg import cross, dot, norm

K_AXIS = np.array([0.0, 0.0, 1.0])
"""Unitary vector along the Z-axis of the fundamental plane."""


@jit(cache=True, fastmath=True)
def get_transfer_angle(r1, r2, is_prograde):
    """Compute the transfer angle of the trajectory.

    Initial and final position vectors are required together with the direction
    of motion.

    Parameters
    ----------
    r1 : ~np.array
        Initial position vector.
    r2 : ~np.array
        Final position vector.
    is_prograde : bool
        ``True`` for prograde motion, ``False`` otherwise.

    Returns
    -------
    dtheta : float
        Transfer angle in radians.

    """
    # Check if both position vectors are collinear. If so, check if the transfer
    # angle is 0 or pi.
    cross_r1r2 = cross(r1, r2)
    if np.all(cross_r1r2 == 0):
        return 0 if np.all(np.sign(r1) == np.sign(r2)) else np.pi

    # Solve for a unitary vector normal to the vector plane. Its direction and
    # sense the one given by the cross product (right-hand) from r1 to r2.
    h = cross_r1r2 / norm(cross_r1r2)

    # Compute the projection of the normal vector onto the reference plane.
    alpha = dot(K_AXIS, h)

    # Get the minimum angle (0 <= dtheta <= pi) between r1 and r2.
    r1_norm = norm(r1)
    r2_norm = norm(r2)
    theta0 = np.arccos(dot(r1, r2) / (r1_norm * r2_norm))

    # Fix the value of theta if necessary
    if is_prograde is True:
        dtheta = theta0 if alpha > 0 else 2 * np.pi - theta0
    else:
        dtheta = theta0 if alpha < 0 else 2 * np.pi - theta0

    return dtheta


@jit(cache=True, fastmath=True)
def get_orbit_normal_vector(r1, r2, is_prograde):
    """
    Computes a unitary normal vector aligned with the specific angular momentum
    one of the orbit.

    Parameters
    ----------
    r1: np.array
        Initial position vector.
    r2: np.array
        Final position vector.
    is_prograde: bool
        If True, it assumes prograde motion, otherwise assumes retrograde.

    Returns
    -------
    i_h: np.array
        Unitary vector aligned with orbit specific angular momentum.

    """
    # Compute the normal vector and its projection onto the vertical axis
    cross_r1r2 = cross(r1, r2)
    i_h = cross_r1r2 / norm(cross_r1r2)

    # Solve the projection onto the positive vertical direction of the
    # fundamental plane.
    alpha = dot(K_AXIS, i_h)

    # An prograde orbit always has a positive vertical component of its specific
    # angular momentum. Therefore, we just need to check for this condition
    if is_prograde is True:
        i_h = i_h if alpha > 0 else -i_h
    else:
        i_h = i_h if alpha < 0 else -i_h

    return i_h


@jit(cache=True, fastmath=True)
def get_orbit_inc_and_raan_from_position_vectors(r1, r2, is_prograde):
    """
    Computes the inclination of the orbit being known an initial and a final
    position vectors together with the sense of motion.

    Parameters
    ----------
    r1: np.array
        Initial position vector.
    r2: np.array
        Final position vector.
    is_prograde: bool
        If True, it assumes prograde motion, otherwise assumes retrograde.

    Returns
    -------
    inc: float
        Inclination of the orbit.
    raan: float
        Right ascension of the ascending node.

    """
    # Get a unitary vector aligned in direction and sense with the specific
    # angular momentum one.
    i_h = get_orbit_normal_vector(r1, r2, is_prograde)

    # If the orbit is coplanar with fundamental plane, neither inc or raan are
    # defined. TODO: use atol and rtol instead of pure floating zero comparison
    if i_h[0] == 0 and i_h[1] == 0:
        inc, raan = 0, 0
    else:
        # Inclination is always bounded between [0, pi], so no correction is
        # needed
        inc = np.arccos(i_h[2] / norm(i_h))

        # Compute the RAAN using a vector in the direction and sense of the line
        # of nodes. Because RAAN is bounded between [0, 2pi], the arctan2
        # function is used.
        n = cross(K_AXIS, i_h)
        raan = np.arctan2(n[1], n[0]) % (2 * np.pi)

    return inc, raan


@jit(cache=True, fastmath=True)
def nu_to_E(nu, ecc):
    """
    Retrieves eccentric anomaly from true one.

    Parameters
    ----------
    nu: float
        True anomaly.
    ecc: float
        Eccentricity of the orbit.

    Returns
    -------
    E: float
        Eccentric anomaly.

    """
    E = 2 * np.arctan(np.sqrt((1 - ecc) / (1 + ecc)) * np.tan(nu / 2))
    return E


@jit(cache=True, fastmath=True)
def E_to_nu(E, ecc):
    """
    Retrieves true anomaly from eccentric one.

    Parameters
    ----------
    E: float
        Eccentric anomaly.
    ecc: float
        Eccentricity of the orbit.

    Returns
    -------
    nu: float
        True anomaly.

    """
    nu = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(E / 2))
    return nu


@jit(cache=True, fastmath=True)
def nu_to_B(nu):
    """
    Retrieves parabolic anomaly from true one.

    Parameters
    ----------
    nu: float
        True anomaly

    Returns
    -------
    B: float
        Parabolic anomaly

    Notes
    -----
    As explained in Vallado's [1], :math:`B` is used instead of :math:`P` just
    to not confuse with the orbital parameter.

    """
    B = np.tan(nu / 2)
    return B


@jit(cache=True, fastmath=True)
def B_to_nu(B):
    """
    Retrieves the true anomaly from parabolic one.

    Parameters
    ----------
    B: float
        Parabolic anomaly

    Returns
    -------
    nu: float
        True anomaly

    Notes
    -----
    As explained in Vallado's [1], :math:`B` is used instead of :math:`P` just
    to not confuse with the orbital parameter.

    """
    nu = 2 * np.arctan(B)
    return nu


@jit(cache=True, fastmath=True)
def nu_to_H(nu, ecc):
    """
    Retrieves hyperbolic anomaly from true one.

    Parameters
    ----------
    nu: float
        True anomaly
    ecc: float
        Eccentricity of the orbit

    Returns
    -------
    H: float
        Hyperbolic anomaly

    """
    H = 2 * np.arctanh(np.sqrt((ecc - 1) / (ecc + 1)) * np.tan(nu / 2))
    return H


@jit(cache=True, fastmath=True)
def H_to_nu(H, ecc):
    """
    Retrieves hyperbolic anomaly from true one.

    Parameters
    ----------
    H: float
        Hyperbolic anomaly
    ecc: float
        Eccentricity of the orbit

    Returns
    -------
    nu: float
        True anomaly

    """
    nu = 2 * np.arctan(np.sqrt((ecc + 1) / (ecc - 1)) * np.tanh(H / 2))
    return nu
