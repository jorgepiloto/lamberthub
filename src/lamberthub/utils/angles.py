""" Utilities related to angles computations """

import numpy as np
from numpy import cross, dot
from numpy.linalg import norm


def get_transfer_angle(r1, r2, prograde):
    """
    Solves for the transfer angle being known the sense of rotation.

    Parameters
    ----------
    r1: np.array
        Initial position vector.
    r2: np.array
        Final position vector.
    prograde: bool
        If True, it assumes prograde motion, otherwise assumes retrograde.

    Returns
    -------
    dtheta: float
        Transfer angle in radians.

    """

    # Check if both position vectors are collinear. If so, check if the transfer
    # angle is 0 or pi.
    if np.all(np.cross(r1, r2) == 0):
        return 0 if np.all(np.sign(r1) == np.sign(r2)) else np.pi

    # Solve for a unitary vector normal to the vector plane. Its direction and
    # sense the one given by the cross product (right-hand) from r1 to r2.
    h = cross(r1, r2) / norm(np.cross(r1, r2))

    # Compute the projection of the normal vector onto the reference plane.
    alpha = dot(np.array([0, 0, 1]), h)

    # Get the minimum angle (0 <= dtheta <= pi) between r1 and r2.
    r1_norm, r2_norm = [norm(vec) for vec in [r1, r2]]
    theta0 = np.arccos(dot(r1, r2) / (r1_norm * r2_norm))

    # Fix the value of theta if necessary
    if prograde is True:
        dtheta = theta0 if alpha > 0 else 2 * np.pi - theta0
    else:
        dtheta = theta0 if alpha < 0 else 2 * np.pi - theta0

    return dtheta


def get_orbit_normal_vector(r1, r2, prograde):
    """
    Computes a unitary normal vector aligned with the specific angular momentum
    one of the orbit.

    Parameters
    ----------
    r1: np.array
        Initial position vector.
    r2: np.array
        Final position vector.
    prograde: bool
        If True, it assumes prograde motion, otherwise assumes retrograde.

    Returns
    -------
    i_h: np.array
        Unitary vector aligned with orbit specific angular momentum.

    """

    # Compute the normal vector and its projection onto the vertical axis
    i_h = np.cross(r1, r2) / norm(np.cross(r1, r2))

    # Solve the projection onto the positive vertical direction of the
    # fundamental plane.
    alpha = dot(np.array([0, 0, 1]), i_h)

    # An prograde orbit always has a positive vertical component of its specific
    # angular momentum. Therefore, we just need to check for this condition
    if prograde is True:
        i_h = i_h if alpha > 0 else -i_h
    else:
        i_h = i_h if alpha < 0 else -i_h

    return i_h
