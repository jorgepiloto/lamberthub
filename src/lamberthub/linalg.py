"""Module containing various routines focused on linear algebra."""

from numba import njit as jit
import numpy as np


@jit
def dot(v1, v2):
    """Compute the dot product between two vectors.

    Parameters
    ----------
    v1 : ~np.array
        First vector.
    v2 : ~np.array
        Second vector.

    Returns
    -------
    float
        Magnitude of the vector.

    Notes
    -----
    This function casts the type of coordinates to double (float64) to avoid
    this issue: https://github.com/numba/numba/issues/8676

    """
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]


@jit
def cross(v1, v2):
    """Compute the cross product between two vectors.

    Parameters
    ----------
    v1 : ~np.array
        First vector.
    v2 : ~np.array
        Second vector.

    Returns
    -------
    ~np.array
        Resultant vector of the cross product between the two vectors.

    """
    return np.cross(v1, v2)


@jit
def norm(vector):
    """Compute the magnitude of a vector.

    Parameters
    ----------
    vector : ~np.array
        Vector whose magnitude is to be computed.

    Returns
    -------
    float
        Magnitude of the vector.

    """
    return dot(vector, vector) ** 0.5
