"""Miscellany utilities"""

import numpy as np

from lamberthub.utils.elements import rotation_matrix


def get_solver_name(solver):
    """
    Returns solver's name.

    Parameters
    ----------
    solver: function
        Solver function.

    Returns
    -------
    name: str
        String representation for the solver.

    """
    # Get the function name and its length
    raw_name = str(solver.__name__.capitalize())
    len_name = len(raw_name)

    # Break into author name and implementation year
    author_name = raw_name[0 : (len_name - 4)]
    year_name = raw_name[(len_name - 4) :]

    # Assemble and return the solver's name
    name = author_name + " " + year_name
    return name


def _get_sample_vectors_from_theta_and_rho(theta, rho):
    """Returns the initial and final position vectors contained in the reference
    XY plane being given the transfer angle and the norm ration between them.

    Parameters
    ----------
    theta: float
        The transfer angle.
    rho: float
        The ratio between the norms of the final and initial vectors.

    Returns
    -------
    r1_vec: ~np.array
        The initial position vector.
    r2_vec: ~np.array
        The final position vector.

    Notes
    -----
    The purpose of this function is to easily generate initial data for the
    Lambert's problem in the sense of position vectors. The function is used by
    the performance plotters and to generate various time of flight curves for
    the different available algorithms.

    """
    # Generate final position vector by rotating initial one a given theta
    r1_vec = np.array([1, 0, 0])
    r2_vec = rho * rotation_matrix(theta, axis=2) @ r1_vec
    return (r1_vec, r2_vec)
