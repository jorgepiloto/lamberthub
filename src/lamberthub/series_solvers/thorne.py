"""This module holds all methods devised by James D. Thorne."""

import time

import numpy as np
from numpy.linalg import norm
from scipy.special import factorial, poch

from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import assert_parameters_are_valid


def thorne2004(
    mu,
    r1,
    r2,
    tof,
    M=0,
    prograde=True,
    low_path=True,
    maxiter=60,
    atol=1e-5,
    rtol=1e-7,
    full_output=False,
):
    r"""
    Lambert's problem solver devised by Carl Friedrich Gauss in 1809. The method
    has been implemented according to Bate's book (see [2]) and extended to the
    hyperbolic case. This method shows poor accuracy, being only suitable for
    low transfer angles.

    Parameters
    ----------
    mu: float
        Gravitational parameter, equivalent to :math:`GM` of attractor body.
    r1: numpy.array
        Initial position vector.
    r2: numpy.array
        Final position vector.
    M: int
        Number of revolutions. Must be equal or greater than 0 value.
    prograde: bool
        If `True`, specifies prograde motion. Otherwise, retrograde motion is imposed.
    low_path: bool
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
    The algorithm originally devised by Gauss exploits the so-called ratio of
    sector to triangle area, which is a numerical value related with the orbital
    parameter. This Algorithm was used to the discovery of the orbit of Ceres by
    the genius and adopted by many other authors of his time due to its
    simplicity. However, the Algorithm is found to be singular for transfer
    angles of 180 degrees and shows a low performance for really small angles.

    References
    ----------
    [1] Thorne, J. D. (2004). Lambert’s theorem—a complete series solution. The
    Journal of the Astronautical Sciences, 52(4), 441-454.

    [2] THORNE, J. (1995, July). Series reversion/inversion of Lambert's time
    function. In Astrodynamics Conference (p. 2886).

    """

    # Check that input parameters are safe
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    # Norm of the initial and final position vectors
    r1_norm, r2_norm, c_norm = [norm(r) for r in [r1, r2, r2 - r1]]
    s = (r1_norm + r2_norm + c_norm) / 2

    # Compute the cosine of the transfer angle and check
    dtheta = get_transfer_angle(r1, r2, prograde)

    # Solve for the parabolic transfer time given by equation (10) from the
    # original report [1]. If the tr angle is greater than 180 degrees, a sign
    # correction needs to be applied
    mp = -1 if dtheta < np.pi else 1
    s_minus_c_over_c = (s - c_norm) / c_norm
    t_p = (
        (np.sqrt(2) / 3)
        * np.sqrt(s ** 3 / mu)
        * (1 + mp * (s_minus_c_over_c) ** (3 / 2))
    )
    t_min = np.sqrt(s ** 3 / (8 * mu)) * (
        np.pi
        + mp
        * (
            2 * np.arcsin(np.sqrt(s_minus_c_over_c))
            - np.sin(2 * np.arcsin(s_minus_c_over_c))
        )
        ** (3 / 2)
    )

    # The non-dimensional transfer time can be obtained using expression (11)
    # from Thorne's article [1]. If the current time of flight is greater than
    # the parabolic one, then the value of T > 0 Otherwise is T < 0. This can be
    # used to apply a particular series solution as the shape of the final orbit
    # can be easily predicted.
    T = tof / t_p - 1
    T_min = tof / t_min - 1

    # Compute the A_array coefficients, Q_matrix and B
    A_array = get_A_array(s, c_norm, maxiter, mp)
    Q_matrix = get_Q_matrix(A_array, maxiter)
    B_array = get_B_array(A_array, Q_matrix, maxiter)

    # Filter out the series expansion to be applied by direct comparison with
    # the current time of flight.

    if (T < 0) or (0 < T < T_min):
        # Hyperbolic case and short elliptic [H, A]
        # a = (s / 2) * np.sum([B_array[i] * T ** (i - 1) for i in range(0, maxiter)])
        print(f"Applying [H,A]")
        a = (s / 2) * np.sum([B_array[i] * T ** (i - 1) for i in range(maxiter)])
    elif T > T_min:
        print(f"Applying [B_inf]")
        a = np.sum([B_array[n] * tof ** ((2 - n) / 3) for n in range(maxiter)])
    else:
        raise ValueError("Not suitable series found!")

    print(f"Computed {a = :.3f}")

    return np.zeros(3), np.zeros(3)


def get_A_array(s, c, maxiter, mp):
    """Computes the values for the A_n coefficients.

    Parameters
    ----------
    s: float
        The semi-permieter.
    c: float
        The norm of the chord vector
    maxiter: int
        Maximum number of coefficients
    mp: float
        Minus/plus sign acording to transfer angle lower or greater than pi.

    Returns
    -------
    A_array: np.array
        An array holding the values for each one of the coefficients.

    Notes
    -----
    The expression is not given in explicit form in neither of the reports [1]
    or [2]. In fact, it is obtained by direct comparison between expressions
    (10) and (11) from original report [1].

    """

    # Allocate the array of coefficients
    A_array = np.zeros(maxiter)

    for n in range(1, maxiter + 1):
        k = (s - c) / s
        numerator = (1 + mp * k ** (n + 3 / 2)) * poch(1 / 2, n) * poch(3 / 2, n)
        denominator = (1 + mp * k ** (3 / 2)) * poch(5 / 2, n) * factorial(n)
        A_array[n - 1] = numerator / denominator

    return A_array


def get_Q_matrix(A_array, maxiter):
    """Computes the Q matrix being given the A_n coefficients.

    Parameters
    ----------
    A_array: np.array
        The array of A_n coefficients.
    maxiter: int
        The maximum number of elements of the series.

    Returns
    -------
    Q_matrix: np.array
        The matrix of Q_n coefficients.

    Notes
    -----
    Expressions (26), (27), (28) and (29) from original report [2] have been
    applied here.

    """

    # Allocate the matrix
    Q_matrix = np.zeros((maxiter, maxiter))

    # Only the lower triangle elements need to be computed, the rest will be
    # zero as the Q matrix is an upper triangular one. The iteration loop needs
    # to be performed straight forward for each one of the rows but reversed in
    # columns as values j depend on j+1.

    for i in range(1, maxiter + 1):
        for j in reversed(range(1, i + 1)):

            # Check if the coefficient to be evaluated is Q(1,1). If so, apply
            # the particular expression, that is (27) from [2]. Be careful about
            # index, here (1,1) is equivalent to (0,0).
            if i == 1 and j == 1:
                Q_matrix[i - 1, j - 1] = A_array[0] ** -1

            # Check if the coefficient belongs to the first column Q(i,1). This
            # is expression (29) from the original report [2]. Again, be careful
            # with index notation, as Q(i,1) is equivalent to Q(i,0).
            elif i != 1 and j == 1:
                Q_matrix[i - 1, j - 1] = np.sum(
                    [
                        (-1 / A_array[0]) * Q_matrix[i - 1, k] * A_array[k]
                        for k in range(1, i)
                    ]
                )
            # If none of previous conditions holds, then apply expression (28)
            # from [2] to evaluate the particular value of the coefficient.
            else:
                Q_matrix[i - 1, j - 1] = np.sum(
                    [
                        Q_matrix[i - k - 1, j - 1 - 1] * Q_matrix[k - 1, 0]
                        for k in range(1, i)
                    ]
                )

    return Q_matrix


def get_B_array(A_array, Q_matrix, maxiter):

    B_array = np.zeros(maxiter)

    B_array[0] = A_array[0]

    for n in range(2, maxiter + 1):
        B_array[n - 1] = np.sum(
            [Q_matrix[n - 1 - 1, m - 1] * A_array[m] for m in range(1, n)]
        )

    return B_array
