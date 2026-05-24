"""This module holds all methods devised by James D. Thorne."""

import time

import numpy as np
from scipy.optimize import brentq
from scipy.special import factorial, poch

from lamberthub.linalg import cross, norm
from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import assert_parameters_are_valid


def thorne2004(
    mu,
    r1,
    r2,
    tof,
    M=0,
    is_prograde=True,
    is_low_path=True,
    maxiter=60,
    atol=1e-5,
    rtol=1e-7,
    full_output=False,
):
    r"""
    Lambert's problem solver devised by James D. Thorne in 2004.

    The method classifies the zero-revolution transfer into Thorne's hyperbolic
    (H), short-way elliptic (A), or long-way elliptic (B) regions and solves the
    corresponding Lagrange time equation for the semimajor axis.

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
    Thorne's paper expresses the Lagrange time equations as direct series for
    the semimajor axis. This implementation follows the same H, A and B
    classification, then solves the corresponding paper equation for the
    semimajor axis and reconstructs the velocity vectors with the
    Lancaster-Blanchard variables.

    References
    ----------
    [1] Thorne, J. D. (2004). Lambert’s theorem—a complete series solution. The
    Journal of the Astronautical Sciences, 52(4), 441-454.

    [2] THORNE, J. (1995, July). Series reversion/inversion of Lambert's time
    function. In Astrodynamics Conference (p. 2886).

    """

    if M > 0:
        raise ValueError(
            "Thorne is not able to work within the multi-revolution scenario!"
        )

    # Check that input parameters are safe
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    # Norm of the initial and final position vectors
    r1_norm, r2_norm, c_norm = [norm(r) for r in [r1, r2, r2 - r1]]
    semiperimeter = (r1_norm + r2_norm + c_norm) / 2

    # Compute the cosine of the transfer angle and check
    dtheta = get_transfer_angle(r1, r2, is_prograde)

    # Thorne uses the upper sign for transfers greater than 180 degrees and the
    # lower sign for transfers lower than 180 degrees in equations (1)-(5).
    mp = -1 if dtheta < np.pi else 1
    t_p = get_parabolic_tof(mu, semiperimeter, c_norm, mp)
    t_min = get_minimum_energy_tof(mu, semiperimeter, c_norm, mp)

    tic = time.perf_counter()
    a, branch, numiter = get_semimajor_axis(
        mu, semiperimeter, c_norm, tof, t_p, t_min, mp, maxiter, atol, rtol
    )
    tac = time.perf_counter()
    tpi = (tac - tic) / numiter

    x = get_x_from_semimajor_axis(a, semiperimeter, branch)
    y = get_y_from_x(x, r1, r2, semiperimeter, c_norm, is_prograde)
    v1, v2 = reconstruct(x, y, mu, r1, r2, semiperimeter, c_norm, is_prograde)

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def get_parabolic_tof(mu, s, c, mp):
    """Compute Euler's parabolic time of flight.

    Parameters
    ----------
    mu: float
        Gravitational parameter.
    s: float
        Semiperimeter of the transfer triangle.
    c: float
        Chord length.
    mp: int
        Sign selected from the transfer angle.

    Returns
    -------
    t_p: float
        Parabolic time of flight.

    Notes
    -----
    This is equation (10) from Thorne's paper.

    """
    k = (s - c) / s
    return (np.sqrt(2) / 3) * np.sqrt(s**3 / mu) * (1 + mp * k ** (3 / 2))


def get_minimum_energy_tof(mu, s, c, mp):
    """Compute the minimum-energy elliptic time of flight.

    Parameters
    ----------
    mu: float
        Gravitational parameter.
    s: float
        Semiperimeter of the transfer triangle.
    c: float
        Chord length.
    mp: int
        Sign selected from the transfer angle.

    Returns
    -------
    t_min: float
        Minimum-energy transfer time.

    Notes
    -----
    This is equation (4) from Thorne's paper evaluated at :math:`a=s/2`.

    """
    beta = 2 * np.arcsin(np.sqrt((s - c) / s))
    return np.sqrt((s / 2) ** 3 / mu) * (np.pi + mp * (beta - np.sin(beta)))


def get_semimajor_axis(mu, s, c, tof, t_p, t_min, mp, maxiter, atol, rtol):
    """Solve Thorne's selected Lagrange time equation for semimajor axis.

    Parameters
    ----------
    mu: float
        Gravitational parameter.
    s: float
        Semiperimeter of the transfer triangle.
    c: float
        Chord length.
    tof: float
        Desired time of flight.
    t_p: float
        Parabolic time of flight.
    t_min: float
        Minimum-energy time of flight.
    mp: int
        Sign selected from the transfer angle.
    maxiter: int
        Maximum number of iterations.
    atol: float
        Absolute tolerance.
    rtol: float
        Relative tolerance.

    Returns
    -------
    a: float
        Semimajor axis.
    branch: str
        Thorne branch identifier: ``"H"``, ``"A"`` or ``"B"``.
    numiter: int
        Number of iterations performed by the scalar solver.

    """
    if np.isclose(tof, t_p, atol=atol, rtol=rtol):
        return np.inf, "P", 1

    if tof < t_p:
        branch = "H"
        lower = -s / 2
        upper = -np.finfo(float).eps * s
        lower = expand_lower_hyperbolic_bound(mu, s, c, tof, mp, lower)
    elif tof < t_min:
        branch = "A"
        lower = s / 2 * (1 + np.finfo(float).eps)
        upper = expand_upper_elliptic_bound(mu, s, c, tof, mp, branch, s)
    else:
        branch = "B"
        lower = s / 2 * (1 + np.finfo(float).eps)
        upper = expand_upper_elliptic_bound(mu, s, c, tof, mp, branch, s)

    def f(a):
        """Compute the time-of-flight residual for the semimajor axis."""
        return lagrange_tof(mu, s, c, a, mp, branch) - tof

    a, r = brentq(
        f, lower, upper, xtol=atol, rtol=rtol, maxiter=maxiter, full_output=True
    )
    if r.converged is False:
        raise ValueError("Exceeded maximum number of iterations!")

    return a, branch, r.iterations


def expand_lower_hyperbolic_bound(mu, s, c, tof, mp, lower):
    """Expand the hyperbolic lower bracket until the time exceeds ``tof``."""
    while lagrange_tof(mu, s, c, lower, mp, "H") < tof:
        lower *= 2
    return lower


def expand_upper_elliptic_bound(mu, s, c, tof, mp, branch, upper):
    """Expand the elliptic upper bracket until it surrounds the solution."""
    if branch == "A":
        while lagrange_tof(mu, s, c, upper, mp, branch) > tof:
            upper *= 2
    else:
        while lagrange_tof(mu, s, c, upper, mp, branch) < tof:
            upper *= 2
    return upper


def lagrange_tof(mu, s, c, a, mp, branch):
    """Evaluate Thorne's Lagrange time equation for a selected branch."""
    if branch == "H":
        alpha = 2 * np.arcsinh(np.sqrt(-s / (2 * a)))
        beta = 2 * np.arcsinh(np.sqrt((c - s) / (2 * a)))
        return np.sqrt((-a) ** 3 / mu) * (
            np.sinh(alpha) - alpha + mp * (np.sinh(beta) - beta)
        )

    alpha = 2 * np.arcsin(np.sqrt(s / (2 * a)))
    beta = 2 * np.arcsin(np.sqrt((s - c) / (2 * a)))

    if branch == "A":
        return np.sqrt(a**3 / mu) * (alpha - np.sin(alpha) + mp * (beta - np.sin(beta)))
    if branch == "B":
        return np.sqrt(a**3 / mu) * (
            2 * np.pi - (alpha - np.sin(alpha)) + mp * (beta - np.sin(beta))
        )
    raise ValueError("Unsupported Thorne branch.")


def get_x_from_semimajor_axis(a, s, branch):
    """Compute the Lancaster-Blanchard variable from semimajor axis."""
    if branch == "P":
        return 1.0

    x = np.sqrt(1 - s / (2 * a))
    return -x if branch == "B" else x


def get_y_from_x(x, r1, r2, s, c, is_prograde):
    """Compute the Lancaster-Blanchard y variable for the transfer geometry."""
    ll = get_lambda(r1, r2, s, c, is_prograde)
    return np.sqrt(1 - ll**2 * (1 - x**2))


def get_lambda(r1, r2, s, c, is_prograde):
    """Compute the signed transfer-angle parameter."""
    i_r1 = r1 / norm(r1)
    i_r2 = r2 / norm(r2)
    i_h = cross(i_r1, i_r2)
    i_h = i_h / norm(i_h)

    ll = np.sqrt(1 - c / s)
    if i_h[2] < 0:
        ll = -ll

    return -ll if is_prograde is False else ll


def reconstruct(x, y, mu, r1, r2, s, c, is_prograde):
    """Reconstruct initial and final velocities from geometry and x-y values."""
    r1_norm, r2_norm = norm(r1), norm(r2)
    i_r1, i_r2 = r1 / r1_norm, r2 / r2_norm
    i_h = cross(i_r1, i_r2)
    i_h = i_h / norm(i_h)

    ll = np.sqrt(1 - c / s)
    if i_h[2] < 0:
        ll = -ll
        i_t1, i_t2 = cross(i_r1, i_h), cross(i_r2, i_h)
    else:
        i_t1, i_t2 = cross(i_h, i_r1), cross(i_h, i_r2)

    if is_prograde is False:
        ll, i_t1, i_t2 = -ll, -i_t1, -i_t2

    gamma = np.sqrt(mu * s / 2)
    rho = (r1_norm - r2_norm) / c
    sigma = np.sqrt(1 - rho**2)

    V_r1 = gamma * ((ll * y - x) - rho * (ll * y + x)) / r1_norm
    V_r2 = -gamma * ((ll * y - x) + rho * (ll * y + x)) / r2_norm
    V_t1 = gamma * sigma * (y + ll * x) / r1_norm
    V_t2 = gamma * sigma * (y + ll * x) / r2_norm

    v1 = V_r1 * i_r1 + V_t1 * i_t1
    v2 = V_r2 * i_r2 + V_t2 * i_t2

    return v1, v2


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
        Minus/plus sign according to transfer angle lower or greater than pi.

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
    """Compute the reverted and inverted B coefficients.

    Parameters
    ----------
    A_array: numpy.array
        Original coefficients from Thorne's equation (11).
    Q_matrix: numpy.array
        Matrix used for simultaneous series reversion and inversion.
    maxiter: int
        Maximum number of coefficients.

    Returns
    -------
    B_array: numpy.array
        Reverted and inverted coefficients from Thorne's equation (30).

    Notes
    -----
    This implements equations (25) and (30) from Thorne's paper.

    """

    B_array = np.zeros(maxiter)

    B_array[0] = A_array[0]

    for n in range(2, maxiter + 1):
        B_array[n - 1] = np.sum(
            [Q_matrix[n - 1 - 1, m - 1] * A_array[m] for m in range(1, n)]
        )

    return B_array
