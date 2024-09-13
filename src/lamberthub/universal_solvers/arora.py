"""A module hosting all algorithms devised by Arora"""

import time

import numpy as np
from numpy.linalg import norm

from lamberthub.utils.angles import get_transfer_angle
from lamberthub.utils.assertions import (
    assert_parameters_are_valid,
    assert_transfer_angle_not_zero,
)


def arora2013(
    mu,
    r1,
    r2,
    tof,
    M=0,
    prograde=True,
    low_path=True,
    maxiter=35,
    atol=1e-5,
    rtol=1e-7,
    full_output=False,
):
    r"""
    Solves Lambert problem using Arora's devised algorithm

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
        Absolute tolerance :math:`abs(x_{i+1} - x_{i})`
    rtol: float
        Relative tolerance :math:`abs(\frac{x_{i+1}}{x_{i}} - 1)`
    full_output: bool
        If True, the number of iterations and time per iteration are also returned.

    Returns
    -------
    v1: numpy.array
        Initial velocity vector
    v2: numpy.array
        Final velocity vector
    numiter: int
        Number of iterations.
    tpi: float
        Time per iteration in seconds.

    Notes
    -----
    Lambert's problem solver using the method proposed by Nitin Arora and Ryan
    P. Rusell in 2013, see [1]. This algorithm exploits the universal formulae
    by defining a new cosine-based transformation and developing a robust
    initial guess. Although based on arbitrary conditions, the algorithm shows a
    high performance.

    References
    ----------
    [1] Arora, N., & Russell, R. P. (2013). A fast and robust multiple
    revolution Lambert algorithm using a cosine transformation. Paper AAS, 13,
    728.

    """
    # TODO: implement solver for the multi-revolution case
    # https://github.com/jorgepiloto/lamberthub/issues/3
    if M > 0:
        raise NotImplementedError(
            "See https://github.com/jorgepiloto/lamberthub/issues/3"
        )

    # Check that input parameters are safe
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    # Norm of the initial and final position vectors
    r1_norm, r2_norm, c_norm = [norm(r) for r in [r1, r2, (r2 - r1)]]

    # Unitary vectors along the radial directions
    i_r1, i_r2 = [r / r_norm for r, r_norm in zip([r1, r2], [r1_norm, r2_norm])]

    # Compute the transfer angle and the transfer angle parameter
    dtheta = get_transfer_angle(r1, r2, prograde)
    assert_transfer_angle_not_zero(dtheta)
    d = 1 if dtheta <= np.pi else -1

    # Compute the characteristic length and time variables
    L_ref = r1_norm
    T_ref = np.sqrt(L_ref**3 / mu)

    # All variables named "hat" have been non-dimensionalized w.r.t. previously
    # computed reference quantities.
    r1_hat, r2_hat = [r_vec / L_ref for r_vec in [r1, r2]]
    r1hat_norm, r2hat_norm = [norm(r_vec) for r_vec in [r1_hat, r2_hat]]
    tof = tof / T_ref
    mu_hat = mu * (T_ref**2 / L_ref**3)

    # Auxiliary variables
    S = np.sqrt((r1hat_norm + r2hat_norm) ** 3 / mu_hat)

    # Solve Lambert's geometry parameter
    tau = (
        d
        * np.sqrt(r1hat_norm * r2hat_norm * (1 + np.cos(dtheta)))
        / (r1hat_norm + r2hat_norm)
    )

    # Compute the parabolic time of flight
    tof_p = S * np.sqrt(1 - np.sqrt(2) * tau) * (tau + np.sqrt(2)) / 3

    # Compare actual and parabolic TOF to generate the initial guess
    if tof <= tof_p:
        # The final orbit is expected to be hyperbolic. The initial guess is
        # determined by the formulae given by (TABLE 2) from official report.
        # Three regions, depending on the transfer angle are possible: one where
        # d = 1 and other two with d = -1. These last two are named H1 or H2,
        # where H refers to Hyperbolic and 1 or 2 to the type of region.

        # Equations (48a and 48b) from official report, used to fix the limits
        # between H1 and H2 regions
        tof20 = S * np.sqrt(1 - 20 * tau) * (tau + 0.04940968903 * (1 - 20 * tau))
        tof100 = S * np.sqrt(1 - 100 * tau) * (tau + 0.00999209404 * (1 - 100 * tau))

        # Apply the initial guess associated to a particular hyperbolic region.

        if d == 1:
            # Value from the first row of (Table 2)
            k_n, k_m = np.sqrt(2), 1 / tau
            k_i = (k_n + k_m) / 2
            Z = 1.0 / np.sqrt(2)
            alpha = 1 / 2
            F_0 = tof_p
            F_1 = 0.0
            W = _get_W(k_i, M)
            F_i = _get_TOF(k_i, tau, S, W)
            F_star = tof

            # Compute the Sundman transformation and the initial guess
            x_star = _get_x(F_0, F_1, F_i, F_star, Z, alpha)
            k = k_n + (k_m - k_n) * x_star

        elif (d == -1) and (tof > tof20):
            # Initial guess falls in hyperbolic region H1
            k_n, k_m = np.sqrt(2), 20
            k_i = (2 * k_n + k_m) / 3
            Z = 1 / 3
            alpha = 1
            F_0 = tof_p
            F_1 = tof20
            W = _get_W(k_i, M)
            F_i = _get_TOF(k_i, tau, S, W)
            F_star = tof

            # Compute the Sundman transformation and the initial guess
            x_star = _get_x(F_0, F_1, F_i, F_star, Z, alpha)
            k = k_n + (k_m - k_n) * x_star

        elif (d == -1) and (tof <= tof20):
            # Initial guess falls in the H2 region. In this case, equation (47)
            # is directly applied.
            t_star, t_0, t_1 = tof, tof20, tof100
            k = (
                (t_1 * (t_0 - t_star) * 10 - t_0 * np.sqrt(20) * (t_1 - t_star))
                / (t_star * (t_0 - t_1))
            ) ** 2

    else:
        # The final orbit is expected to be and ellipse. However, these orbits
        # can also be found for multi-revolution problems, so it is required to
        # distinguish between them.

        if M == 0:
            # For the zero-revolution initial guess, there exist a total of four
            # regions named E1, E2, E3 and E4, limited by the following

            # A collection of auxiliary independent variable values
            k_set = -1.41, -1.38, -1.00, -1 / 2, 0, 1 / np.sqrt(2)

            # Precomputed values of W(k) required for the zero-rev initial guess.
            # Can be obtained via _get_W(k, 0)
            W_set = (
                4839.684497246,  # k = -1.41
                212.087279879,  # k = -1.38
                5.712388981,  # k = -1.00
                1.954946607,  # k = -0.50
                1.110720735,  # k = 0.00
                0.6686397730,  # k = 1 / sqrt(2)
            )

            # Time of flight for each one of the previous auxiliary values.
            # These values could also be precomputed but they are not provided
            # in the original report.
            t_set = [_get_TOF(k, tau, S, W) for k, W in zip(k_set, W_set)]
            tof_m141, tof_m138, tof_m1, tof_m1half, tof_0, tof_1oversq2 = t_set

            # Filter out the region. Associated values can be found within
            # (TABLE 3) and (TABLE 4) from official report.

            if tof <= tof_0:
                # Region E1 applies
                k_n, k_m, k_i = 0, np.sqrt(2), 1 / np.sqrt(2)
                Z, alpha = 1 / 2, 1
                F_0, F_1, F_i, F_star = tof_0, tof_p, tof_1oversq2, tof
                x_star = _get_x(F_0, F_1, F_i, F_star, Z, alpha)
                k = k_n + (k_m - k_n) * x_star

            elif tof_0 <= tof <= tof_m1:
                # Region E2 applies
                k_n, k_m, k_i = 0, -1, -1 / 2
                Z, alpha = 1 / 2, 1
                F_0, F_1, F_i, F_star = tof_0, tof_m1, tof_m1half, tof
                x_star = _get_x(F_0, F_1, F_i, F_star, Z, alpha)
                k = k_n + (k_m - k_n) * x_star

            elif tof_m1 <= tof <= tof_m138:
                # Region E3 applies
                k_n, k_m, k_i = -1, -np.sqrt(2), -1.38
                c1, c2, c3, c4, alpha = 540649 / 3125, 256, 1, 1, 16
                F_n, F_i, F_star = tof_m1**-1, tof_m138**-1, tof**-1
                gamma1, gamma2, gamma3 = _get_gammas(F_i, F_n, F_star)
                k = -c4 * (
                    ((gamma1 * c1 - c3 * gamma3) * c2 + c3 * c1 * gamma2)
                    / (gamma3 * c1 - c3 * gamma1 - gamma2 * c2)
                ) ** (1 / alpha)

            else:
                # Region E4 applies, that is tof >= tof_m138
                k_n, k_m, k_i = -1.38, -np.sqrt(2), -1.41
                c1, c2, c3, c4, alpha = (
                    49267 / 27059,
                    67286 / 17897,
                    2813 / 287443,
                    4439 / 3156,
                    243,
                )
                F_n, F_i, F_star = tof_m138**-1, tof_m141**-1, tof**-1
                gamma1, gamma2, gamma3 = _get_gammas(F_i, F_n, F_star)
                k = -c4 * (
                    ((gamma1 * c1 - c3 * gamma3) * c2 + c3 * c1 * gamma2)
                    / (gamma3 * c1 - c3 * gamma1 - gamma2 * c2)
                ) ** (1 / alpha)

        else:
            # Multirevolutions apply. For this last case, Arora also divided the
            # solution space into four regions named M1, M2, M3 and M4, being th
            # letter M associated with Multi-revolutions. The procedure, as
            # explained in the official report, consists into two parts: compute
            # k_bi and using this value to get the final initial guess.
            raise NotImplementedError("Still need to implement Arora's multirev.")

    # Now that the initial guess has been performed, it is possible to start the
    # iterative process. Initialize the timer also.
    tic = time.perf_counter()
    for numiter in range(1, maxiter + 1):
        # Evaluate the auxiliary function, its first and second derivative
        # w.r.t. to the independent variable
        W = _get_W(k, M)
        Wp = _get_Wprime(k, W)
        Wpp = _get_W2prime(k, W, Wp)

        # Evaluate the predicted time of flight with the current value of k
        c = (1 - k * tau) / tau
        tofc = _get_TOF(k, tau, S, W)

        # Check computed time of flight matches target one
        if np.abs(tof - tofc) <= atol:
            tac = time.perf_counter()
            break

        # Compute the time derivatives to proceed by using Halley's method
        tofc_p = (-tofc / (2 * c)) + S * tau * np.sqrt(c * tau) * (Wp * c - W)
        tofc_pp = (-tofc / (4 * c**2)) + S * tau * np.sqrt(c * tau) * (
            W / c + c * Wpp - 3 * Wp
        )

        # Solve Halley's step and check if convergence was achieved
        deltak = -(tofc - tof) / (tofc_p - (tofc - tof) * tofc_pp / (2.0 * tofc_p))

        # Update the value of the independent variable and carry a new iteration
        k += deltak

        # Bound the value of k if required
        if k < -np.sqrt(2):
            k = -np.sqrt(2) + 1e-12

        if k < np.sqrt(2) and tof < tof_p:
            k = np.sqrt(2) + 1e-12

        if k > np.sqrt(2) and tof > tof_p:
            k = np.sqrt(2) - 1e-12

        if tof < tof_p and d > 0 and (1 - tau * k) < 0:
            k = 1 / tau - 1e-12

    # Compute the time per iteration
    tpi = (tac - tic) / numiter

    # Evaluate f and g functions. These are equations (32) and (33)
    f = 1 - (1 - k * tau) * (r1hat_norm + r2hat_norm) / r1hat_norm
    g = S * tau * np.sqrt((1 - k * tau) * mu_hat)
    g_dot = 1 - (1 - k * tau) * (r1hat_norm + r2hat_norm) / r2hat_norm

    # Finally, compute the initial and final velocity vectors using equations
    # (35) and (36)
    v1 = (r2_hat - f * r1_hat) / g * (L_ref / T_ref)
    v2 = (g_dot * r2_hat - r1_hat) / g * (L_ref / T_ref)

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def _get_gammas(F_i, F_n, F_star):
    """Compute different gamma values"""
    gamma1, gamma2, gamma3 = (
        F_i * (F_star - F_n),
        F_star * (F_n - F_i),
        F_n * (F_star - F_i),
    )
    return gamma1, gamma2, gamma3


def _get_x(F_0, F_1, F_i, F_star, Z, alpha):
    """Computes Sundman transformation variable.

    Parameters
    ----------
    F_0: float
        First boundary coefficient.
    F_1: float
        Second boundary coefficient.
    F_i: float
        Third boundary coefficient.
    F_star: float
        Last boundary coefficient.
    Z: float
        Auxiliary constant.
    alpha: float
        Auxiliary constant.

    Returns
    -------
    x: float
        Sundman transformation value.

    Notes
    -----
    This is equation (44) from original report.

    """
    x = (
        (Z * (F_0 - F_star) * (F_1 - F_i))
        / ((F_i - F_star) * (F_1 - F_0) * Z + (F_0 - F_i) * (F_1 - F_star))
    ) ** (1 / alpha)
    return x


def _get_W(k, M, epsilon=2e-2):
    """
    Evaluates the auxiliary function at particular value of the independent
    variable.

    Parameters
    ----------
    k: float
        Independent variable.
    M: int
        Number of revolutions
    epsilon: float
        Tolerance parameter. Default value as in the original report.

    Returns
    -------
    W: float
        Value of the auxiliary function.

    Notes
    -----
    This is equation (27) from official report.

    """
    # Evaluate the sign of k
    m = 2 - k**2
    sgn_k = np.sign(k)
    sq2 = np.sqrt(2)

    # Apply a particular formulae depending on the case

    if -sq2 <= k < (sq2 - epsilon):
        # Elliptical orbits
        W = ((1 - sgn_k) * np.pi + sgn_k * np.arccos(1 - m) + 2 * np.pi * M) / (
            np.sqrt(m**3)
        ) - k / m

    elif k > sq2 + epsilon:
        # Hyperbolic orbits
        W = -np.arccosh(1 - m) / np.sqrt(-(m**3)) - k / m

    elif sq2 - epsilon <= k <= sq2 + epsilon:
        # Direct transfer, no complete revolutions (M = 0)

        # Allocate auxiliary variables
        v = k - sq2
        v2, v3, v4, v5, v6, v7, v8 = [v**i for i in range(2, 9)]

        W = (
            (np.sqrt(2) / 3)
            - (1 / 5) * v
            + (2 / 35) * sq2 * v2
            - (2 / 63) * v3
            + (2 / 231) * sq2 * v4
            - (2 / 429) * v5
            + (8 / 6435) * sq2 * v6
            - (8 / 12155) * v7
            + (8 / 46189) * sq2 * v8
        )
    else:
        raise ValueError("Did not found a suitable equaiton to find W!")

    return W


def _get_Wsprime(k):
    """Evaluate the first derivative of Ws w.r.t. independent variable k.

    Parameters
    ----------
    k: float
        The independent variable.

    Returns
    -------
    Ws_prime: float
        Value of the first derivative w.r.t. to k.

    Notes
    -----
    This equation was not provided in the original report, probably because
    author assumed it was trivial.

    """
    # Allocate auxiliary variables
    sq2 = np.sqrt(2)
    v = k - sq2
    v2, v3, v4, v5, v6, v7 = [v**i for i in range(2, 8)]

    Ws_prime = (
        -1 / 5
        + sq2 * v * (4 / 35)
        - v2 * (6 / 63)
        + sq2 * v3 * (8 / 231)
        - v4 * (10 / 429)
        + sq2 * v5 * (48 / 6435)
        - v6 * (56 / 12155)
        + sq2 * v7 * (64 / 46189)
    )

    return Ws_prime


def _get_Ws2prime(k):
    """Evaluate the second derivative of Ws w.r.t. independent variable k.

    Parameters
    ----------
    k: float
        The independent variable.

    Returns
    -------
    Ws_2prime: float
        Value of the second derivative w.r.t. to k.

    Notes
    -----
    This equation was not provided in the original report, probably because
    author assumed it was trivial.

    """
    # Allocate auxiliary variables
    sq2 = np.sqrt(2)
    v = k - sq2
    v2, v3, v4, v5, v6 = [v**i for i in range(2, 7)]

    Ws_2prime = (
        sq2 * (4 / 35)
        - v * (12 / 63)
        + sq2 * v2 * (24 / 231)
        - v3 * (40 / 429)
        + sq2 * v4 * (240 / 6435)
        - v5 * (336 / 12155)
        + sq2 * v6 * (448 / 46189)
    )

    return Ws_2prime


def _get_Wprime(k, W, epsilon=2e-2):
    """
    Evaluates the first derivative of the auxiliary function w.r.t. the
    independent variable k.

    Parameters
    ----------
    k: float
        The independent variable.
    W: float
        The auxiliary function value.

    Returns
    -------
    W_prime: float
        The value of the first derivative of the auxiliary function.

    Notes
    -----
    This is equation set (38) from official report.

    """
    # Evaluate m
    m = 2 - k**2

    # Filter case
    if k < np.sqrt(2) - epsilon:
        W_prime = (-2 + 3 * W * k) / m

    elif np.sqrt(2) - epsilon < k < np.sqrt(2) + epsilon:
        # The partial derivative of Ws/k was not provided in the original
        # report. This is probably because it is trivial.
        W_prime = _get_Wsprime(k)

    else:
        # TODO: a minus sign before "m" is imposed in original report
        # https://github.com/jorgepiloto/lamberthub/issues/151

        # W_prime = (-2 + 3 * W * k) / -m
        W_prime = (-2 + 3 * W * k) / m

    return W_prime


def _get_W2prime(k, W, W_prime, epsilon=2e-2):
    """
    Evaluates the second derivative of the auxiliary function w.r.t. the
    independent variable k.

    Parameters
    ----------
    k: float
        The independent variable.
    W: float
        The auxiliary function value.

    Returns
    -------
    W_2prime: float
        The value of the second derivative of the auxiliary function.

    Notes
    -----
    This is equation set (39) from official report.

    """
    # Evaluate m
    m = 2 - k**2

    # Filter case
    if k < np.sqrt(2) - epsilon:
        W_2prime = (5 * W_prime * k + 3 * W) / m

    elif np.sqrt(2) - epsilon < k < np.sqrt(2) + epsilon:
        # The partial derivative of Ws/k was not provided in the original
        # report. This is probably because it is trivial.
        W_2prime = _get_Ws2prime(k)

    else:
        # TODO: a minus sign before "m" is imposed in original report
        # https://github.com/jorgepiloto/lamberthub/issues/151

        # W_2prime = (5 * W_prime * k + 3 * W) / -m
        W_2prime = (5 * W_prime * k + 3 * W) / m

    return W_2prime


def _get_TOF(k, tau, S, W):
    """Evaluates the time of flight at a particular value of the independent variable.

    Parameters
    ----------
    k: float
        The independent variable.
    tau: float
        Lambert's geometry parameter.
    S: float
        Auxiliary variable.

    Returns
    -------
    TOF: float
        Computed time of flight.

    Notes
    -----
    This is equation (26) form official report.

    """
    TOF = S * np.sqrt(1 - k * tau) * (tau + (1 - k * tau) * W)
    return TOF
