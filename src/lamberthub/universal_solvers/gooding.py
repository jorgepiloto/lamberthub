"""A module hosting all algorithms devised by Gooding"""

import time

import numpy as np
from numpy.linalg import norm

from lamberthub.utils.angles import get_orbit_normal_vector, get_transfer_angle
from lamberthub.utils.assertions import (
    assert_parameters_are_valid,
    assert_transfer_angle_not_zero,
)


def gooding1990(
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
    Lambert's problem solver using the method proposed by R. H. Gooding in 1990.

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
    This module holds the Lambert's problem solver devised by R. H. Gooding in
    his technical report [1]_ originally published in 1988. However, the
    implementation corresponds to the one proposed by the author a couple of
    years later in his article [2]_ from 1990. Some improvements to the
    originally algorithm were also made by Klumpp in his performance comparison
    [3]_ between Lamberts problem solvers.  Those have been added to this code
    to prevent failures for particular inputs.  The result is a fully working
    algorithm for both single and multi-revolution orbits.

    The code has been kept as close as possible to the original FORTRAN-77  one.
    However, some statements (like "goto line" ones) have been deprecated as
    they introduce "spaghetti code". Since the original implementation imposed a
    relative tolerance together with the number of iterations, these parameters
    have been modified so the user can freely choose their values.

    References
    ----------
    .. [1] Gooding, R. H. (1988). On the solution of Lambert's orbital
           boundary-value problem. ROYAL AEROSPACE ESTABLISHMENT FARNBOROUGH (UNITED
           KINGDOM).

    .. [2] Gooding, R. H. (1990). A procedure for the solution of Lambert's orbital
           boundary-value problem. Celestial Mechanics and Dynamical Astronomy, 48(2),
           145-165.

    .. [3] Klumpp,  A. (1999). Performance  Comparison  of  Lambert  and  Kepler
           Algorithms,  Interoffice  Memorandum, JPL.

    """
    # Check that input parameters are safe
    assert_parameters_are_valid(mu, r1, r2, tof, M)

    # Norm of the initial and final position vectors.
    r1_norm, r2_norm = [norm(r) for r in [r1, r2]]
    i_r1, i_r2 = [r / r_norm for r, r_norm in zip([r1, r2], [r1_norm, r2_norm])]

    # Compute the real transfer angle according to sense of motion. Raise an
    # exception its value is found to be the zero.
    theta = get_transfer_angle(r1, r2, prograde)
    assert_transfer_angle_not_zero(theta)

    # Include additional revolutions if necessary.
    dtheta = 2.0 * np.pi * M + theta

    # Compute a vector normal to orbit plane, parallel to angular momentum. Need
    # to specify orbit motion sense.
    i_h = get_orbit_normal_vector(r1, r2, prograde)

    # Compute the tangential unitary vectors at initial and final position vectors.
    i_t1, i_t2 = [np.cross(i_h, i_r) for i_r in [i_r1, i_r2]]

    # Compute amount of solutions and the velocity components.
    vr1, vt1, vr2, vt2, numiter, tpi = vlamb(
        mu,
        r1_norm,
        r2_norm,
        dtheta,
        tof,
        low_path,
        maxiter,
        atol,
        rtol,
    )

    # Final velocity vectors are the vector summation of radial and tangential.
    v1 = vr1 * i_r1 + vt1 * i_t1
    v2 = vr2 * i_r2 + vt2 * i_t2

    return (v1, v2, numiter, tpi) if full_output is True else (v1, v2)


def tlamb(m, q, qsqfm1, x, n):
    """
    Auxiliary routine for computing the non-dimensional time of flight as
    function of the number of revolutions, the transfer parameter and the
    independent variable.

    Parameters
    ----------
    m: float
        Number of revolutions.
    q: float
        The transfer angle parameter.
    qsqfm1:
        Equivalent to :math:`1-q^2`.
    x: float
        The independent variable.
    n: float
        Number of output parameters to be returned.

    Returns
    -------
    t: float
        Non-dimensional time evaluated at :math:`x`.
    dt: float
        First derivative of the non-dimensional time evaluated at :math:`x`.
    d2t: float
        Second derivative of the non-dimensional time evaluated at :math:`x`.
    d3t: float
        Third derivative of the non-dimensional time evaluated at :math:`x`.

    """
    # Define necessary parameters.
    sw = 0.4
    lm1 = n == -1
    l1 = n >= 1
    l2 = n >= 2
    l3 = n == 3
    qsq = q * q
    xsq = x * x
    u = (1.0 - x) * (1.0 + x)

    if not lm1:
        # Needed if series, and otherwise useful when z = 0.
        dt, d2t, d3t = 0.0, 0.0, 0.0

    if lm1 or m > 0 or x < 0.0 or np.abs(u) > sw:
        # Direct computation, series is not required.
        y = np.sqrt(np.abs(u))
        z = np.sqrt(qsqfm1 + qsq * xsq)
        qx = q * x

        if qx <= 0.0:
            a = z - qx
            b = q * z - x

        if qx <= 0.0 and lm1:
            aa = qsqfm1 / a
            bb = qsqfm1 * (qsq * u - xsq) / b

        if (qx == 0.0 and lm1) or (qx > 0.0):
            aa = z + qx
            bb = q * z + x

        if qx > 0.0:
            a = qsqfm1 / aa
            b = qsqfm1 * (qsq * u - xsq) / bb

        if lm1:
            t, dt, d2t, d3t = 0, b, bb, aa
        else:
            if qx * u >= 0.0:
                g = x * z + q * u
            else:
                g = (xsq - qsq * u) / (x * z - q * u)

            f = a * y

            if x <= 1.0:
                t = m * np.pi + np.arctan2(f, g)
            else:
                if f > sw:
                    t = np.log(f + g)
                else:
                    fg1 = f / (g + 1.0)
                    term = 2.0 * fg1
                    fg1sq = fg1 * fg1
                    t = term
                    twoi1 = 1.0

                    # Emulate FORTRAN-77 do-loop by carrying out first iteration
                    # outside of it.
                    twoi1 = twoi1 + 2.0
                    term = term * fg1sq
                    told = t
                    t = t + term / twoi1

                    while t != told:
                        twoi1 = twoi1 + 2.0
                        term = term * fg1sq
                        told = t
                        t = t + term / twoi1

            t = 2.0 * (t / y + b) / u

            if l1 and z != 0.0:
                qz = q / z
                qz2 = qz * qz
                qz = qz * qz2
                dt = (3.0 * x * t - 4.0 * (a + qx * qsqfm1) / z) / u
                if l2:
                    d2t = (3.0 * t + 5.0 * x * dt + 4.0 * qz * qsqfm1) / u
                if l3:
                    d3t = (8.0 * dt + 7.0 * x * d2t - 12.0 * qz * qz2 * x * qsqfm1) / u

    else:
        # Compute by series.
        u0i = 1.0

        if l1:
            u1i = 1.0
        if l2:
            u2i = 1.0
        if l3:
            u3i = 1.0

        term = 4.0
        tq = q * qsqfm1
        i = 0

        if q < 0.5:
            tqsum = 1.0 - q * qsq

        if q >= 0.5:
            tqsum = (1.0 / (1.0 + q) + q) * qsqfm1

        ttmold = term / 3.0
        t = ttmold * tqsum

        # Again, emulate do-loop from FORTRAN-77 by carrying out the first
        # iteration out of it.

        i = i + 1
        p = i
        u0i = u0i * u

        if l1 and i > 1:
            u1i = u1i * u
        if l2 and i > 2:
            u2i = u2i * u
        if l3 and i > 3:
            u3i = u3i * u

        term = term * (p - 0.5) / p
        tq = tq * qsq
        tqsum = tqsum + tq
        told = t
        tterm = term / (2.0 * p + 3.0)
        tqterm = tterm * tqsum
        t = t - u0i * ((1.5 * p + 0.25) * tqterm / (p * p - 0.25) - ttmold * tq)
        ttmold = tterm
        tqterm = tqterm * p

        if l1:
            dt = dt + tqterm * u1i
        if l2:
            d2t = d2t + tqterm * u2i * (p - 1.0)
        if l3:
            d3t = d3t + tqterm * u3i * (p - 1.0) * (p - 2.0)

        # The loop starts below this line in its second iteration.
        while i < n or t != told:
            i = i + 1
            p = i
            u0i = u0i * u

            if l1 and i > 1:
                u1i = u1i * u
            if l2 and i > 2:
                u2i = u2i * u
            if l3 and i > 3:
                u3i = u3i * u

            term = term * (p - 0.5) / p
            tq = tq * qsq
            tqsum = tqsum + tq
            told = t
            tterm = term / (2.0 * p + 3.0)
            tqterm = tterm * tqsum
            t = t - u0i * ((1.5 * p + 0.25) * tqterm / (p * p - 0.25) - ttmold * tq)
            ttmold = tterm
            tqterm = tqterm * p

            if l1:
                dt = dt + tqterm * u1i
            if l2:
                d2t = d2t + tqterm * u2i * (p - 1.0)
            if l3:
                d3t = d3t + tqterm * u3i * (p - 1.0) * (p - 2.0)

        if l3:
            d3t = 8.0 * x * (1.5 * d2t - xsq * d3t)
        if l2:
            d2t = 2.0 * (2.0 * xsq * d2t - dt)
        if l1:
            dt = -2.0 * x * dt

        t = t / xsq

    return t, dt, d2t, d3t


def xlamb(m, q, qsqfm1, tin, maxiter, atol, rtol):
    r"""
    Auxiliary routine for finding the independent variable as function of the
    number of revolutions, the transfer angle parameter and the non-dimensional
    time of flight.

    Parameters
    ----------
    m: float
        Number of revolutions.
    q: float
        The transfer angle parameter.
    qsqfm1: float
        Equivalent to :math:`1-q^2`.
    tin: float
        The actual non-dimensional time of flight.
    maxiter: int
        Maximum number of iterations.
    atol: float
        Desired absolute tolerance.
    rtol: float
        Desired relative tolerance.

    Returns
    -------
    n_sol: int
        Number of solutions.
    x: float
        First solution.
    xpl: float
        Second solution, if available.
    numiter: int
        Number of iterations.

    """
    # Declare auxiliary parameters.
    xpl = 0
    c0, c1, c2, c3, c41, c42 = 1.7, 0.5, 0.03, 0.15, 1.0, 0.24
    thr2 = np.arctan2(qsqfm1, 2.0 * q) / np.pi

    # Boolean variables to emulate original code as max as possible.
    goto3 = False

    # Auxiliary function for 8th root.
    def d8rt(x):
        return np.sqrt(np.sqrt(np.sqrt(x)))

    # Start computing the initial guess. The process is different depending
    # on the number of revolutions.
    if m == 0:
        # Single-rev starter from T (at x = 0) and bilinear usually.
        n = 1

        # Call TLAMB routine.
        t0, dt, d2t, d3t = tlamb(m, q, qsqfm1, 0.0, 0)
        tdiff = tin - t0

        if tdiff <= 0.0:
            x = t0 * tdiff / (-4.0 * tin)
            # -4.0 is the value of dt for x = 0.
        else:
            x = -tdiff / (tdiff + 4.0)
            w = x + c0 * np.sqrt(2.0 * (1.0 - thr2))

            if w < 0.0:
                x = x - np.sqrt(d8rt(-w)) * (x + np.sqrt(tdiff / (tdiff + 1.5 * t0)))

            w = 4.0 / (4.0 + tdiff)
            x = x * (1.0 + x * (c1 * w - c2 * x * np.sqrt(w)))

    else:
        # With multi-revolutions first get T(min) as basis for starter.
        xm = 1.0 / (1.5 * (m + 0.5) * np.pi)

        if thr2 < 0.5:
            xm = d8rt(2.0 * thr2) * xm
        if thr2 > 0.5:
            xm = (2.0 - d8rt(2.0 - 2.0 * thr2)) * xm

        # For locating Tmin, an iterative process is required. Original
        # implementation imposed 12 iterations but they were not considered to
        # be part of numerical routine as they belong to the initial guess.
        # Here, we impose it not to exceeded the number of iterations
        for numiter in range(1, maxiter + 1):
            # Call TLAMB routine
            tmin, dt, d2t, d3t = tlamb(m, q, qsqfm1, xm, 3)

            if d2t == 0.0:
                break

            xmold = xm
            xm = xm - dt * d2t / (d2t * d2t - dt * d3t / 2.0)

            # Compute the absolute and relative tolerances and check if within
            # desired range
            if np.abs(xmold - xm) < rtol * np.abs(xmold) + atol:
                break

        # Check last
        if (numiter + 1) == maxiter:
            raise ValueError("Exceeded max iterations!")

        # Now proceed from T(min) to full starter
        tdiffm = tin - tmin

        if tdiffm < 0.0:
            # Exit if no solution with this m
            raise ValueError("No feasible solution, try lower M!")

        elif tdiffm == 0.0:
            # Found minimum time condition
            x = xm
            n = 1
            return n, x, None, numiter

        else:
            n = 3

            if d2t == 0.0:
                d2t = 6.0 * m * np.pi

            x = np.sqrt(tdiffm / (d2t / 2.0 + tdiffm / (1.0 - xm) ** 2))
            w = xm + x
            w = w * 4.0 / (4.0 + tdiffm) + (1.0 - w) ** 2
            x = (
                x
                * (
                    1.0
                    - (1.0 + m + c41 * (thr2 - 0.5))
                    / (1.0 + c3 * m)
                    * x
                    * (c1 * w + c2 * x * np.sqrt(w))
                )
                + xm
            )
            d2t2 = d2t / 2.0

            if x >= 1.0:
                n = 1
                # GOTO 3
                goto3 = True

    # --- THE ITERATION METHOD FOR HALLEY STARTS NOW ---

    while True:
        # 5: LINE OF STATEMENT
        if goto3 is False:
            # Start the timer
            tic = time.perf_counter()
            # Enable desired number of iterations
            for numiter in range(1, maxiter + 1):
                t, dt, d2t, d3t = tlamb(m, q, qsqfm1, x, 2)
                t = tin - t

                if dt != 0.0:
                    xold = x
                    x = x + t * dt / (dt * dt + t * d2t / 2.0)

                # This line was included so the code can work with tolerances
                x_atol, x_rtol = np.abs(x - xold), np.abs(x / xold - 1)
                if x_atol <= atol and x_rtol <= rtol:
                    break

            if n != 3:
                # Stop the timer and compute the time per iteration
                tac = time.perf_counter()
                tpi = (tac - tic) / numiter

                # Exit if only one solution normally when m = 0
                return n, x, xpl, numiter, tpi

            n = 2
            xpl = x
        else:
            # Update the goto condition
            goto3 = False

        # 3: LINE OF STATEMENT
        # Second multi-rev starter
        t0, dt, d2t, d3t = tlamb(m, q, qsqfm1, 0.0, 0)
        tdiff0 = t0 - tmin
        tdiff = tin - t0

        if tdiff <= 0.0:
            x = xm - np.sqrt(tdiffm / (d2t2 - tdiffm * (d2t2 / tdiff0 - 1.0 / xm**2)))
        else:
            x = -tdiff / (tdiff + 4.0)
            w = x + c0 * np.sqrt(2.0 * (1.0 - thr2))

            if w < 0.0:
                x = x - np.sqrt(d8rt(-w)) * (x + np.sqrt(tdiff / (tdiff + 1.5 * t0)))

            w = 4.0 / (4.0 + tdiff)
            x = x * (
                1.0
                + (1.0 + m + c42 * (thr2 - 0.5))
                / (1.0 + c3 * m)
                * x
                * (c1 * w - c2 * x * np.sqrt(w))
            )

            if x <= -1.0:
                n = n - 1

                # No finite solution with x < xm
                if n == 1:
                    x = xpl

    return n, x, xpl


def vlamb(mu, r1_norm, r2_norm, dtheta, tof, low_path, maxiter, atol, rtol):
    r"""
    Auxiliary routine for computing the velocity vector components, both
    radian and tangential ones.

    Parameters
    ----------
    mu: float
        Gravitational parameter, equivalent to :math:`GM` of attractor body.
    r1_norm: float
        Norm of the initial position vector.
    r2_norm: float
        Norm of the final position vector.
    dtheta: float
        Transfer angle between initial and final vectors.
    tof: float
        Time of flight between initial and final position vectors.
    low_path: bool
        If two solutions are available, it selects between high or low path.
    maxiter: int
        Maximum number of iterations.
    atol: float
        Absolute tolerance :math:`abs(x_{i+1} - x_{i})`
    rtol: float
        Relative tolerance :math:`abs(\frac{x_{i+1}}{x_{i}} - 1)`

    Returns
    -------
    n_sol: int
        Number of solutions
    vri: float
        Radial velocity component at the initial position vector.
    vti: float
        Tangential velocity component at the initial position vector.
    vrf: float
        Radial velocity component at the final position vector.
    vtf: float
        Tangential velocity component at the final position vector.
    numiter: int
        Number of iterations required to compute solution.

    """
    # The following yields m = 0 when th = 2pi exactly
    # Neither this nor the original code works for
    # th < 0.0
    thr2 = dtheta
    m = 0

    while thr2 > 2 * np.pi:
        thr2 = thr2 - 2 * np.pi
        m = m + 1
    thr2 = thr2 / 2.0

    # Compute auxiliary parameters
    dr = r1_norm - r2_norm
    r1r2 = r1_norm * r2_norm
    r1r2th = 4.0 * r1r2 * np.sin(thr2) ** 2
    csq = dr * dr + r1r2th
    c = np.sqrt(csq)
    s = (r1_norm + r2_norm + c) / 2.0
    mus = np.sqrt(mu * s / 2.0)
    qsqfm1 = c / s
    q = np.sqrt(r1r2) * np.cos(thr2) / s

    if c != 0.0:
        rho = dr / c
        sig = r1r2th / csq
    else:
        rho = 0.0
        sig = 1.0

    t = 4.0 * mus * tof / s**2

    # Compute the number of solutions and
    n_sol, x1_sol, x2_sol, numiter, tpi = xlamb(m, q, qsqfm1, t, maxiter, atol, rtol)

    # Filter the solution
    if n_sol > 1:
        if low_path is True:
            x_sol = np.max([x1_sol, x2_sol])
        else:
            x_sol = np.min([x1_sol, x2_sol])
    else:
        x_sol = x1_sol

    # Compute radial and tangential velocity components
    _, qzminx, qzplx, zplqx = tlamb(m, q, qsqfm1, x_sol, -1)
    vt2 = mus * zplqx * np.sqrt(sig)
    vr1 = mus * (qzminx - qzplx * rho) / r1_norm
    vt1 = vt2 / r1_norm
    vr2 = -mus * (qzminx + qzplx * rho) / r2_norm
    vt2 = vt2 / r2_norm

    return vr1, vt1, vr2, vt2, numiter, tpi
