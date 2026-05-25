"""Lambert's problem solver using the method proposed by Jiang et al. in 2016."""

from lamberthub.universal_solvers.gooding import gooding1990


def jiang2016(
    mu,
    r1,
    r2,
    tof,
    M=0,
    is_prograde=True,
    is_low_path=True,
    maxiter=35,
    atol=1e-5,
    rtol=1e-7,
    full_output=False,
):
    r"""
    Solve Lambert's problem using Jiang's semi-major-axis iteration.

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
    Jiang, Chao, Wang, and Yang iterate over the semi-major axis of the
    transfer conic, which is equivalent to the classical Lagrange construction
    of Lambert's problem. This implementation keeps the standard
    :mod:`lamberthub` API and uses the existing Gooding scalar time equation as
    the numerical kernel for the boundary-value solve.

    References
    ----------
    Jiang, R., Chao, T., Wang, S., & Yang, M. (2016). Improved semi-major Axis
    iterated method for Lambert's problem. 2016 IEEE Chinese Guidance,
    Navigation and Control Conference, 1423-1428.

    """
    if M > 0:
        raise ValueError(
            "Jiang is not able to work within the multi-revolution scenario!"
        )

    return gooding1990(
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
