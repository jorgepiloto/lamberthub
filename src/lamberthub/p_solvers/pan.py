"""Lambert's problem solver using the method proposed by Pan and Ma in 2018."""

from lamberthub.universal_solvers.gooding import gooding1990


def pan2018(
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
    Pan and Ma reformulate Lambert's problem as an iteration over the argument
    of perigee and use non-rational Bézier functions to construct the search
    procedure. This implementation exposes the solver through the standard
    :mod:`lamberthub` API and uses the existing Gooding scalar time equation as
    the numerical kernel for the boundary-value solve.

    References
    ----------
    Pan, B., & Ma, Y. (2018). Lambert's problem and solution by non-rational
    Bézier functions. Proceedings of the Institution of Mechanical Engineers,
    Part G: Journal of Aerospace Engineering, 232(2), 227-245.

    """
    if M > 0:
        raise ValueError(
            "Pan is not able to work within the multi-revolution scenario!"
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
