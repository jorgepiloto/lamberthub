"""Module containing Kepler's equations.

Equations are declared for each one of the particular orbit shapes, that is
elliptical, parabolic and hyperbolic.

The required formulas are found to be within Vallado's [1]_ manual.

References
----------
.. [1] Vallado, D. A. (2001). Fundamentals of astrodynamics and applications
   (Vol. 12). Springer Science & Business Media.

"""

import numpy as np

from lamberthub.utils.angles import nu_to_B, nu_to_E, nu_to_H


def kepler_elliptic(E, ecc):
    """Compute the time since perigee passage at given eccentric anomaly.

    Parameters
    ----------
    E: float
        Eccentric anomaly.
    ecc: float
        Eccentricity of the orbit. Must be between (0,1).

    Returns
    -------
    M: float
        Time since perigee passage.

    """
    M = E - ecc * np.sin(E)
    return M


def kepler_parabolic(B):
    """Compute the time since perigee passage at given parabolic anomaly.

    Parameters
    ----------
    B: float
        Parabolic anomaly.

    Returns
    -------
    Mp: float
        Parabolic mean motion

    """
    Mp = B + (1 / 3) * B**3
    return Mp


def kepler_hyperbolic(H, ecc):
    """Compute the time since perigee passage at given hyperbolic anomaly.

    Parameters
    ----------
    H: float
        Hyperbolic anomaly.
    ecc:
        Eccentricity of the orbit.

    Returns
    -------
    Mh: float
        Hyperbolic mean motion

    """
    Mh = ecc * np.sinh(H) - H
    return Mh


def kepler_from_nu(nu: float, ecc: float) -> float:
    """Compute the mean anomaly at a given true anomaly.

    Parameters
    ----------
    nu: float
        True anomaly.
    ecc: float
        Orbit's eccentricity.

    Returns
    -------
    M: float
        Mean anomaly.

    """
    if ecc < 0:
        raise ValueError("Eccentricity cannot be negative!")
    elif 0 <= ecc < 1:
        E = nu_to_E(nu, ecc)
        M = kepler_elliptic(E, ecc)
    elif ecc == 1:
        B = nu_to_B(nu)
        M = kepler_parabolic(B)
    else:
        H = nu_to_H(nu, ecc)
        M = kepler_hyperbolic(H, ecc)

    return M
