"""
Holds routines for converting between different orbital elements sets.

.. note::
    Some of these functions were directly taken from the poliastro[1] library,
    see the references section for more information about this powerful
    software.

    Copyright (c) 2012-2021 Juan Luis Cano Rodríguez and the poliastro
    development team.

References
----------
[1] poliastro: https://github.com/poliastro/poliastro/

"""

import numpy as np
from numpy import cos, cross, sin, sqrt
from numpy.linalg import norm

from lamberthub.utils.angles import E_to_nu, H_to_nu


def rv_pqw(k, p, ecc, nu):
    r"""Return r and v vectors in perifocal frame.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    p : float
        Semi-latus rectum or parameter (km).
    ecc : float
        Eccentricity.
    nu: float
        True anomaly (rad).

    Returns
    -------
    r: ndarray
        Position. Dimension 3 vector
    v: ndarray
        Velocity. Dimension 3 vector

    Notes
    -----
    These formulas can be checked at Curtis 3rd. Edition, page 110. Also the
    example proposed is 2.11 of Curtis 3rd Edition book.

    .. math::

        \vec{r} = \frac{h^2}{\mu}\frac{1}{1 + e\cos(\theta)}\begin{bmatrix}
        \cos(\theta)\\
        \sin(\theta)\\
        0
        \end{bmatrix} \\\\\\

        \vec{v} = \frac{h^2}{\mu}\begin{bmatrix}
        -\sin(\theta)\\
        e+\cos(\theta)\\
        0
        \end{bmatrix}

    """
    pqw = np.array([[cos(nu), sin(nu), 0], [-sin(nu), ecc + cos(nu), 0]]) * np.array(
        [[p / (1 + ecc * cos(nu))], [sqrt(k / p)]]
    )
    return pqw


def coe_rotation_matrix(inc, raan, argp):
    """Create a rotation matrix for coe transformation."""
    r = rotation_matrix(raan, 2)
    r = r @ rotation_matrix(inc, 0)
    r = r @ rotation_matrix(argp, 2)
    return r


def coe2rv(k, p, ecc, inc, raan, argp, nu):
    r"""Convert from classical orbital to state vectors.

    Classical orbital elements are converted into position and velocity
    vectors by `rv_pqw` algorithm. A rotation matrix is applied to position
    and velocity vectors to get them expressed in terms of an IJK basis.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2).
    p : float
        Semi-latus rectum or parameter (km).
    ecc : float
        Eccentricity.
    inc : float
        Inclination (rad).
    omega : float
        Longitude of ascending node (rad).
    argp : float
        Argument of perigee (rad).
    nu : float
        True anomaly (rad).

    Returns
    -------
    r_ijk: np.array
        Position vector in basis ijk.
    v_ijk: np.array
        Velocity vector in basis ijk.

    Notes
    -----

    .. math::
        \begin{align}
            \vec{r}_{IJK} &= [ROT3(-\Omega)][ROT1(-i)][ROT3(-\omega)]\vec{r}_{PQW}
                               = \left [ \frac{IJK}{PQW} \right ]\vec{r}_{PQW}\\
            \vec{v}_{IJK} &= [ROT3(-\Omega)][ROT1(-i)][ROT3(-\omega)]\vec{v}_{PQW}
                               = \left [ \frac{IJK}{PQW} \right ]\vec{v}_{PQW}\\
        \end{align}

    Previous rotations (3-1-3) can be expressed in terms of a single rotation matrix:

    .. math::
        \left [ \frac{IJK}{PQW} \right ]

    .. math::
        \begin{bmatrix}
        \cos(\Omega)\cos(\omega) - \sin(\Omega)\sin(\omega)\cos(i) & -\cos(\Omega)\sin(\omega) - \sin(\Omega)\cos(\omega)\cos(i) & \sin(\Omega)\sin(i)\\
        \sin(\Omega)\cos(\omega) + \cos(\Omega)\sin(\omega)\cos(i) & -\sin(\Omega)\sin(\omega) + \cos(\Omega)\cos(\omega)\cos(i) & -\cos(\Omega)\sin(i)\\
        \sin(\omega)\sin(i) & \cos(\omega)\sin(i) & \cos(i)
        \end{bmatrix}

    """
    pqw = rv_pqw(k, p, ecc, nu)
    rm = coe_rotation_matrix(inc, raan, argp)

    ijk = pqw @ rm.T

    return ijk


def rv2coe(k, r, v, tol=1e-8):
    r"""Convert from vectors to classical orbital elements.

    Parameters
    ----------
    k : float
        Standard gravitational parameter (km^3 / s^2)
    r : array
        Position vector (km)
    v : array
        Velocity vector (km / s)
    tol : float, optional
        Tolerance for eccentricity and inclination checks, default to 1e-8

    Returns
    -------
    p : float
        Semi-latus rectum of parameter (km)
    ecc: float
        Eccentricity
    inc: float
        Inclination (rad)
    raan: float
        Right ascension of the ascending nod (rad)
    argp: float
        Argument of Perigee (rad)
    nu: float
        True Anomaly (rad)

    Notes
    -----
    This example is a real exercise from Orbital Mechanics for Engineering
    students by Howard D.Curtis. This exercise is 4.3 of 3rd. Edition, page 200.

    1. First the angular momentum is computed:

    .. math::
        \vec{h} = \vec{r} \times \vec{v}

    2. With it the eccentricity can be solved:

    .. math::
        \begin{align}
        \vec{e} &= \frac{1}{\mu}\left [ \left ( v^{2} - \frac{\mu}{r}\right ) \vec{r}  - (\vec{r} \cdot \vec{v})\vec{v} \right ] \\
        e &= \sqrt{\vec{e}\cdot\vec{e}} \\
        \end{align}

    3. The node vector line is solved:

    .. math::
        \begin{align}
        \vec{N} &= \vec{k} \times \vec{h} \\
        N &= \sqrt{\vec{N}\cdot\vec{N}}
        \end{align}

    4. The right ascension node is computed:

    .. math::
        \Omega = \left\{ \begin{array}{lcc}
         cos^{-1}{\left ( \frac{N_{x}}{N} \right )} &   if  & N_{y} \geq  0 \\
         \\ 360^{o} -cos^{-1}{\left ( \frac{N_{x}}{N} \right )} &  if & N_{y} < 0 \\
         \end{array}
        \right.

    5. The argument of perigee:

    .. math::
        \omega  = \left\{ \begin{array}{lcc}
         cos^{-1}{\left ( \frac{\vec{N}\vec{e}}{Ne} \right )} &   if  & e_{z} \geq  0 \\
         \\ 360^{o} -cos^{-1}{\left ( \frac{\vec{N}\vec{e}}{Ne} \right )} &  if & e_{z} < 0 \\
         \end{array}
        \right.

    6. And finally the true anomaly:

    .. math::
        \nu  = \left\{ \begin{array}{lcc}
         cos^{-1}{\left ( \frac{\vec{e}\vec{r}}{er} \right )} &   if  & v_{r} \geq  0 \\
         \\ 360^{o} -cos^{-1}{\left ( \frac{\vec{e}\vec{r}}{er} \right )} &  if & v_{r} < 0 \\
         \end{array}
        \right.

    """
    h = cross(r, v)
    n = cross([0, 0, 1], h)
    e = ((v.dot(v) - k / (norm(r))) * r - r.dot(v) * v) / k
    ecc = norm(e)
    p = h.dot(h) / k
    inc = np.arccos(h[2] / norm(h))

    circular = ecc < tol
    equatorial = abs(inc) < tol

    if equatorial and not circular:
        raan = 0
        argp = np.arctan2(e[1], e[0]) % (2 * np.pi)  # Longitude of periapsis
        nu = np.arctan2(h.dot(cross(e, r)) / norm(h), r.dot(e))
    elif not equatorial and circular:
        raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
        argp = 0
        # Argument of latitude
        nu = np.arctan2(r.dot(cross(h, n)) / norm(h), r.dot(n))
    elif equatorial and circular:
        raan = 0
        argp = 0
        nu = np.arctan2(r[1], r[0]) % (2 * np.pi)  # True longitude
    else:
        a = p / (1 - (ecc**2))
        ka = k * a
        if a > 0:
            e_se = r.dot(v) / sqrt(ka)
            e_ce = norm(r) * v.dot(v) / k - 1
            nu = E_to_nu(np.arctan2(e_se, e_ce), ecc)
        else:
            e_sh = r.dot(v) / sqrt(-ka)
            e_ch = norm(r) * (norm(v) ** 2) / k - 1
            nu = H_to_nu(np.log((e_ch + e_sh) / (e_ch - e_sh)) / 2, ecc)

        raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
        px = r.dot(n)
        py = r.dot(cross(h, n)) / norm(h)
        argp = (np.arctan2(py, px) - nu) % (2 * np.pi)

    nu = (nu + np.pi) % (2 * np.pi) - np.pi

    return p, ecc, inc, raan, argp, nu


def rotation_matrix(angle, axis):
    """Return a rotation matrix for a certain angle around an axis.

    Parameters
    ----------
    angle: float
        Angle of rotation.
    axis: int
        Rotation relating 0, 1, and 2 with x, y, and z respectively.

    """
    c = cos(angle)
    s = sin(angle)
    if axis == 0:
        return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
    elif axis == 1:
        return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [s, 0.0, c]])
    elif axis == 2:
        return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    else:
        raise ValueError("Invalid axis: must be one of 'x', 'y' or 'z'")
