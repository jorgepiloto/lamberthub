"""
Holds auxiliary functions used by eccentricity based Lambert's problem
solvers. The majority of the routines hosted within this module were directly
taken from the following publications [1] [2] [3].

[1] Avanzini, G. (2008). A simple Lambert algorithm. Journal of guidance,
    control, and dynamics, 31(6), 1587-1594.

[2] He, Q., Li, J., & Han, C. (2010). Multiple-revolution solutions of the
    transverse-eccentricity-based Lambert problem. Journal of guidance,
    control, and dynamics, 33(1), 265-269.

[3] Wen, Changxuan, Yushan Zhao, and Peng Shi. "Derivative analysis and
    algorithm modification of transverse-eccentricity-based Lambert problem."
    Journal of Guidance, Control, and Dynamics 37.4 (2014): 1195-1201.

"""

import numpy as np
from numpy import cross
from numpy.linalg import norm

from lamberthub.utils.angles import (
    get_orbit_inc_and_raan_from_position_vectors,
    get_orbit_normal_vector,
    get_transfer_angle,
)
from lamberthub.utils.assertions import assert_transfer_angle_not_zero
from lamberthub.utils.kepler import kepler_from_nu


def get_geometry(r1, r2, prograde):
    """
    Computes associated problem geometry.

    Parameters
    ----------
    r1: np.array
        Initial position vector.
    r2: np.array
        Final position vector.
    prograde: bool
        If True, assumes prograde motion, otherwise retrograde is applied.

    Returns
    -------
    r1_norm: float
        Norm of the initial position vector.
    r1_norm: float
        Norm of the final position vector.
    c_norm: float
        Norm of the chord vector.
    dtheta: float
        Transfer angle.
    w_c: float
        Angle between initial position and chord vectors.

    """
    # Solve for the norms
    r1_norm, r2_norm, c_norm = [np.linalg.norm(vec) for vec in [r1, r2, (r2 - r1)]]

    # Compute angles
    dtheta = get_transfer_angle(r1, r2, prograde)
    assert_transfer_angle_not_zero(dtheta)
    w_c = get_transfer_angle(r1, (r2 - r1), prograde)

    return r1_norm, r2_norm, c_norm, dtheta, w_c


def get_eccF(r1_norm, r2_norm, c_norm):
    """
    Computes the eccentricity component along the chord. This value is kept
    constant for all the problem as long as the boundary conditions are not
    changed.

    Parameters
    ----------
    r1_norm: float
        Norm of the initial vector position.
    r2_norm: float
        Norm of the final vector position.
    c_norm: float
        Norm of the chord vector.

    Returns
    -------
    ecc_F: float
        Eccentricity component along the chord direction.

    Notes
    -----
    Equation (3) from Avanzini's report [1].

    """
    ecc_F = (r1_norm - r2_norm) / c_norm
    return ecc_F


def get_aF(r1_norm, r2_norm):
    """
    Computes the semi-major axis of the fundamental ellipse. This value is
    kept constant for all the problem as long as the boundary conditions are not
    changed.

    Parameters
    ----------
    r1_norm: float
        Norm of the initial vector position.
    r2_norm: float
        Norm of the final vector position.

    Returns
    -------
    a_F: float
        Semi-major axis of the fundamental ellipse.

    Notes
    -----
    No labeled equation (appears between [3] and [4]) from Avanzini's report
    [1].

    """
    a_F = (r1_norm + r2_norm) / 2
    return a_F


def get_pF(a_F, ecc_F):
    """
    Computes the orbital parameter (semi-latus) rectum of the fundamental
    ellipse. This value is kept constant for all the problem as long as the
    boundary conditions are not changed.

    Parameters
    ----------
    a_F: float
        Semi-major axis of the fundamental ellipse.
    ecc_F: float
        Eccentricity of the fundamental ellipse.

    Returns
    -------
    p_F: float
        Orbital parameter / semi-latus rectum of the fundamental ellipse.

    Notes
    -----
    No labeled equation (appears between [3] and [4]) from Avanzini's report

    """
    p_F = a_F * (1 - ecc_F**2)
    return p_F


def get_fundamental_ellipse_properties(r1_norm, r2_norm, c_norm):
    """
    Computes the fundamental ellipse properties. Those are the eccentricity,
    semi-major axis and the orbital parameter.

    Parameters
    ----------
    r1_norm: float
        Norm of the initial vector position.
    r2_norm: float
        Norm of the final vector position.
    c_norm: float
        Norm of the chord vector.

    Returns
    -------
    ecc_F: float
        Eccentricity component along the chord direction.
    a_F: float
        Semi-major axis of the fundamental ellipse.
    p_F: float
        Orbital parameter / semi-latus rectum of the fundamental ellipse.

    """
    # Compute the fundamental ellipse parameters
    ecc_F = get_eccF(r1_norm, r2_norm, c_norm)
    a_F = get_aF(r1_norm, r2_norm)
    p_F = get_pF(a_F, ecc_F)

    return ecc_F, a_F, p_F


def ecc_at_eccT(ecc_T, ecc_F):
    """
    Computes transfer orbit eccentricity from transverse and fundamental
    components.

    Parameters
    ----------
    ecc_T: float
        Eccentricity transverse component.
    ecc_F: float
        Eccentricity of the fundamental ellipse.

    Returns
    -------
    ecc: float
        Eccentricity of the transfer orbit.

    """
    ecc = np.sqrt(ecc_T**2 + ecc_F**2)
    return ecc


def p_at_eccT(ecc_T, r1_norm, r2_norm, c_norm, dtheta, p_F):
    """
    Computes the orbital parameter or semi-latus rectum of the transfer orbit.

    Parameters
    ----------
    ecc_T: float
        Eccentricity transverse component.
    r1_norm: float
        Norm of the initial vector position.
    r2_norm: float
        Norm of the final vector position.
    c_norm: float
        Norm of the chord vector.
    dtheta: float
        Transfer angle.
    p_F: float
        Orbital parameter or semi-latus rectum of the fundamental ellipse.

    Returns
    -------
    p: float
        Orbital parameter or semi-lactus rectum.

    """
    p = p_F - ecc_T * r1_norm * r2_norm * np.sin(dtheta) / c_norm
    return p


def a_at_eccT(ecc_T, ecc_F, p):
    """
    Computes the semi-major axis of the transfer orbit.

    Parameters
    ----------
    ecc_T: float
        Eccentricity transverse component.
    ecc_F: float
        Eccentricity of the fundamental ellipse.
    p: float
        Transfer orbit parameter or semi-latus rectum.

    Returns
    -------
    a: float
        Semi-major axis of the transfer orbit.

    """
    a = p / (1 - ecc_F**2 - ecc_T**2)
    return a


def eap_from_eccT(ecc_T, geometry):
    """
    Solves for transfer orbit eccentricity, semi-major axis and orbital
    parameter.

    Parameters
    ----------
    ecc_T: float
        Eccentricity component along transverse direction.

    geometry: tuple
        A tuple hosting r1_norm, r2_norm, c_norm, dtheta and w_c geometry values.

    Returns
    -------
    ecc: float
        Absolute eccentricity of the transfer orbit.
    a: float
        Semi-major axis of the transfer orbit.
    p: float
        Semi-latus rectum of the transfer orbit.

    """
    # Unpack useful parameters
    r1_norm, r2_norm, c_norm, dtheta, w_c = geometry

    # Solve for the fundamental ellipse properties
    ecc_F, a_F, p_F = get_fundamental_ellipse_properties(r1_norm, r2_norm, c_norm)

    # Compute the transfer orbit eccentricity, semi-latus rectum and
    # semi-major axis.
    ecc = ecc_at_eccT(ecc_T, ecc_F)
    p = p_at_eccT(ecc_T, r1_norm, r2_norm, c_norm, dtheta, p_F)
    a = a_at_eccT(ecc_T, ecc_F, p)

    return ecc, a, p


def w_at_eccT(ecc_T, ecc_F, w_c):
    """
    Compute the true anomalies for the initial and final position vectors
    with respect to the transfer orbit.

    Parameters
    ----------
    ecc_T: float
        Eccentricity transverse component.
    ecc_F: float
        Eccentricity of the fundamental ellipse.
    dtheta: float
        Transfer angle.
    w_c: float
        Angle between the initial and chord vector.

    Returns
    -------
    nu_1: float
        True anomaly of the initial position vector w.r.t. transfer orbit.
    nu_2: float
        True anomaly of the final position vector w.r.t. transfer orbit.

    Notes
    -----
    This is equation (6) from Quan He's report [2].

    """
    # Compute the coordinates
    y = ecc_F * np.sin(w_c) + ecc_T * np.cos(w_c)
    x = ecc_F * np.cos(w_c) - ecc_T * np.sin(w_c)
    w = np.arctan2(y, x)

    return w


def get_true_anomalies(w, dtheta):
    """
    Compute the initial and final true anomalies.

    Parameters
    ----------
    w: float
        Argument of periapsis.
    dtheta: float
        Transfer angle.

    Returns
    -------
    nu_1: float
        Initial true anomaly.
    nu_2: float
        Final true anomaly.

    """
    nu_1 = -w
    nu_2 = nu_1 + dtheta
    return nu_1, nu_2


def kepler_tof_at_eccT(ecc_T, mu, geometry):
    """
    Computes the time of flight at particular value of transverse eccentricity
    and problem boundary conditions.

    Parameters
    ----------
    ecc_T: float
        Eccentricity component along transverse direction.
    mu: float
        The gravitational parameter.
    geometry: tuple
        A tuple hosting r1_norm, r2_norm, c_norm, dtheta and w_c geometry values.

    Returns
    -------
    tof: float
        Dimensional time of flight from Kepler's equation.

    """
    # Unpack useful parameters
    r1_norm, r2_norm, c_norm, dtheta, w_c = geometry

    # Compute the transfer orbit eccentricity, semi-latus rectum and
    # semi-major axis.
    ecc, a, p = eap_from_eccT(ecc_T, geometry)
    ecc_F, _, _ = get_fundamental_ellipse_properties(r1_norm, r2_norm, c_norm)

    # Last step before solving Kepler's equation is to compute the true
    # anomalies. Hence, the argument of periapsis is required to do so.
    w = w_at_eccT(ecc_T, ecc_F, w_c)
    nu_1, nu_2 = get_true_anomalies(w, dtheta)

    # Compute the change in mean motion from corresponding Kepler's
    # equation, that is elliptical, parabolic or hyperbolic.
    M_1, M_2 = [kepler_from_nu(nu, ecc) for nu in [nu_1, nu_2]]
    deltaM = M_2 - M_1

    # Compute the time of flight for that particular value of mean motion
    # variation.
    if 0 <= ecc < 1:
        tof12 = np.sqrt(a**3 / mu) * deltaM
        if tof12 < 0:
            # Add an orbital period so the time of flight becomes positive
            tof12 += 2 * np.pi * np.sqrt(a**3 / mu)
    elif ecc == 1:
        tof12 = (1 / 2) * np.sqrt(p**3 / mu) * deltaM
        if tof12 < 0:
            return np.inf
    else:
        tof12 = np.sqrt(-(a**3) / mu) * deltaM

    return tof12


def _f(x, ecc_T_at_x, mu, geometry, tof12_s):
    """
    Returns a zero once the value of x makes the numerically compute time of
    flight to be exactly the desired one.

    Parameters
    ----------
    x: float
        The independent variable to be solved.
    ecc_T_at_x: function
        Function to solve the eccentricity component along transverse direction.
    mu: float
        The gravitational parameter.
    geometry: tuple
        A tuple hosting r1_norm, r2_norm, c_norm, dtheta and w_c geometry values.
    tof12_s: float
        Desired time of flight to be achieved.

    Notes
    -----
    This is equation (14) from Avanzini's report [1].

    """
    # Compute the predicted time of flight at particular ecc_T
    ecc_T = ecc_T_at_x(x)
    tof12 = kepler_tof_at_eccT(ecc_T, mu, geometry)

    # Compute the non-dimensional times of flight
    t_ref = np.sqrt(geometry[0] ** 3 / mu)
    tau12_s, tau12 = [t / t_ref for t in [tof12_s, tof12]]

    # Compute the values of y
    y, y_s = [np.log(tau) for tau in [tau12, tau12_s]]

    return y - y_s


def coe_at_eccT(ecc_T, r1, r2, sense):
    """
    Computes the classical orbita elements at particular value of transverse
    eccentricity.

    Parameters
    ----------
    ecc_T: float
        Transverse eccentricity.
    r1: np.array
        Initial position vecor.
    r1: np.array
        Final position vector.
    sense: bool
        If True assumes prograde motion, otherwise retrograde is applied.

    Returns
    -------
    p:float
        Orbital parameter or semi-lactus rectum.
    ecc: float
        Absolute orbit eccentricity.
    inc: float
        Inclination of the orbit.
    raan: float
        Right ascension of the ascending node.
    argp: float
        Argument of periapsis.
    nu_1: float
        True anomaly at the first initial position vector.
    nu_2: float
        True anomaly at the first final position vector.
    """
    # Retrieve auxiliary geometry parameters
    geometry = get_geometry(r1, r2, sense)
    r1_norm, r2_norm, c_norm, dtheta, w_c = geometry

    # Orbit parameters as function of actual ecc_T
    ecc, a, p = eap_from_eccT(ecc_T, geometry)
    ecc_F, _, _ = get_fundamental_ellipse_properties(r1_norm, r2_norm, c_norm)
    w = w_at_eccT(ecc_T, ecc_F, w_c)
    nu_1, nu_2 = get_true_anomalies(w, dtheta)

    # Compute the orbit inclination and RAAN
    i_r1, i_r2 = [r / r_norm for r, r_norm in zip([r1, r2], [r1_norm, r2_norm])]

    # Compute a normal vector normal to orbit plane with proper sense
    i_h = get_orbit_normal_vector(r1, r2, sense)

    # Solve for inclination and RAAN
    inc, raan = get_orbit_inc_and_raan_from_position_vectors(r1, r2, sense)

    # Get a vector in the direction of the line of nodes
    n = np.cross(np.array([0, 0, 1]), i_h)
    px = r1.dot(n)
    py = r1.dot(cross(i_h, n)) / norm(i_h)
    argp = (np.arctan2(py, px) - nu_1) % (2 * np.pi)

    return p, ecc, inc, raan, argp, nu_1, nu_2
