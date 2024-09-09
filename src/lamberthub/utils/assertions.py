"""A collection of checkers."""

import numpy as np


def assert_parameters_are_valid(mu, r1, r2, tof, M):
    """Check if solver input parameters are safe.

    Parameters
    ----------
    mu : float
        Gravitational parameter.
    r1 : np.array
        Initial position vector.
    r2 : np.array
        Final position vector.
    tof: float
        Time of flight.
    M : int
        Number of revolutions

    """
    assert_gravitational_parameter_is_positive(mu)
    assert_position_vectors_are_valid(r1, r2)
    assert_time_of_flight_is_positive(tof)
    assert_number_of_revolutions_not_negative(M)

    return True


def assert_gravitational_parameter_is_positive(mu):
    """Check if the gravitational parameter is positive.

    Parameters
    ----------
    mu: float
        Gravitational parameter

    Raises
    ------
    ValueError

    """
    # Check positive gravitational parameter
    if mu <= 0:
        raise ValueError("Gravitational parameter must be positive!")
    else:
        return True


def assert_position_vector_is_valid(r):
    """Check if position vector has proper dimensions and is not the null one.

    Parameters
    ----------
    r: np.array
        Initial position vector.

    Raises
    ------
    ValueError

    """
    if r.shape != (3,):
        raise ValueError("Vector must be three-dimensional!")

    if np.all(r == 0):
        raise ValueError("Position vector cannot be the null vector [0,0,0]!")

    return True


def assert_position_vectors_are_valid(r1, r2):
    """Check if position vectors are safe in dimension and values.

    Parameters
    ----------
    r1: np.array
        Initial position vector.
    r2: np.array
        Final position vector.

    Raises
    ------
    ValueError

    """
    # Check if position vectors have proper dimensions
    for r in [r1, r2]:
        assert_position_vector_is_valid(r1)

    # Check both vectors are different
    if np.all(np.equal(r1, r2)):
        raise ValueError("Initial and final position vectors cannot be equal!")

    return True


def assert_time_of_flight_is_positive(tof):
    """Check if time of flight is positive.

    Parameters
    ----------
    tof: float
        Time of flight.

    Raises
    ------
    ValueError

    """
    if tof <= 0:
        raise ValueError("Time of flight must be positive!")
    else:
        return True


def assert_number_of_revolutions_not_negative(M):
    """Check if the number of revolutions is not negative.

    Parameters
    ----------
    M: int
        Number of revolutions

    Raises
    ------
    ValueError

    """
    if M < 0:
        raise ValueError("Number of revolutions must be equal or greater than zero!")
    else:
        return True


def assert_transfer_angle_not_zero(dtheta):
    """Check if the transfer angle is the null value.

    Parameters
    ----------
    dtheta: float
        Transfer angle value.

    Raises
    ------
    ValueError

    """
    if dtheta == 0:
        raise ValueError("Transfer angle was found to be zero!")
    else:
        return True


def assert_transfer_angle_not_pi(dtheta):
    """Check if the transfer angle is pi radians or 180 degrees.

    Parameters
    ----------
    dtheta: float
        Transfer angle value.

    Raises
    ------
    ValueError

    """
    if dtheta == np.pi:
        raise ValueError("Transfer angle was found to be 180 degrees!")
    else:
        return True
