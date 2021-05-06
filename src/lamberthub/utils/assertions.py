""" A collection of checkers for raising custom exceptions if required """

import numpy as np


def assert_parameters_are_valid(mu, r1, r2, tof, M):
    """
    Checks if solver input parameters are safe.

    """

    # Run all required checks
    assert_gravitational_parameter_is_positive(mu)
    assert_position_vectors_are_valid(r1, r2)
    assert_time_of_flight_is_positive(tof)
    assert_number_of_revolutions_not_negative(M)

    return True


def assert_gravitational_parameter_is_positive(mu):
    """
    Checks if the gravitational parameter is positive.

    Parameters
    ----------
    mu: float
        Gravitational parameter

    """

    # Check positive gravitational parameter
    if mu <= 0:
        raise ValueError("Gravitational parameter must be positive!")
    else:
        return True


def assert_position_vector_is_valid(r):
    """
    Checks if position vector has proper dimensions and is not the null one.

    Parameters
    ----------
    r: np.array
        Initial position vector.

    """

    # Check that vector belongs to three-dimensional space
    if r.shape != (3,):
        raise ValueError("Vector must be three-dimensional!")

    if np.all(r == 0):
        raise ValueError("Position vector cannot be the null vector [0,0,0]!")

    return True


def assert_position_vectors_are_valid(r1, r2):
    """
    Checks if position vectors are safe in dimension and values.

    Parameters
    ----------
    r1: np.array
        Initial position vector.
    r2: np.array
        Final position vector.

    """

    # Check if position vectors have proper dimensions
    for r in [r1, r2]:
        assert_position_vector_is_valid(r1)

    # Check both vectors are different
    if np.all(np.equal(r1, r2)):
        raise ValueError("Initial and final position vectors cannot be equal!")

    return True


def assert_time_of_flight_is_positive(tof):
    """
    Checks if time of flight is positive.

    Parameters
    ----------
    tof: float
        Time of flight.

    """

    if tof <= 0:
        raise ValueError("Time of flight must be positive!")
    else:
        return True


def assert_number_of_revolutions_not_negative(M):
    """
    Checks if the number of revolutions is zero or positive, that is, it does
    not have a negative value.

    Parameters
    ----------
    M: int
        Number of revolutions

    """

    if M < 0:
        raise ValueError("Number of revolutions must be equal or greater than zero!")
    else:
        return True


def assert_transfer_angle_not_zero(dtheta):
    """
    Checks if the transfer angle is the null value, if so, raises an exception.

    Parameters
    ----------
    dtheta: float
        Transfer angle value.

    """

    if dtheta == 0:
        raise ValueError("Transfer angle was found to be zero!")
    else:
        return True
