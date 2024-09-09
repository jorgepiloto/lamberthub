import numpy as np
import pytest

from lamberthub.utils.assertions import (
    assert_gravitational_parameter_is_positive,
    assert_number_of_revolutions_not_negative,
    assert_position_vector_is_valid,
    assert_position_vectors_are_valid,
    assert_time_of_flight_is_positive,
)


def test_exception_if_mu_is_negative():
    with pytest.raises(ValueError) as excinfo:
        assert_gravitational_parameter_is_positive(-1.00)
    assert "Gravitational parameter must be positive!" in excinfo.exconly()


def test_exception_if_vec_not_three_dimensional():
    with pytest.raises(ValueError) as excinfo:
        assert_position_vector_is_valid(np.zeros(4))
    assert "Vector must be three-dimensional!" in excinfo.exconly()


def test_exception_if_vec_is_the_null_vector():
    with pytest.raises(ValueError) as excinfo:
        assert_position_vector_is_valid(np.zeros(3))
    assert "Position vector cannot be the null vector [0,0,0]!" in excinfo.exconly()


def test_exception_position_vectors_have_same_components():
    with pytest.raises(ValueError) as excinfo:
        assert_position_vectors_are_valid(np.ones(3), np.ones(3))
    assert "Initial and final position vectors cannot be equal!" in excinfo.exconly()


def test_exception_if_time_of_flight_is_negative():
    with pytest.raises(ValueError) as excinfo:
        assert_time_of_flight_is_positive(-1.00)
    assert "Time of flight must be positive!" in excinfo.exconly()


def test_exception_if_number_of_revolutions_is_negative():
    with pytest.raises(ValueError) as excinfo:
        assert_number_of_revolutions_not_negative(-1)
    assert (
        "Number of revolutions must be equal or greater than zero!" in excinfo.exconly()
    )
