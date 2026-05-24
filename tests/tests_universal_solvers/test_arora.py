import numpy as np
import pytest

from lamberthub import arora2013


def test_transfer_angle_is_zero_raises_exception():
    with pytest.raises(ValueError) as excinfo:
        r1, r2 = [i * np.ones(3) for i in range(1, 3)]
        arora2013(1.00, r1, r2, 1.00)
    assert "Transfer angle was found to be zero!" in excinfo.exconly()


def test_exception_try_lower_M_multi_revolution_scenario():
    with pytest.raises(ValueError) as excinfo:
        r1 = np.array([1.0, 0.0, 0.0])
        r2 = np.array([0.0, 1.0, 0.0])
        arora2013(1.00, r1, r2, 1.00, M=1)
    assert "No feasible solution, try lower M!" in excinfo.exconly()
