import numpy as np
import pytest

from lamberthub import arora2013


def test_transfer_angle_is_zero_raises_exception():
    with pytest.raises(ValueError) as excinfo:
        r1, r2 = [i * np.ones(3) for i in range(1, 3)]
        arora2013(1.00, r1, r2, 1.00)
    assert "Transfer angle was found to be zero!" in excinfo.exconly()


def test_not_implemented_error_multi_revolution_scenario():
    with pytest.raises(NotImplementedError) as excinfo:
        r1, r2 = [i * np.ones(3) for i in range(1, 3)]
        arora2013(1.00, r1, r2, 1.00, M=1)
    assert "See https://github.com/jorgepiloto/lamberthub/issues/3" in excinfo.exconly()
