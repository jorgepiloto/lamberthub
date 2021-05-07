import numpy as np
import pytest

from lamberthub import gooding1990


def test_transfer_angle_is_zero_raises_exception():
    with pytest.raises(ValueError) as excinfo:
        r1, r2 = [i * np.ones(3) for i in range(1, 3)]
        gooding1990(1.00, r1, r2, 1.00)
    assert "Transfer angle was found to be zero!" in excinfo.exconly()
