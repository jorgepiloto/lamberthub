import numpy as np
import pytest

from lamberthub import avanzini2008


def test_transfer_angle_is_zero_raises_exception():
    with pytest.raises(ValueError) as excinfo:
        r1, r2 = [i * np.ones(3) for i in range(1, 3)]
        avanzini2008(1.00, r1, r2, 1.00)
    assert "Transfer angle was found to be zero!" in excinfo.exconly()


def test_multi_revolution_not_possible():
    with pytest.raises(ValueError) as excinfo:
        r1, r2 = [i * np.ones(3) for i in range(1, 3)]
        avanzini2008(1.00, r1, r2, 1.00, M=2)
    assert (
        "Avanzini is not able to work within the multi-revolution scenario!"
        in excinfo.exconly()
    )
