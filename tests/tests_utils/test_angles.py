import numpy as np
import pytest
from numpy.testing import assert_allclose

from lamberthub.utils.angles import get_transfer_angle


@pytest.mark.parametrize("sense", [True, False])
def test_get_transfer_angle_collinear_vectors_same_sense(sense):

    # Build two vectors with the same direction and sense using scalar
    # proportion
    r1 = np.array([1, 0, 0])
    r2 = r1 if sense is True else -r1

    # Compute the transfer angle
    dtheta = get_transfer_angle(r1, r2, True)
    expected_dtheta = 0 if sense is True else np.pi

    # Check that transfer angle is zero
    assert_allclose(dtheta, expected_dtheta)
