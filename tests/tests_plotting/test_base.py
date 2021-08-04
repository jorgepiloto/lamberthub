import numpy as np
from numpy.testing import assert_allclose

from lamberthub.p_solvers.gauss import gauss1809
from lamberthub.plotting._base import _measure_performance


def test_measure_performance_handles_iteration_exceptions():
    results = _measure_performance(gauss1809, np.pi / 2, 2 * np.pi)
    for value in results:
        assert_allclose(value, 0.0)
