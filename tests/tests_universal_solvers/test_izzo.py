import pytest

from lamberthub.universal_solvers import izzo


@pytest.mark.parametrize("M", [1, 2, 3])
def test_minimum_time_of_flight_convergence(M):
    ll = -1
    x_T_min_expected, T_min_expected = izzo._compute_T_min(ll, M, maxiter=10, rtol=1e-8)
    y = izzo._compute_y(x_T_min_expected, ll)
    T_min = izzo._tof_equation_y(x_T_min_expected, y, 0.0, ll, M)
    assert T_min_expected == T_min
