"""pytest configuration and shared fixtures."""

import numpy as np
import pytest

from lamberthub import ALL_SOLVERS


def pytest_collection_modifyitems(config, items):
    """Skip benchmark tests unless explicitly requested with ``-m benchmark``."""
    markexpr = getattr(config.option, "markexpr", "")
    if markexpr.strip() == "benchmark":
        return

    skip_mark = pytest.mark.skip(
        reason="Benchmark tests are opt-in: pass -m benchmark to run"
    )
    for item in items:
        if item.get_closest_marker("benchmark"):
            item.add_marker(skip_mark)


@pytest.fixture(scope="session", autouse=True)
def _warmup_jit():
    """Pre-compile all Numba JIT functions before any benchmark timing begins.

    Each solver is called once with a representative set of inputs so that
    JIT compilation (and any Numba cache misses) happen during setup rather
    than during the timed benchmark rounds.
    """
    mu = 3.986004418e5
    r1 = np.array([15945.34, 0.0, 0.0])
    r2 = np.array([12214.83899, 10249.46731, 0.0])
    tof = 76.0 * 60

    # Multi-rev inputs for solvers that support M > 0
    r1_mr = np.array([22592.145603, -1599.915239, -19783.950506])
    r2_mr = np.array([1922.067697, 4054.157051, -8925.727465])
    tof_mr = 36000.0

    for solver in ALL_SOLVERS:
        try:
            solver(mu, r1, r2, tof, M=0, is_prograde=True)
        except Exception:
            pass
        try:
            solver(mu, r1_mr, r2_mr, tof_mr, M=1, is_prograde=True)
        except Exception:
            pass
