"""pytest configuration and shared fixtures."""

import pytest


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
