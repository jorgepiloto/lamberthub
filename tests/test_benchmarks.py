"""Benchmarks for representative Lambert solver test cases."""

import numpy as np
import pytest

from lamberthub import (
    ALL_SOLVERS as ALL_SOLVERS_LAMBERTHUB,
    MULTI_REV_SOLVERS,
    NON_ROBUST_SOLVERS,
    ROBUST_SOLVERS,
)

ALL_SOLVERS = [
    solver for solver in ALL_SOLVERS_LAMBERTHUB if solver not in NON_ROBUST_SOLVERS
]


def make_case(name, mu, r1, r2, tof, M=0, is_prograde=True, is_low_path=True):
    """Build a benchmark case using the solver call signature."""
    return {
        "name": name,
        "mu": mu,
        "r1": np.array(r1),
        "r2": np.array(r2),
        "tof": tof,
        "M": M,
        "is_prograde": is_prograde,
        "is_low_path": is_low_path,
    }


def make_solver_cases(cases, solvers):
    """Flatten cases and solvers into readable pytest parameters."""
    params = []
    for case in cases:
        for solver in solvers:
            params.append(
                pytest.param(case, solver, id=f"{case['name']}-{solver.__name__}")
            )
    return params


ZERO_REV_ALL_SOLVERS_CASES = [
    make_case(
        "vallado-book",
        3.986004418e5,
        [15945.34, 0.0, 0.0],
        [12214.83899, 10249.46731, 0.0],
        76.0 * 60,
    ),
    make_case(
        "battin-book",
        39.47692641,
        [0.159321004, 0.579266185, 0.052359607],
        [0.057594337, 0.605750797, 0.068345246],
        0.010794065,
    ),
]

ZERO_REV_ROBUST_ENOUGH_CASES = [
    make_case(
        "curtiss-book",
        3.986004418e5,
        [5000.0, 10000.0, 2100.0],
        [-14600.0, 2500.0, 7000.0],
        3600,
    ),
    make_case(
        "gmat-hyperbolic-prograde",
        3.986004418e5,
        [7100.0, 200.0, 1300.0],
        [-38113.5870, 67274.1946, 29309.5799],
        12000.0,
    ),
    make_case(
        "gmat-hyperbolic-retrograde",
        3.986004418e5,
        [7100.0, 200.0, 1300.0],
        [-47332.7499, -54840.2027, -37100.17067],
        12000.0,
        is_prograde=False,
    ),
    make_case(
        "der-article-i-prograde-low",
        3.986004418e5,
        [2249.171260, 1898.007100, 5639.599193],
        [1744.495443, -4601.556054, 4043.864391],
        1618.50,
        is_low_path=True,
    ),
    make_case(
        "der-article-i-retrograde-high",
        3.986004418e5,
        [2249.171260, 1898.007100, 5639.599193],
        [1744.495443, -4601.556054, 4043.864391],
        1618.50,
        is_prograde=False,
        is_low_path=False,
    ),
    make_case(
        "der-article-ii-prograde-high",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        is_low_path=False,
    ),
    make_case(
        "der-article-ii-retrograde-high",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        is_prograde=False,
        is_low_path=False,
    ),
]

MULTI_REV_CASES = [
    make_case(
        "m1-prograde-high",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        M=1,
        is_low_path=False,
    ),
    make_case(
        "m1-prograde-low",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        M=1,
        is_low_path=True,
    ),
    make_case(
        "m1-retrograde-high",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        M=1,
        is_prograde=False,
        is_low_path=False,
    ),
    make_case(
        "m1-retrograde-low",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        M=1,
        is_prograde=False,
        is_low_path=True,
    ),
]

ROBUST_CASES = [
    make_case(
        "hard-m0-prograde-high",
        3.986004418e5,
        [7231.58074563487, 218.02523761425, 11.79251215952],
        [7357.06485698842, 253.55724281562, 38.81222241557],
        12300,
        is_low_path=False,
    ),
    make_case(
        "hard-m1-prograde-high",
        3.986004418e5,
        [7231.58074563487, 218.02523761425, 11.79251215952],
        [7357.06485698842, 253.55724281562, 38.81222241557],
        12300,
        M=1,
        is_low_path=False,
    ),
    make_case(
        "hard-m1-prograde-low",
        3.986004418e5,
        [7231.58074563487, 218.02523761425, 11.79251215952],
        [7357.06485698842, 253.55724281562, 38.81222241557],
        12300,
        M=1,
        is_low_path=True,
    ),
    make_case(
        "hard-m2-prograde-high",
        3.986004418e5,
        [7231.58074563487, 218.02523761425, 11.79251215952],
        [7357.06485698842, 253.55724281562, 38.81222241557],
        12300,
        M=2,
        is_low_path=False,
    ),
    make_case(
        "hard-m2-prograde-low",
        3.986004418e5,
        [7231.58074563487, 218.02523761425, 11.79251215952],
        [7357.06485698842, 253.55724281562, 38.81222241557],
        12300,
        M=2,
        is_low_path=True,
    ),
]

ZERO_REV_PARAMS = [
    *make_solver_cases(ZERO_REV_ALL_SOLVERS_CASES, ALL_SOLVERS_LAMBERTHUB),
    *make_solver_cases(ZERO_REV_ROBUST_ENOUGH_CASES, ALL_SOLVERS),
]
MULTI_REV_PARAMS = make_solver_cases(MULTI_REV_CASES, MULTI_REV_SOLVERS)
ROBUST_PARAMS = make_solver_cases(ROBUST_CASES, ROBUST_SOLVERS)


def solve_case(solver, case):
    """Run a solver with one benchmark case."""
    return solver(
        case["mu"],
        case["r1"],
        case["r2"],
        case["tof"],
        M=case["M"],
        is_prograde=case["is_prograde"],
        is_low_path=case["is_low_path"],
    )


def assert_solver_result_is_finite(result):
    """Check the benchmarked solver returned finite velocity vectors."""
    v1, v2 = result
    assert np.all(np.isfinite(v1))
    assert np.all(np.isfinite(v2))


@pytest.mark.benchmark(group="zero-revolution-general")
@pytest.mark.parametrize(("case", "solver"), ZERO_REV_PARAMS)
def test_zero_revolution_general(benchmark, case, solver):
    """Benchmark all general zero-revolution solver test cases."""
    assert_solver_result_is_finite(solve_case(solver, case))

    result = benchmark(solve_case, solver, case)

    assert_solver_result_is_finite(result)


@pytest.mark.benchmark(group="multi-revolution-general")
@pytest.mark.parametrize(("case", "solver"), MULTI_REV_PARAMS)
def test_multi_revolution_general(benchmark, case, solver):
    """Benchmark all general multi-revolution solver test cases."""
    assert_solver_result_is_finite(solve_case(solver, case))

    result = benchmark(solve_case, solver, case)

    assert_solver_result_is_finite(result)


@pytest.mark.benchmark(group="robust-general")
@pytest.mark.parametrize(("case", "solver"), ROBUST_PARAMS)
def test_robust_general(benchmark, case, solver):
    """Benchmark all robust solver test cases."""
    assert_solver_result_is_finite(solve_case(solver, case))

    result = benchmark(solve_case, solver, case)

    assert_solver_result_is_finite(result)
