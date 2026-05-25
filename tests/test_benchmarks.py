"""Benchmarks for fair Lambert solver performance comparisons."""

import numpy as np
import pytest

from lamberthub import (
    ALL_SOLVERS,
    MULTI_REV_SOLVERS,
    ROBUST_SOLVERS,
    ZERO_REV_SOLVERS,
)
from lamberthub.universal_solvers.gooding import gooding1990

REFERENCE_SOLVER = gooding1990
REFERENCE_RTOL = 2e-2


def make_case(
    regime,
    name,
    mu,
    r1,
    r2,
    tof,
    M=0,
    is_prograde=True,
    is_low_path=True,
):
    """Build a benchmark case using the solver call signature."""
    return {
        "regime": regime,
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
                pytest.param(
                    case,
                    solver,
                    id=f"{case['regime']}-{case['name']}-{solver.__name__}",
                )
            )
    return params


ZERO_REV_NOMINAL_CASES = [
    make_case(
        "zero-rev-nominal",
        "vallado-book",
        3.986004418e5,
        [15945.34, 0.0, 0.0],
        [12214.83899, 10249.46731, 0.0],
        76.0 * 60,
    ),
    make_case(
        "zero-rev-nominal",
        "battin-book",
        39.47692641,
        [0.159321004, 0.579266185, 0.052359607],
        [0.057594337, 0.605750797, 0.068345246],
        0.010794065,
    ),
]

ZERO_REV_GENERAL_CASES = [
    make_case(
        "zero-rev-general",
        "curtiss-book",
        3.986004418e5,
        [5000.0, 10000.0, 2100.0],
        [-14600.0, 2500.0, 7000.0],
        3600,
    ),
    make_case(
        "zero-rev-general",
        "der-article-i-prograde-low",
        3.986004418e5,
        [2249.171260, 1898.007100, 5639.599193],
        [1744.495443, -4601.556054, 4043.864391],
        1618.50,
        is_low_path=True,
    ),
    make_case(
        "zero-rev-general",
        "der-article-i-retrograde-high",
        3.986004418e5,
        [2249.171260, 1898.007100, 5639.599193],
        [1744.495443, -4601.556054, 4043.864391],
        1618.50,
        is_prograde=False,
        is_low_path=False,
    ),
    make_case(
        "zero-rev-general",
        "der-article-ii-prograde-high",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        is_low_path=False,
    ),
    make_case(
        "zero-rev-general",
        "der-article-ii-retrograde-high",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        is_prograde=False,
        is_low_path=False,
    ),
]

ZERO_REV_HYPERBOLIC_CASES = [
    make_case(
        "zero-rev-hyperbolic",
        "gmat-prograde",
        3.986004418e5,
        [7100.0, 200.0, 1300.0],
        [-38113.5870, 67274.1946, 29309.5799],
        12000.0,
    ),
    make_case(
        "zero-rev-hyperbolic",
        "gmat-retrograde",
        3.986004418e5,
        [7100.0, 200.0, 1300.0],
        [-47332.7499, -54840.2027, -37100.17067],
        12000.0,
        is_prograde=False,
    ),
]

MULTI_REV_M1_CASES = [
    make_case(
        "multi-rev-m1",
        "prograde-high",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        M=1,
        is_low_path=False,
    ),
    make_case(
        "multi-rev-m1",
        "prograde-low",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        M=1,
        is_low_path=True,
    ),
    make_case(
        "multi-rev-m1",
        "retrograde-high",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        M=1,
        is_prograde=False,
        is_low_path=False,
    ),
    make_case(
        "multi-rev-m1",
        "retrograde-low",
        3.986004418e5,
        [22592.145603, -1599.915239, -19783.950506],
        [1922.067697, 4054.157051, -8925.727465],
        36000,
        M=1,
        is_prograde=False,
        is_low_path=True,
    ),
]

MULTI_REV_HARD_CASES = [
    make_case(
        "multi-rev-hard",
        "m1-prograde-high",
        3.986004418e5,
        [7231.58074563487, 218.02523761425, 11.79251215952],
        [7357.06485698842, 253.55724281562, 38.81222241557],
        12300,
        M=1,
        is_low_path=False,
    ),
    make_case(
        "multi-rev-hard",
        "m1-prograde-low",
        3.986004418e5,
        [7231.58074563487, 218.02523761425, 11.79251215952],
        [7357.06485698842, 253.55724281562, 38.81222241557],
        12300,
        M=1,
        is_low_path=True,
    ),
    make_case(
        "multi-rev-hard",
        "m2-prograde-high",
        3.986004418e5,
        [7231.58074563487, 218.02523761425, 11.79251215952],
        [7357.06485698842, 253.55724281562, 38.81222241557],
        12300,
        M=2,
        is_low_path=False,
    ),
    make_case(
        "multi-rev-hard",
        "m2-prograde-low",
        3.986004418e5,
        [7231.58074563487, 218.02523761425, 11.79251215952],
        [7357.06485698842, 253.55724281562, 38.81222241557],
        12300,
        M=2,
        is_low_path=True,
    ),
]

ROBUST_ZERO_REV_CASES = [
    make_case(
        "robust-zero-rev",
        "near-tangent-prograde-high",
        3.986004418e5,
        [7231.58074563487, 218.02523761425, 11.79251215952],
        [7357.06485698842, 253.55724281562, 38.81222241557],
        12300,
        is_low_path=False,
    ),
]

ZERO_REV_NOMINAL_PARAMS = make_solver_cases(ZERO_REV_NOMINAL_CASES, ALL_SOLVERS)
ZERO_REV_GENERAL_PARAMS = make_solver_cases(ZERO_REV_GENERAL_CASES, ZERO_REV_SOLVERS)
ZERO_REV_HYPERBOLIC_PARAMS = make_solver_cases(
    ZERO_REV_HYPERBOLIC_CASES, ZERO_REV_SOLVERS
)
MULTI_REV_M1_PARAMS = make_solver_cases(MULTI_REV_M1_CASES, MULTI_REV_SOLVERS)
MULTI_REV_HARD_PARAMS = make_solver_cases(MULTI_REV_HARD_CASES, ROBUST_SOLVERS)
ROBUST_ZERO_REV_PARAMS = make_solver_cases(ROBUST_ZERO_REV_CASES, ROBUST_SOLVERS)


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


def assert_solver_matches_reference(result, case):
    """Check a benchmarked solver against the reference solution."""
    expected_v1, expected_v2 = solve_case(REFERENCE_SOLVER, case)
    v1, v2 = result

    assert np.all(np.isfinite(v1))
    assert np.all(np.isfinite(v2))

    actual = np.concatenate((v1, v2))
    expected = np.concatenate((expected_v1, expected_v2))
    relative_error = np.linalg.norm(actual - expected) / np.linalg.norm(expected)
    assert relative_error <= REFERENCE_RTOL


def benchmark_solver(benchmark, case, solver):
    """Benchmark a solver only after confirming it solves the same case."""
    assert_solver_matches_reference(solve_case(solver, case), case)

    result = benchmark(solve_case, solver, case)

    assert_solver_matches_reference(result, case)


@pytest.mark.benchmark(group="zero-rev-nominal")
@pytest.mark.parametrize(("case", "solver"), ZERO_REV_NOMINAL_PARAMS)
def test_zero_revolution_nominal(benchmark, case, solver):
    """Benchmark nominal zero-revolution cases supported by all solvers."""
    benchmark_solver(benchmark, case, solver)


@pytest.mark.benchmark(group="zero-rev-general")
@pytest.mark.parametrize(("case", "solver"), ZERO_REV_GENERAL_PARAMS)
def test_zero_revolution_general(benchmark, case, solver):
    """Benchmark general zero-revolution cases using zero-revolution solvers."""
    benchmark_solver(benchmark, case, solver)


@pytest.mark.benchmark(group="zero-rev-hyperbolic")
@pytest.mark.parametrize(("case", "solver"), ZERO_REV_HYPERBOLIC_PARAMS)
def test_zero_revolution_hyperbolic(benchmark, case, solver):
    """Benchmark hyperbolic zero-revolution cases."""
    benchmark_solver(benchmark, case, solver)


@pytest.mark.benchmark(group="multi-rev-m1")
@pytest.mark.parametrize(("case", "solver"), MULTI_REV_M1_PARAMS)
def test_multi_revolution_m1(benchmark, case, solver):
    """Benchmark first-revolution cases across path and direction branches."""
    benchmark_solver(benchmark, case, solver)


@pytest.mark.benchmark(group="multi-rev-hard")
@pytest.mark.parametrize(("case", "solver"), MULTI_REV_HARD_PARAMS)
def test_multi_revolution_hard(benchmark, case, solver):
    """Benchmark hard multi-revolution cases with robust solvers."""
    benchmark_solver(benchmark, case, solver)


@pytest.mark.benchmark(group="robust-zero-rev")
@pytest.mark.parametrize(("case", "solver"), ROBUST_ZERO_REV_PARAMS)
def test_robust_zero_revolution(benchmark, case, solver):
    """Benchmark difficult zero-revolution cases with robust solvers."""
    benchmark_solver(benchmark, case, solver)
