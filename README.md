# lamberthub: a hub of Lambert's problem solvers

<img align="left" width=350px src="https://github.com/jorgepiloto/lamberthub/raw/main/doc/source/_static/lamberts_problem_geometry.png"/>

A Python library designed to provide solutions to Lambert's problem, a
classical problem in astrodynamics that involves determining the orbit of a
spacecraft given two points in space and the time of flight between them. The
problem is essential for trajectory planning, particularly for interplanetary
missions.

This library implements multiple algorithms, each named after its author and
publication year, for solving different variations of Lambert's problem. These
algorithms can handle different types of orbits, including multi-revolution
paths and direct transfers.

<br>

<!-- vale off -->
[![Python](https://img.shields.io/pypi/pyversions/lamberthub?logo=pypi)](https://pypi.org/project/lamberthub/)
[![PyPI](https://img.shields.io/pypi/v/lamberthub.svg?logo=python&logoColor=white)](https://pypi.org/project/lamberthub/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CI](https://github.com/jorgepiloto/lamberthub/actions/workflows/ci_cd.yml/badge.svg?branch=main)](https://github.com/jorgepiloto/lamberthub/actions/workflows/ci_cd.yml)
[![Coverage](https://codecov.io/gh/jorgepiloto/lamberthub/branch/main/graph/badge.svg?token=3BY2J5AB8D)](https://codecov.io/gh/jorgepiloto/lamberthub)
[![DOI](https://zenodo.org/badge/364482782.svg)](https://zenodo.org/badge/latestdoi/364482782)
<!-- vale on -->

## Installation

Multiple installation methods are supported:

|                             **Logo**                              | **Platform** |                                    **Command**                                    |
|:-----------------------------------------------------------------:|:------------:|:---------------------------------------------------------------------------------:|
|       ![PyPI logo](https://simpleicons.org/icons/pypi.svg)        |     PyPI     |                        ``python -m pip install lamberthub``                        |
|     ![GitHub logo](https://simpleicons.org/icons/github.svg)      |    GitHub    | ``python -m pip install https://github.com/jorgepiloto/lamberthub/archive/main.zip`` |

## Available solvers

<!-- vale off -->

| Algorithm     | Reference                                                                                                                                               |
|---------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|
| `gauss1809`   | C. F. Gauss, *Theoria motus corporum coelestium in sectionibus conicis solem ambientium*. 1809.                                                         |
| `battin1984`  | R. H. Battin and R. M. Vaughan, “An elegant lambert algorithm,” *Journal of Guidance, Control, and Dynamics*, vol. 7, no. 6, pp. 662–670, 1984.         |
| `gooding1990` | R. Gooding, “A procedure for the solution of lambert’s orbital boundary-value problem,” *Celestial Mechanics and Dynamical Astronomy*, vol. 48, no. 2, pp. 145–165, 1990. |
| `der2011`     | G. J. Der, “The superior Lambert algorithm,” AMOS, Maui, Hawaii, 2011.                                                                                 |
| `thorne2004`  | J. D. Thorne, “Lambert’s theorem: a complete series solution,” *The Journal of the Astronautical Sciences*, vol. 52, no. 4, pp. 441–454, 2004.        |
| `avanzini2008`| G. Avanzini, “A simple lambert algorithm,” *Journal of Guidance, Control, and Dynamics*, vol. 31, no. 6, pp. 1587–1594, 2008.                          |
| `arora2013`   | N. Arora and R. P. Russell, “A fast and robust multiple revolution lambert algorithm using a cosine transformation,” Paper AAS, vol. 13, p. 728, 2013.  |
| `vallado2013` | D. A. Vallado, *Fundamentals of astrodynamics and applications*. Springer Science & Business Media, 2013, vol. 12.                                       |
| `izzo2015`    | D. Izzo, “Revisiting lambert’s problem,” *Celestial Mechanics and Dynamical Astronomy*, vol. 121, no. 1, pp. 1–15, 2015.                                |
| `jiang2016`   | R. Jiang, T. Chao, S. Wang, and M. Yang, “Improved semi-major Axis iterated method for Lambert's problem,” 2016 IEEE Chinese Guidance, Navigation and Control Conference, pp. 1423–1428, 2016. |
| `pan2016`     | B. Pan and Y. Ma, “Lambert’s problem and solution by non-rational Bézier functions,” *Proceedings of the Institution of Mechanical Engineers, Part G: Journal of Aerospace Engineering*, first published online Nov. 16, 2016. |
| `delatorre2018` | D. De La Torre, R. Flores, and E. Fantino, “On the solution of Lambert's problem by regularization,” *Acta Astronautica*, vol. 153, pp. 26–38, 2018. |
| `negrete2024` | A. Negrete and O. Abdelkhalik, “An Exact Solution to Lambert's Problem Using Contour Integrals.”                                                          |
| `mcelreath2025` | J. McElreath, I. M. Down, and M. Majji, “A universal approach for solving the multi-revolution Lambert's problem,” *Celestial Mechanics and Dynamical Astronomy*, vol. 137, 22, 2025. |

<!-- vale on -->

## Using a solver

Any Lambert's problem algorithm implemented in `lamberthub` is a Python function
which accepts the following parameters:

```python
from lamberthub import authorYYYY


v1, v2 = authorYYYY(
    mu, r1, r2, tof, M=0, is_prograde=True, is_low_path=True,  # Type of solution
    maxiter=35, atol=1e-5, rtol=1e-7, full_output=False  # Iteration config
)
```

where `author` is the name of the author which developed the solver and `YYYY`
the year of publication. Any of the solvers hosted by the `ALL_SOLVERS` list.

### Parameters

| Parameters    | Type      | Description |
|---------------|-----------|-------------|
| `mu`          | `float`   | The gravitational parameter, that is, the mass of the attracting body times the gravitational constant. |
| `r1`          | `np.array`| Initial position vector. |
| `r2`          | `np.array`| Final position vector. |
| `tof`         | `float`   | Time of flight between initial and final vectors. |
| `M`           | `int`     | The number of revolutions. If zero (default), direct transfer is assumed. |
| `is_prograde`    | `bool`    | Controls the inclination of the final orbit. If `True`, inclination between 0 and 90 degrees. If `False`, inclination between 90 and 180 degrees. |
| `is_low_path`    | `bool`    | Selects the type of path when more than two solutions are available. No specific advantage unless there are mission constraints. |
| `maxiter`     | `int`     | Maximum number of iterations allowed when computing the solution. |
| `atol`        | `float`   | Absolute tolerance for the iterative method. |
| `rtol`        | `float`   | Relative tolerance for the iterative method. |
| `full_output` | `bool`    | If `True`, returns additional information such as the number of iterations. |

### Returns

| Returns       | Type       | Description |
|---------------|------------|-------------|
| `v1`          | `np.array` | Initial velocity vector. |
| `v2`          | `np.array` | Final velocity vector. |
| `numiter`     | `int`      | Number of iterations (only if `full_output` is `True`). |
| `tpi`         | `float`    | Time per iteration (only if `full_output` is `True`). |

## Examples

### Example: solving for a direct and prograde transfer orbit

**Problem statement**

Suppose you want to solve for the orbit of an interplanetary vehicle (that is
Sun is the main attractor) form which you know that the initial and final
positions are given by:

```math
\vec{r_1} = \begin{bmatrix} 0.159321004 \\ 0.579266185 \\ 0.052359607 \end{bmatrix} \text{ [AU]} \quad \quad
\vec{r_2} = \begin{bmatrix} 0.057594337 \\ 0.605750797 \\ 0.068345246 \end{bmatrix} \text{ [AU]} \quad \quad
```

<br>

The time of flight is $\Delta t = 0.010794065$ years. The orbit is
prograde and direct, thus $M=0$. Remember that when $M=0$, there is only one
possible solution, so the `is_low_path` flag does not play any role in this
problem.

**Solution**

For this problem, `gooding1990` is used. Any other solver would work too. Next,
the parameters of the problem are instantiated. Finally, the initial and final
velocity vectors are computed.

```python
from lamberthub import gooding1990
import numpy as np


mu_sun = 39.47692641
r1 = np.array([0.159321004, 0.579266185, 0.052359607])
r2 = np.array([0.057594337, 0.605750797, 0.068345246])
tof = 0.010794065

v1, v2 = gooding1990(mu_sun, r1, r2, tof, M=0, is_prograde=True)
print(f"Initial velocity: {v1} [AU / years]")
print(f"Final velocity:   {v2} [AU / years]")
```

**Result**

```console
Initial velocity: [-9.303608  3.01862016  1.53636008] [AU / years]
Final velocity:   [-9.511186  1.88884006  1.42137810] [AU / years]
```

Directly taken from *An Introduction to the Mathematics and Methods of
Astrodynamics, revised edition, by R.H. Battin, problem 7-12*.

## Performance comparison

Benchmarks run with
[pytest-benchmark](https://pytest-benchmark.readthedocs.io). All Numba JIT
functions are pre-compiled once before timing begins, compilation cost is never
included in the measurements. Six test groups cover the main solver regimes:

| Case group | What it tests |
|------------|---------------|
| Zero-revolution nominal cases | Classic textbook cases, all solvers |
| Zero-revolution general cases | General elliptic, zero-revolution solvers |
| Zero-revolution hyperbolic cases | Hyperbolic transfers |
| One-revolution branch cases | First multi-revolution branch across prograde/retrograde and low/high paths |
| Near-minimum-energy multi-revolution cases | Near-minimum-energy cases with one and two revolutions |
| Near-tangent robustness cases | Near-tangent transfers that stress robustness |

<!-- performance-comparison:start -->

_This section is auto-generated from CI benchmark artifacts._

- Generated: 2026-06-14 06:39 UTC
- Commit: `5fff8c99d908`
- Environment: Linux, Python 3.12.13, AMD EPYC 9V74 80-Core Processor

Times are in microseconds (lower is better).
**Speedup** is relative to the slowest solver for each case.
All solvers are JIT-warmed before timing begins.

### Near-minimum-energy multi-revolution cases

#### Near-minimum-energy

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 1 | Yes | high | 46.8 | 47.8 | 0.5 | 13.9x | 4332 |
| 2 | `izzo2015` | 1 | Yes | high | 61.0 | 62.2 | 1.0 | 10.7x | 3351 |
| 3 | `gooding1990` | 1 | Yes | high | 102.4 | 104.7 | 2.0 | 6.4x | 1992 |
| 4 | `negrete2024` | 1 | Yes | high | 194.7 | 197.4 | 3.8 | 3.3x | 1034 |
| 5 | `arora2013` | 1 | Yes | high | 213.2 | 215.6 | 7.1 | 3.1x | 969 |
| 6 | `der2011` | 1 | Yes | high | 314.0 | 320.9 | 10.1 | 2.1x | 649 |
| 7 | `delatorre2018` | 1 | Yes | high | 651.5 | 651.9 | 9.9 | 1.0x | 315 |
| 1 | `mcelreath2025` | 1 | Yes | low | 47.3 | 48.1 | 0.5 | 14.9x | 4292 |
| 2 | `izzo2015` | 1 | Yes | low | 61.5 | 62.6 | 0.9 | 11.5x | 3319 |
| 3 | `gooding1990` | 1 | Yes | low | 102.5 | 104.6 | 1.9 | 6.9x | 1981 |
| 4 | `negrete2024` | 1 | Yes | low | 194.9 | 197.1 | 2.8 | 3.6x | 1033 |
| 5 | `arora2013` | 1 | Yes | low | 233.1 | 238.1 | 9.9 | 3.0x | 872 |
| 6 | `der2011` | 1 | Yes | low | 421.0 | 423.1 | 10.5 | 1.7x | 486 |
| 7 | `delatorre2018` | 1 | Yes | low | 704.6 | 707.3 | 7.7 | 1.0x | 291 |
| 1 | `mcelreath2025` | 2 | Yes | high | 47.5 | 48.4 | 0.5 | 13.7x | 4300 |
| 2 | `izzo2015` | 2 | Yes | high | 61.4 | 62.8 | 1.0 | 10.6x | 3343 |
| 3 | `gooding1990` | 2 | Yes | high | 104.0 | 106.3 | 2.0 | 6.3x | 1954 |
| 4 | `negrete2024` | 2 | Yes | high | 196.1 | 202.7 | 8.4 | 3.3x | 1027 |
| 5 | `arora2013` | 2 | Yes | high | 215.3 | 221.5 | 10.9 | 3.0x | 957 |
| 6 | `der2011` | 2 | Yes | high | 267.9 | 271.3 | 9.9 | 2.4x | 758 |
| 7 | `delatorre2018` | 2 | Yes | high | 650.2 | 651.0 | 9.4 | 1.0x | 316 |
| 1 | `mcelreath2025` | 2 | Yes | low | 47.2 | 48.4 | 0.5 | 14.4x | 4315 |
| 2 | `izzo2015` | 2 | Yes | low | 61.9 | 63.0 | 0.9 | 11.0x | 3303 |
| 3 | `gooding1990` | 2 | Yes | low | 103.8 | 106.5 | 2.0 | 6.6x | 1960 |
| 4 | `negrete2024` | 2 | Yes | low | 196.3 | 201.4 | 7.2 | 3.5x | 1027 |
| 5 | `arora2013` | 2 | Yes | low | 229.3 | 236.1 | 10.9 | 3.0x | 895 |
| 6 | `der2011` | 2 | Yes | low | 267.9 | 271.4 | 9.6 | 2.5x | 764 |
| 7 | `delatorre2018` | 2 | Yes | low | 680.2 | 681.3 | 9.4 | 1.0x | 299 |

### One-revolution branch cases

#### Der article II

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 1 | Yes | high | 47.5 | 48.7 | 0.7 | 13.5x | 4324 |
| 2 | `izzo2015` | 1 | Yes | high | 61.4 | 62.7 | 1.0 | 10.5x | 3319 |
| 3 | `gooding1990` | 1 | Yes | high | 104.3 | 106.7 | 2.3 | 6.2x | 1950 |
| 4 | `arora2013` | 1 | Yes | high | 158.7 | 161.3 | 6.2 | 4.1x | 1298 |
| 5 | `negrete2024` | 1 | Yes | high | 197.8 | 200.3 | 2.4 | 3.3x | 1018 |
| 6 | `der2011` | 1 | Yes | high | 273.8 | 277.8 | 10.3 | 2.3x | 747 |
| 7 | `delatorre2018` | 1 | Yes | high | 643.2 | 646.1 | 10.8 | 1.0x | 320 |
| 1 | `mcelreath2025` | 1 | Yes | low | 47.5 | 48.4 | 0.6 | 13.0x | 4269 |
| 2 | `izzo2015` | 1 | Yes | low | 61.9 | 63.0 | 0.9 | 10.0x | 3309 |
| 3 | `gooding1990` | 1 | Yes | low | 104.1 | 106.4 | 2.2 | 5.9x | 1940 |
| 4 | `arora2013` | 1 | Yes | low | 160.5 | 163.9 | 5.1 | 3.9x | 1273 |
| 5 | `negrete2024` | 1 | Yes | low | 197.8 | 200.4 | 3.0 | 3.1x | 1019 |
| 6 | `der2011` | 1 | Yes | low | 251.6 | 271.4 | 14.3 | 2.5x | 812 |
| 7 | `delatorre2018` | 1 | Yes | low | 618.8 | 621.7 | 11.4 | 1.0x | 331 |
| 1 | `mcelreath2025` | 1 | No | high | 47.0 | 47.8 | 0.7 | 13.4x | 4291 |
| 2 | `izzo2015` | 1 | No | high | 61.9 | 63.1 | 1.0 | 10.2x | 3282 |
| 3 | `gooding1990` | 1 | No | high | 104.0 | 106.4 | 2.2 | 6.1x | 1960 |
| 4 | `arora2013` | 1 | No | high | 146.6 | 149.5 | 3.8 | 4.3x | 1398 |
| 5 | `negrete2024` | 1 | No | high | 195.3 | 197.6 | 2.9 | 3.2x | 1031 |
| 6 | `der2011` | 1 | No | high | 310.6 | 313.7 | 10.1 | 2.0x | 658 |
| 7 | `delatorre2018` | 1 | No | high | 630.6 | 632.3 | 9.9 | 1.0x | 324 |
| 1 | `mcelreath2025` | 1 | No | low | 47.7 | 48.6 | 0.6 | 12.9x | 4287 |
| 2 | `izzo2015` | 1 | No | low | 61.9 | 62.9 | 1.0 | 10.0x | 3284 |
| 3 | `gooding1990` | 1 | No | low | 103.4 | 105.7 | 2.3 | 6.0x | 1969 |
| 4 | `arora2013` | 1 | No | low | 149.3 | 152.0 | 3.8 | 4.1x | 1375 |
| 5 | `negrete2024` | 1 | No | low | 195.4 | 197.7 | 2.6 | 3.2x | 1030 |
| 6 | `der2011` | 1 | No | low | 290.8 | 294.4 | 10.4 | 2.1x | 706 |
| 7 | `delatorre2018` | 1 | No | low | 616.2 | 617.8 | 9.5 | 1.0x | 331 |

### Near-tangent robustness cases

#### Near tangent

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 0 | Yes | high | 48.5 | 49.6 | 0.6 | 6.6x | 4241 |
| 2 | `arora2013` | 0 | Yes | high | 54.1 | 55.0 | 0.8 | 5.9x | 3784 |
| 3 | `izzo2015` | 0 | Yes | high | 61.4 | 62.5 | 1.0 | 5.2x | 3334 |
| 4 | `gooding1990` | 0 | Yes | high | 85.2 | 87.2 | 1.8 | 3.8x | 2383 |
| 5 | `negrete2024` | 0 | Yes | high | 129.3 | 131.0 | 1.3 | 2.5x | 1558 |
| 6 | `der2011` | 0 | Yes | high | 252.3 | 255.1 | 10.1 | 1.3x | 816 |
| 7 | `delatorre2018` | 0 | Yes | high | 321.1 | 324.2 | 9.2 | 1.0x | 634 |

### Zero-revolution general cases

#### Curtiss book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 33.4 | 34.1 | 0.4 | 975.5x | 6123 |
| 2 | `mcelreath2025` | 0 | Yes | low | 46.9 | 47.8 | 0.6 | 694.7x | 4337 |
| 3 | `arora2013` | 0 | Yes | low | 56.5 | 57.4 | 0.8 | 576.7x | 3627 |
| 4 | `izzo2015` | 0 | Yes | low | 60.0 | 61.0 | 1.0 | 543.0x | 3422 |
| 5 | `gooding1990` | 0 | Yes | low | 83.6 | 85.2 | 1.5 | 389.7x | 2422 |
| 6 | `thorne2004` | 0 | Yes | low | 97.6 | 99.7 | 2.2 | 333.8x | 2086 |
| 7 | `der2011` | 0 | Yes | low | 101.1 | 103.2 | 2.3 | 322.3x | 2016 |
| 8 | `jiang2016` | 0 | Yes | low | 119.5 | 121.8 | 2.7 | 272.7x | 1698 |
| 9 | `avanzini2008` | 0 | Yes | low | 120.0 | 123.5 | 3.3 | 271.5x | 1702 |
| 10 | `negrete2024` | 0 | Yes | low | 127.4 | 128.9 | 1.2 | 255.7x | 1587 |
| 11 | `delatorre2018` | 0 | Yes | low | 322.3 | 325.8 | 9.7 | 101.1x | 636 |
| 12 | `pan2016` | 0 | Yes | low | 32582.3 | 32569.6 | 130.7 | 1.0x | 7 |

#### Der article I

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 36.4 | 37.1 | 0.4 | 854.6x | 5602 |
| 2 | `mcelreath2025` | 0 | Yes | low | 46.7 | 47.9 | 0.6 | 666.1x | 4373 |
| 3 | `arora2013` | 0 | Yes | low | 53.7 | 54.7 | 0.7 | 579.3x | 3789 |
| 4 | `izzo2015` | 0 | Yes | low | 60.1 | 61.2 | 1.0 | 517.6x | 3415 |
| 5 | `gooding1990` | 0 | Yes | low | 82.6 | 84.2 | 1.5 | 376.6x | 2474 |
| 6 | `thorne2004` | 0 | Yes | low | 106.4 | 108.4 | 2.4 | 292.4x | 1929 |
| 7 | `avanzini2008` | 0 | Yes | low | 126.0 | 129.1 | 3.5 | 246.9x | 1621 |
| 8 | `negrete2024` | 0 | Yes | low | 129.0 | 130.8 | 1.4 | 241.1x | 1563 |
| 9 | `der2011` | 0 | Yes | low | 137.5 | 139.5 | 3.6 | 226.2x | 1488 |
| 10 | `jiang2016` | 0 | Yes | low | 197.2 | 200.2 | 6.6 | 157.7x | 1038 |
| 11 | `delatorre2018` | 0 | Yes | low | 284.9 | 288.4 | 9.5 | 109.2x | 714 |
| 12 | `pan2016` | 0 | Yes | low | 31106.2 | 31144.7 | 154.0 | 1.0x | 7 |
| 1 | `battin1984` | 0 | No | high | 35.3 | 35.9 | 0.4 | 877.6x | 5759 |
| 2 | `mcelreath2025` | 0 | No | high | 46.8 | 47.8 | 0.6 | 662.0x | 4358 |
| 3 | `arora2013` | 0 | No | high | 56.1 | 57.1 | 0.7 | 552.2x | 3616 |
| 4 | `izzo2015` | 0 | No | high | 60.1 | 61.2 | 1.0 | 515.5x | 3374 |
| 5 | `gooding1990` | 0 | No | high | 82.9 | 84.5 | 1.5 | 373.7x | 2433 |
| 6 | `thorne2004` | 0 | No | high | 100.5 | 102.3 | 2.2 | 308.3x | 2038 |
| 7 | `avanzini2008` | 0 | No | high | 122.9 | 125.7 | 3.1 | 252.1x | 1656 |
| 8 | `negrete2024` | 0 | No | high | 128.2 | 130.0 | 1.4 | 241.7x | 1570 |
| 9 | `der2011` | 0 | No | high | 176.1 | 178.9 | 5.8 | 175.9x | 1161 |
| 10 | `jiang2016` | 0 | No | high | 248.1 | 250.9 | 8.9 | 124.9x | 826 |
| 11 | `delatorre2018` | 0 | No | high | 316.6 | 319.4 | 9.1 | 97.9x | 637 |
| 12 | `pan2016` | 0 | No | high | 30979.6 | 30990.9 | 106.4 | 1.0x | 7 |

#### Der article II

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | high | 36.9 | 37.5 | 0.4 | 969.9x | 5490 |
| 2 | `mcelreath2025` | 0 | Yes | high | 47.1 | 48.1 | 0.6 | 759.8x | 4355 |
| 3 | `arora2013` | 0 | Yes | high | 53.1 | 54.1 | 0.7 | 674.0x | 3804 |
| 4 | `izzo2015` | 0 | Yes | high | 60.4 | 61.3 | 1.0 | 592.5x | 3371 |
| 5 | `gooding1990` | 0 | Yes | high | 84.0 | 85.6 | 1.4 | 426.0x | 2415 |
| 6 | `thorne2004` | 0 | Yes | high | 88.3 | 90.1 | 2.7 | 405.3x | 2316 |
| 7 | `jiang2016` | 0 | Yes | high | 107.6 | 109.3 | 2.1 | 332.6x | 1889 |
| 8 | `der2011` | 0 | Yes | high | 118.5 | 120.6 | 2.8 | 302.0x | 1721 |
| 9 | `avanzini2008` | 0 | Yes | high | 126.4 | 129.8 | 3.5 | 283.1x | 1610 |
| 10 | `negrete2024` | 0 | Yes | high | 129.7 | 133.0 | 1.9 | 275.9x | 1556 |
| 11 | `delatorre2018` | 0 | Yes | high | 313.1 | 320.7 | 9.6 | 114.3x | 649 |
| 12 | `pan2016` | 0 | Yes | high | 35787.9 | 35899.0 | 176.0 | 1.0x | 6 |
| 1 | `battin1984` | 0 | No | high | 35.4 | 35.9 | 0.4 | 1024.9x | 5719 |
| 2 | `mcelreath2025` | 0 | No | high | 46.7 | 47.7 | 0.6 | 776.9x | 4350 |
| 3 | `izzo2015` | 0 | No | high | 60.6 | 61.6 | 0.9 | 598.7x | 3356 |
| 4 | `arora2013` | 0 | No | high | 63.3 | 64.3 | 0.8 | 573.2x | 3207 |
| 5 | `gooding1990` | 0 | No | high | 85.1 | 86.6 | 1.4 | 426.3x | 2370 |
| 6 | `thorne2004` | 0 | No | high | 89.8 | 92.0 | 1.9 | 404.0x | 2286 |
| 7 | `jiang2016` | 0 | No | high | 106.9 | 108.8 | 2.1 | 339.4x | 1909 |
| 8 | `avanzini2008` | 0 | No | high | 120.8 | 123.6 | 3.1 | 300.3x | 1682 |
| 9 | `negrete2024` | 0 | No | high | 129.5 | 131.3 | 1.4 | 280.2x | 1558 |
| 10 | `der2011` | 0 | No | high | 138.2 | 140.3 | 3.3 | 262.5x | 1469 |
| 11 | `delatorre2018` | 0 | No | high | 299.4 | 302.1 | 9.1 | 121.2x | 677 |
| 12 | `pan2016` | 0 | No | high | 36280.6 | 36334.5 | 206.3 | 1.0x | 6 |

### Zero-revolution hyperbolic cases

#### GMAT hyperbolic

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 34.8 | 35.3 | 0.4 | 1066.1x | 5850 |
| 2 | `mcelreath2025` | 0 | Yes | low | 47.9 | 49.0 | 0.7 | 774.6x | 4233 |
| 3 | `arora2013` | 0 | Yes | low | 52.7 | 53.6 | 0.7 | 704.0x | 3886 |
| 4 | `izzo2015` | 0 | Yes | low | 61.9 | 62.9 | 1.1 | 599.4x | 3312 |
| 5 | `gooding1990` | 0 | Yes | low | 85.2 | 86.6 | 1.8 | 435.5x | 2366 |
| 6 | `thorne2004` | 0 | Yes | low | 104.4 | 106.6 | 2.5 | 355.4x | 1990 |
| 7 | `avanzini2008` | 0 | Yes | low | 126.4 | 129.3 | 3.6 | 293.5x | 1619 |
| 8 | `jiang2016` | 0 | Yes | low | 131.6 | 133.9 | 3.0 | 281.9x | 1544 |
| 9 | `negrete2024` | 0 | Yes | low | 133.5 | 135.3 | 1.3 | 277.9x | 1513 |
| 10 | `der2011` | 0 | Yes | low | 141.8 | 143.7 | 3.2 | 261.6x | 1434 |
| 11 | `delatorre2018` | 0 | Yes | low | 334.4 | 336.8 | 9.4 | 111.0x | 616 |
| 12 | `pan2016` | 0 | Yes | low | 37101.9 | 37101.5 | 229.7 | 1.0x | 6 |
| 1 | `battin1984` | 0 | No | low | 35.2 | 35.8 | 0.4 | 1052.4x | 5743 |
| 2 | `mcelreath2025` | 0 | No | low | 47.6 | 51.5 | 1.0 | 778.2x | 4238 |
| 3 | `arora2013` | 0 | No | low | 62.0 | 62.8 | 0.8 | 597.5x | 3277 |
| 4 | `izzo2015` | 0 | No | low | 63.0 | 64.0 | 1.0 | 588.0x | 3233 |
| 5 | `gooding1990` | 0 | No | low | 84.3 | 86.1 | 1.4 | 439.4x | 2402 |
| 6 | `thorne2004` | 0 | No | low | 106.7 | 108.5 | 2.4 | 347.2x | 1953 |
| 7 | `der2011` | 0 | No | low | 115.6 | 117.9 | 3.1 | 320.5x | 1785 |
| 8 | `avanzini2008` | 0 | No | low | 128.0 | 131.4 | 3.4 | 289.4x | 1593 |
| 9 | `jiang2016` | 0 | No | low | 132.5 | 134.6 | 2.8 | 279.6x | 1541 |
| 10 | `negrete2024` | 0 | No | low | 133.4 | 135.1 | 1.4 | 277.7x | 1509 |
| 11 | `delatorre2018` | 0 | No | low | 329.2 | 333.4 | 9.7 | 112.5x | 617 |
| 12 | `pan2016` | 0 | No | low | 37044.4 | 37066.7 | 413.7 | 1.0x | 6 |

### Zero-revolution nominal cases

#### Battin book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `gauss1809` | 0 | Yes | low | 33.3 | 34.5 | 0.6 | 1095.4x | 6074 |
| 2 | `battin1984` | 0 | Yes | low | 36.6 | 37.8 | 0.5 | 996.6x | 5600 |
| 3 | `mcelreath2025` | 0 | Yes | low | 61.2 | 63.7 | 1.2 | 596.0x | 3344 |
| 4 | `arora2013` | 0 | Yes | low | 78.1 | 80.7 | 1.8 | 467.1x | 2661 |
| 5 | `izzo2015` | 0 | Yes | low | 83.7 | 86.9 | 1.9 | 435.8x | 2492 |
| 6 | `jiang2016` | 0 | Yes | low | 105.3 | 109.1 | 3.8 | 346.4x | 1960 |
| 7 | `gooding1990` | 0 | Yes | low | 109.1 | 113.2 | 3.8 | 334.3x | 1854 |
| 8 | `thorne2004` | 0 | Yes | low | 124.1 | 129.1 | 6.4 | 293.9x | 1644 |
| 9 | `negrete2024` | 0 | Yes | low | 129.1 | 133.2 | 5.2 | 282.5x | 1567 |
| 10 | `vallado2013` | 0 | Yes | low | 131.4 | 134.8 | 4.6 | 277.6x | 1531 |
| 11 | `der2011` | 0 | Yes | low | 156.2 | 161.5 | 8.3 | 233.5x | 1317 |
| 12 | `avanzini2008` | 0 | Yes | low | 156.9 | 162.5 | 6.5 | 232.5x | 1300 |
| 13 | `delatorre2018` | 0 | Yes | low | 310.9 | 315.6 | 15.4 | 117.3x | 658 |
| 14 | `pan2016` | 0 | Yes | low | 36477.1 | 36498.6 | 196.2 | 1.0x | 6 |

#### Vallado book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 38.2 | 39.5 | 0.5 | 10.5x | 5341 |
| 2 | `gauss1809` | 0 | Yes | low | 51.9 | 53.5 | 0.8 | 7.7x | 3915 |
| 3 | `mcelreath2025` | 0 | Yes | low | 62.1 | 66.3 | 1.7 | 6.5x | 3330 |
| 4 | `vallado2013` | 0 | Yes | low | 67.1 | 69.3 | 0.9 | 6.0x | 2988 |
| 5 | `arora2013` | 0 | Yes | low | 78.3 | 80.9 | 1.9 | 5.1x | 2610 |
| 6 | `izzo2015` | 0 | Yes | low | 84.0 | 87.2 | 2.1 | 4.8x | 2493 |
| 7 | `gooding1990` | 0 | Yes | low | 109.8 | 113.4 | 3.6 | 3.7x | 1854 |
| 8 | `pan2016` | 0 | Yes | low | 115.2 | 119.6 | 3.7 | 3.5x | 1777 |
| 9 | `thorne2004` | 0 | Yes | low | 122.8 | 127.2 | 5.0 | 3.3x | 1675 |
| 10 | `negrete2024` | 0 | Yes | low | 130.5 | 134.0 | 5.2 | 3.1x | 1545 |
| 11 | `avanzini2008` | 0 | Yes | low | 167.0 | 172.9 | 7.2 | 2.4x | 1236 |
| 12 | `der2011` | 0 | Yes | low | 193.4 | 198.6 | 10.6 | 2.1x | 1061 |
| 13 | `delatorre2018` | 0 | Yes | low | 356.2 | 360.9 | 15.9 | 1.1x | 583 |
| 14 | `jiang2016` | 0 | Yes | low | 401.0 | 406.7 | 16.2 | 1.0x | 516 |

<!-- performance-comparison:end -->
