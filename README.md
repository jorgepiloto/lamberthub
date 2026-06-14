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

- Generated: 2026-06-14 15:44 UTC
- Commit: `a3af78b1af9a`
- Environment: Linux, Python 3.12.13, AMD EPYC 7763 64-Core Processor

Times are in microseconds (lower is better).
**Speedup** is relative to the slowest solver for each case.
All solvers are JIT-warmed before timing begins.

### Near-minimum-energy multi-revolution cases

#### Near-minimum-energy

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `der2011` | 1 | Yes | high | 60.6 | 63.0 | 1.1 | 11.9x | 3386 |
| 2 | `mcelreath2025` | 1 | Yes | high | 61.2 | 63.3 | 1.2 | 11.7x | 3366 |
| 3 | `izzo2015` | 1 | Yes | high | 84.1 | 86.7 | 1.7 | 8.5x | 2436 |
| 4 | `gooding1990` | 1 | Yes | high | 129.1 | 134.0 | 3.9 | 5.6x | 1499 |
| 5 | `negrete2024` | 1 | Yes | high | 197.1 | 201.5 | 5.3 | 3.6x | 1022 |
| 6 | `arora2013` | 1 | Yes | high | 251.5 | 257.1 | 20.3 | 2.9x | 821 |
| 7 | `delatorre2018` | 1 | Yes | high | 718.6 | 717.5 | 12.0 | 1.0x | 287 |
| 1 | `mcelreath2025` | 1 | Yes | low | 62.1 | 64.5 | 1.2 | 12.4x | 3338 |
| 2 | `izzo2015` | 1 | Yes | low | 85.4 | 88.2 | 1.7 | 9.0x | 2415 |
| 3 | `der2011` | 1 | Yes | low | 86.3 | 89.2 | 1.6 | 8.9x | 2362 |
| 4 | `gooding1990` | 1 | Yes | low | 126.8 | 131.8 | 3.7 | 6.0x | 1607 |
| 5 | `negrete2024` | 1 | Yes | low | 197.1 | 201.9 | 5.5 | 3.9x | 1020 |
| 6 | `arora2013` | 1 | Yes | low | 274.0 | 279.0 | 19.9 | 2.8x | 748 |
| 7 | `delatorre2018` | 1 | Yes | low | 766.9 | 768.9 | 6.4 | 1.0x | 263 |
| 1 | `der2011` | 2 | Yes | high | 55.9 | 57.9 | 1.0 | 12.9x | 3650 |
| 2 | `mcelreath2025` | 2 | Yes | high | 61.5 | 63.9 | 1.3 | 11.7x | 3335 |
| 3 | `izzo2015` | 2 | Yes | high | 84.4 | 87.1 | 1.7 | 8.5x | 2455 |
| 4 | `gooding1990` | 2 | Yes | high | 128.7 | 133.5 | 4.0 | 5.6x | 1598 |
| 5 | `negrete2024` | 2 | Yes | high | 197.4 | 202.3 | 6.2 | 3.6x | 1021 |
| 6 | `arora2013` | 2 | Yes | high | 252.5 | 258.8 | 20.7 | 2.9x | 817 |
| 7 | `delatorre2018` | 2 | Yes | high | 720.2 | 721.7 | 13.2 | 1.0x | 285 |
| 1 | `der2011` | 2 | Yes | low | 55.9 | 57.8 | 0.9 | 13.4x | 3642 |
| 2 | `mcelreath2025` | 2 | Yes | low | 61.5 | 63.7 | 1.3 | 12.2x | 3367 |
| 3 | `izzo2015` | 2 | Yes | low | 84.1 | 87.2 | 2.0 | 8.9x | 2446 |
| 4 | `gooding1990` | 2 | Yes | low | 129.2 | 133.9 | 3.5 | 5.8x | 1588 |
| 5 | `negrete2024` | 2 | Yes | low | 197.6 | 202.4 | 6.3 | 3.8x | 1020 |
| 6 | `arora2013` | 2 | Yes | low | 269.7 | 276.2 | 20.9 | 2.8x | 767 |
| 7 | `delatorre2018` | 2 | Yes | low | 750.3 | 753.1 | 9.5 | 1.0x | 272 |

### One-revolution branch cases

#### Der article II

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 1 | No | high | 48.3 | 49.2 | 0.6 | 13.2x | 4245 |
| 2 | `der2011` | 1 | No | high | 50.0 | 51.0 | 0.8 | 12.8x | 4086 |
| 3 | `izzo2015` | 1 | No | high | 60.8 | 61.9 | 1.0 | 10.5x | 3338 |
| 4 | `gooding1990` | 1 | No | high | 94.1 | 96.2 | 1.9 | 6.8x | 2162 |
| 5 | `arora2013` | 1 | No | high | 148.6 | 151.3 | 3.9 | 4.3x | 1396 |
| 6 | `negrete2024` | 1 | No | high | 195.4 | 197.9 | 3.5 | 3.3x | 1032 |
| 7 | `delatorre2018` | 1 | No | high | 639.2 | 645.1 | 11.6 | 1.0x | 317 |
| 1 | `mcelreath2025` | 1 | No | low | 48.3 | 49.3 | 0.7 | 13.0x | 4201 |
| 2 | `der2011` | 1 | No | low | 49.5 | 50.4 | 0.8 | 12.6x | 4141 |
| 3 | `izzo2015` | 1 | No | low | 60.7 | 62.1 | 1.1 | 10.3x | 3378 |
| 4 | `gooding1990` | 1 | No | low | 94.6 | 97.0 | 2.1 | 6.6x | 2158 |
| 5 | `arora2013` | 1 | No | low | 150.8 | 153.6 | 3.9 | 4.1x | 1364 |
| 6 | `negrete2024` | 1 | No | low | 195.4 | 197.9 | 3.1 | 3.2x | 1031 |
| 7 | `delatorre2018` | 1 | No | low | 625.3 | 625.2 | 10.9 | 1.0x | 325 |
| 1 | `mcelreath2025` | 1 | Yes | high | 48.0 | 48.9 | 0.6 | 13.6x | 4207 |
| 2 | `der2011` | 1 | Yes | high | 49.5 | 50.6 | 1.0 | 13.2x | 4171 |
| 3 | `izzo2015` | 1 | Yes | high | 60.0 | 61.2 | 1.0 | 10.9x | 3425 |
| 4 | `gooding1990` | 1 | Yes | high | 97.2 | 100.2 | 2.5 | 6.7x | 2108 |
| 5 | `arora2013` | 1 | Yes | high | 160.3 | 163.5 | 5.0 | 4.1x | 1282 |
| 6 | `negrete2024` | 1 | Yes | high | 198.3 | 201.1 | 4.0 | 3.3x | 1018 |
| 7 | `delatorre2018` | 1 | Yes | high | 652.4 | 654.4 | 10.7 | 1.0x | 313 |
| 1 | `mcelreath2025` | 1 | Yes | low | 47.5 | 48.5 | 0.6 | 13.3x | 4235 |
| 2 | `der2011` | 1 | Yes | low | 49.3 | 50.2 | 0.8 | 12.8x | 4144 |
| 3 | `izzo2015` | 1 | Yes | low | 61.5 | 62.6 | 1.0 | 10.3x | 3343 |
| 4 | `gooding1990` | 1 | Yes | low | 94.9 | 97.3 | 2.1 | 6.7x | 2143 |
| 5 | `arora2013` | 1 | Yes | low | 163.8 | 166.7 | 4.6 | 3.9x | 1263 |
| 6 | `negrete2024` | 1 | Yes | low | 198.2 | 200.5 | 3.0 | 3.2x | 1018 |
| 7 | `delatorre2018` | 1 | Yes | low | 632.3 | 633.2 | 10.9 | 1.0x | 325 |

### Near-tangent robustness cases

#### Near tangent

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 0 | Yes | high | 48.3 | 49.2 | 0.6 | 6.8x | 4217 |
| 2 | `arora2013` | 0 | Yes | high | 54.8 | 55.6 | 0.6 | 6.0x | 3690 |
| 3 | `izzo2015` | 0 | Yes | high | 60.2 | 61.5 | 1.0 | 5.4x | 3388 |
| 4 | `der2011` | 0 | Yes | high | 66.4 | 67.7 | 0.9 | 4.9x | 3081 |
| 5 | `gooding1990` | 0 | Yes | high | 81.9 | 83.6 | 1.3 | 4.0x | 2504 |
| 6 | `negrete2024` | 0 | Yes | high | 129.9 | 131.5 | 1.3 | 2.5x | 1552 |
| 7 | `delatorre2018` | 0 | Yes | high | 327.3 | 330.4 | 8.8 | 1.0x | 616 |

### Zero-revolution general cases

#### Curtiss book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 25.9 | 26.3 | 0.3 | 991.5x | 7847 |
| 2 | `der2011` | 0 | Yes | low | 30.3 | 30.7 | 0.3 | 849.1x | 6761 |
| 3 | `jiang2016` | 0 | Yes | low | 31.7 | 32.2 | 0.3 | 810.0x | 6463 |
| 4 | `mcelreath2025` | 0 | Yes | low | 36.9 | 37.4 | 0.4 | 696.9x | 5539 |
| 5 | `arora2013` | 0 | Yes | low | 43.5 | 45.7 | 0.5 | 591.0x | 4647 |
| 6 | `izzo2015` | 0 | Yes | low | 46.1 | 46.8 | 0.7 | 557.8x | 4447 |
| 7 | `gooding1990` | 0 | Yes | low | 62.1 | 63.2 | 0.9 | 413.8x | 3250 |
| 8 | `thorne2004` | 0 | Yes | low | 76.2 | 77.7 | 1.2 | 337.4x | 2669 |
| 9 | `avanzini2008` | 0 | Yes | low | 93.9 | 96.3 | 1.9 | 273.7x | 2148 |
| 10 | `negrete2024` | 0 | Yes | low | 105.2 | 106.4 | 0.8 | 244.2x | 1914 |
| 11 | `delatorre2018` | 0 | Yes | low | 245.2 | 251.4 | 7.1 | 104.8x | 802 |
| 12 | `pan2016` | 0 | Yes | low | 25698.4 | 25710.4 | 81.5 | 1.0x | 8 |

#### Der article I

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | high | 27.2 | 27.6 | 0.2 | 906.2x | 7433 |
| 2 | `jiang2016` | 0 | No | high | 33.2 | 33.7 | 0.3 | 740.9x | 6160 |
| 3 | `mcelreath2025` | 0 | No | high | 36.7 | 37.2 | 0.4 | 671.3x | 5519 |
| 4 | `der2011` | 0 | No | high | 43.6 | 44.5 | 0.8 | 564.5x | 4719 |
| 5 | `arora2013` | 0 | No | high | 44.1 | 45.1 | 0.5 | 557.8x | 4605 |
| 6 | `izzo2015` | 0 | No | high | 46.8 | 47.4 | 0.7 | 526.6x | 4359 |
| 7 | `gooding1990` | 0 | No | high | 63.0 | 64.1 | 0.8 | 390.9x | 3225 |
| 8 | `thorne2004` | 0 | No | high | 79.4 | 80.8 | 1.5 | 309.9x | 2598 |
| 9 | `avanzini2008` | 0 | No | high | 98.2 | 100.4 | 1.9 | 250.7x | 2071 |
| 10 | `negrete2024` | 0 | No | high | 106.2 | 107.3 | 0.8 | 231.7x | 1895 |
| 11 | `delatorre2018` | 0 | No | high | 246.2 | 248.2 | 6.6 | 100.0x | 794 |
| 12 | `pan2016` | 0 | No | high | 24621.3 | 24631.7 | 54.0 | 1.0x | 9 |
| 1 | `battin1984` | 0 | Yes | low | 28.3 | 28.7 | 0.3 | 864.3x | 7171 |
| 2 | `der2011` | 0 | Yes | low | 30.7 | 31.2 | 0.4 | 795.0x | 6650 |
| 3 | `jiang2016` | 0 | Yes | low | 32.2 | 32.7 | 0.4 | 758.7x | 6235 |
| 4 | `mcelreath2025` | 0 | Yes | low | 36.3 | 37.3 | 0.4 | 673.5x | 5615 |
| 5 | `arora2013` | 0 | Yes | low | 41.8 | 42.4 | 0.4 | 584.7x | 4855 |
| 6 | `izzo2015` | 0 | Yes | low | 46.0 | 46.6 | 0.7 | 531.6x | 4474 |
| 7 | `gooding1990` | 0 | Yes | low | 61.8 | 63.1 | 0.9 | 395.4x | 3275 |
| 8 | `thorne2004` | 0 | Yes | low | 81.8 | 83.3 | 1.5 | 298.8x | 2509 |
| 9 | `avanzini2008` | 0 | Yes | low | 98.8 | 100.8 | 1.9 | 247.5x | 2053 |
| 10 | `negrete2024` | 0 | Yes | low | 105.9 | 107.1 | 0.9 | 230.8x | 1902 |
| 11 | `delatorre2018` | 0 | Yes | low | 221.1 | 223.2 | 5.5 | 110.5x | 889 |
| 12 | `pan2016` | 0 | Yes | low | 24443.2 | 24452.4 | 74.6 | 1.0x | 9 |

#### Der article II

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | high | 27.5 | 27.8 | 0.3 | 1016.1x | 7447 |
| 2 | `der2011` | 0 | No | high | 30.5 | 31.0 | 0.3 | 914.4x | 6655 |
| 3 | `jiang2016` | 0 | No | high | 31.4 | 31.9 | 0.3 | 888.8x | 6470 |
| 4 | `mcelreath2025` | 0 | No | high | 36.5 | 37.0 | 0.4 | 766.0x | 5572 |
| 5 | `izzo2015` | 0 | No | high | 46.4 | 47.1 | 0.6 | 601.7x | 4428 |
| 6 | `arora2013` | 0 | No | high | 49.2 | 54.0 | 1.5 | 567.8x | 4137 |
| 7 | `gooding1990` | 0 | No | high | 63.2 | 64.5 | 0.8 | 441.9x | 3205 |
| 8 | `thorne2004` | 0 | No | high | 70.3 | 71.6 | 1.4 | 397.3x | 2945 |
| 9 | `avanzini2008` | 0 | No | high | 95.4 | 97.5 | 1.7 | 292.6x | 2114 |
| 10 | `negrete2024` | 0 | No | high | 106.5 | 107.6 | 0.9 | 262.2x | 1896 |
| 11 | `delatorre2018` | 0 | No | high | 228.5 | 230.6 | 5.9 | 122.2x | 858 |
| 12 | `pan2016` | 0 | No | high | 27923.2 | 27943.7 | 86.3 | 1.0x | 8 |
| 1 | `battin1984` | 0 | Yes | high | 28.7 | 29.1 | 0.3 | 987.4x | 7040 |
| 2 | `der2011` | 0 | Yes | high | 30.2 | 30.9 | 0.3 | 937.4x | 6684 |
| 3 | `jiang2016` | 0 | Yes | high | 31.3 | 31.7 | 0.3 | 907.1x | 6405 |
| 4 | `mcelreath2025` | 0 | Yes | high | 36.5 | 37.3 | 0.5 | 776.2x | 5519 |
| 5 | `arora2013` | 0 | Yes | high | 41.3 | 41.9 | 0.4 | 686.0x | 4863 |
| 6 | `izzo2015` | 0 | Yes | high | 46.9 | 47.5 | 0.7 | 604.4x | 4402 |
| 7 | `gooding1990` | 0 | Yes | high | 63.2 | 64.3 | 0.8 | 448.5x | 3205 |
| 8 | `thorne2004` | 0 | Yes | high | 69.1 | 71.2 | 1.4 | 410.3x | 2961 |
| 9 | `avanzini2008` | 0 | Yes | high | 99.7 | 101.7 | 2.1 | 284.5x | 2034 |
| 10 | `negrete2024` | 0 | Yes | high | 106.8 | 107.9 | 0.9 | 265.4x | 1889 |
| 11 | `delatorre2018` | 0 | Yes | high | 242.8 | 244.9 | 6.2 | 116.8x | 818 |
| 12 | `pan2016` | 0 | Yes | high | 28351.7 | 28440.6 | 512.1 | 1.0x | 8 |

### Zero-revolution hyperbolic cases

#### GMAT hyperbolic

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | low | 35.1 | 35.6 | 0.3 | 1046.8x | 5797 |
| 2 | `der2011` | 0 | No | low | 39.4 | 40.1 | 0.5 | 931.8x | 5190 |
| 3 | `jiang2016` | 0 | No | low | 39.6 | 40.3 | 0.5 | 927.3x | 5112 |
| 4 | `mcelreath2025` | 0 | No | low | 48.3 | 49.3 | 0.6 | 759.7x | 4223 |
| 5 | `izzo2015` | 0 | No | low | 60.1 | 61.1 | 0.9 | 610.7x | 3392 |
| 6 | `arora2013` | 0 | No | low | 61.8 | 62.9 | 0.7 | 594.0x | 3285 |
| 7 | `gooding1990` | 0 | No | low | 79.4 | 81.0 | 1.3 | 462.5x | 2546 |
| 8 | `thorne2004` | 0 | No | low | 101.5 | 104.3 | 3.0 | 361.5x | 2039 |
| 9 | `avanzini2008` | 0 | No | low | 125.7 | 129.1 | 3.9 | 291.9x | 1618 |
| 10 | `negrete2024` | 0 | No | low | 133.2 | 136.0 | 1.5 | 275.6x | 1514 |
| 11 | `delatorre2018` | 0 | No | low | 337.0 | 340.1 | 9.2 | 108.9x | 605 |
| 12 | `pan2016` | 0 | No | low | 36702.9 | 36703.9 | 83.8 | 1.0x | 6 |
| 1 | `battin1984` | 0 | Yes | low | 34.7 | 35.3 | 0.4 | 1068.9x | 5845 |
| 2 | `der2011` | 0 | Yes | low | 39.5 | 40.3 | 0.6 | 940.2x | 5173 |
| 3 | `jiang2016` | 0 | Yes | low | 39.6 | 40.5 | 0.6 | 937.3x | 5169 |
| 4 | `mcelreath2025` | 0 | Yes | low | 48.1 | 49.1 | 0.6 | 772.2x | 4228 |
| 5 | `arora2013` | 0 | Yes | low | 53.3 | 54.3 | 0.7 | 696.7x | 3812 |
| 6 | `izzo2015` | 0 | Yes | low | 59.6 | 60.8 | 1.0 | 622.7x | 3447 |
| 7 | `gooding1990` | 0 | Yes | low | 78.3 | 80.1 | 1.4 | 474.4x | 2559 |
| 8 | `thorne2004` | 0 | Yes | low | 100.0 | 110.2 | 12.1 | 371.2x | 2088 |
| 9 | `avanzini2008` | 0 | Yes | low | 125.5 | 129.4 | 4.5 | 295.9x | 1618 |
| 10 | `negrete2024` | 0 | Yes | low | 133.5 | 136.0 | 1.5 | 278.2x | 1511 |
| 11 | `delatorre2018` | 0 | Yes | low | 335.2 | 338.4 | 9.3 | 110.8x | 606 |
| 12 | `pan2016` | 0 | Yes | low | 37135.4 | 37196.9 | 175.9 | 1.0x | 6 |

### Zero-revolution nominal cases

#### Battin book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `gauss1809` | 0 | Yes | low | 30.0 | 30.5 | 0.5 | 1064.5x | 6832 |
| 2 | `battin1984` | 0 | Yes | low | 33.7 | 34.3 | 0.3 | 947.0x | 6048 |
| 3 | `jiang2016` | 0 | Yes | low | 44.1 | 44.7 | 0.8 | 723.6x | 4686 |
| 4 | `der2011` | 0 | Yes | low | 44.3 | 44.9 | 0.7 | 720.6x | 4645 |
| 5 | `mcelreath2025` | 0 | Yes | low | 51.8 | 52.4 | 0.8 | 615.7x | 3977 |
| 6 | `izzo2015` | 0 | Yes | low | 71.0 | 71.9 | 1.1 | 449.4x | 2904 |
| 7 | `arora2013` | 0 | Yes | low | 72.4 | 76.5 | 1.3 | 440.7x | 2860 |
| 8 | `gooding1990` | 0 | Yes | low | 91.6 | 92.8 | 1.2 | 348.5x | 2233 |
| 9 | `thorne2004` | 0 | Yes | low | 109.8 | 111.4 | 1.7 | 290.5x | 1878 |
| 10 | `vallado2013` | 0 | Yes | low | 120.2 | 121.4 | 1.5 | 265.5x | 1690 |
| 11 | `negrete2024` | 0 | Yes | low | 123.5 | 124.9 | 1.2 | 258.4x | 1628 |
| 12 | `avanzini2008` | 0 | Yes | low | 137.2 | 139.3 | 1.8 | 232.6x | 1486 |
| 13 | `delatorre2018` | 0 | Yes | low | 285.0 | 287.1 | 6.3 | 111.9x | 719 |
| 14 | `pan2016` | 0 | Yes | low | 31904.7 | 31922.4 | 109.1 | 1.0x | 7 |

#### Vallado book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 35.4 | 35.8 | 0.4 | 9.0x | 5767 |
| 2 | `gauss1809` | 0 | Yes | low | 48.9 | 49.0 | 1.6 | 6.5x | 4280 |
| 3 | `jiang2016` | 0 | Yes | low | 49.1 | 49.8 | 0.9 | 6.5x | 4215 |
| 4 | `mcelreath2025` | 0 | Yes | low | 51.9 | 52.5 | 0.8 | 6.2x | 3961 |
| 5 | `der2011` | 0 | Yes | low | 60.7 | 61.5 | 0.9 | 5.3x | 3398 |
| 6 | `vallado2013` | 0 | Yes | low | 62.1 | 62.8 | 0.8 | 5.1x | 3272 |
| 7 | `izzo2015` | 0 | Yes | low | 70.5 | 71.5 | 1.0 | 4.5x | 2924 |
| 8 | `arora2013` | 0 | Yes | low | 72.8 | 73.8 | 0.9 | 4.4x | 2836 |
| 9 | `gooding1990` | 0 | Yes | low | 92.5 | 93.6 | 1.2 | 3.5x | 2207 |
| 10 | `pan2016` | 0 | Yes | low | 101.6 | 103.0 | 1.6 | 3.1x | 2007 |
| 11 | `thorne2004` | 0 | Yes | low | 108.3 | 110.1 | 1.8 | 3.0x | 1890 |
| 12 | `negrete2024` | 0 | Yes | low | 124.1 | 125.4 | 1.0 | 2.6x | 1620 |
| 13 | `avanzini2008` | 0 | Yes | low | 146.1 | 155.1 | 10.5 | 2.2x | 1408 |
| 14 | `delatorre2018` | 0 | Yes | low | 319.6 | 321.4 | 6.7 | 1.0x | 641 |

<!-- performance-comparison:end -->
