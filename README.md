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

- Generated: 2026-08-03 18:53 UTC
- Commit: `20a8135dd7d8`
- Environment: Linux, Python 3.12.13, AMD EPYC 9V74 80-Core Processor

Times are in microseconds (lower is better).
**Speedup** is relative to the slowest solver for each case.
All solvers are JIT-warmed before timing begins.

### Near-minimum-energy multi-revolution cases

#### Near-minimum-energy

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 1 | Yes | high | 47.2 | 48.0 | 0.6 | 13.9x | 4323 |
| 2 | `der2011` | 1 | Yes | high | 53.4 | 54.3 | 0.8 | 12.3x | 3842 |
| 3 | `izzo2015` | 1 | Yes | high | 60.1 | 61.0 | 0.9 | 10.9x | 3431 |
| 4 | `gooding1990` | 1 | Yes | high | 94.6 | 96.7 | 1.8 | 6.9x | 2191 |
| 5 | `negrete2024` | 1 | Yes | high | 192.7 | 195.0 | 3.0 | 3.4x | 1046 |
| 6 | `arora2013` | 1 | Yes | high | 211.1 | 214.5 | 8.6 | 3.1x | 975 |
| 7 | `delatorre2018` | 1 | Yes | high | 655.1 | 658.2 | 9.4 | 1.0x | 317 |
| 1 | `mcelreath2025` | 1 | Yes | low | 47.2 | 48.0 | 0.6 | 15.0x | 4326 |
| 2 | `izzo2015` | 1 | Yes | low | 61.4 | 62.5 | 1.0 | 11.5x | 3361 |
| 3 | `der2011` | 1 | Yes | low | 79.2 | 80.6 | 1.3 | 8.9x | 2614 |
| 4 | `gooding1990` | 1 | Yes | low | 94.2 | 96.3 | 1.9 | 7.5x | 2164 |
| 5 | `negrete2024` | 1 | Yes | low | 192.4 | 194.9 | 3.0 | 3.7x | 1048 |
| 6 | `arora2013` | 1 | Yes | low | 231.9 | 235.2 | 9.0 | 3.1x | 878 |
| 7 | `delatorre2018` | 1 | Yes | low | 707.8 | 709.1 | 8.4 | 1.0x | 289 |
| 1 | `mcelreath2025` | 2 | Yes | high | 47.6 | 48.7 | 0.7 | 13.8x | 4263 |
| 2 | `der2011` | 2 | Yes | high | 49.2 | 50.0 | 0.8 | 13.4x | 4177 |
| 3 | `izzo2015` | 2 | Yes | high | 61.0 | 62.0 | 0.9 | 10.8x | 3368 |
| 4 | `gooding1990` | 2 | Yes | high | 97.0 | 99.0 | 1.7 | 6.8x | 2142 |
| 5 | `negrete2024` | 2 | Yes | high | 193.8 | 196.5 | 3.5 | 3.4x | 1038 |
| 6 | `arora2013` | 2 | Yes | high | 214.7 | 218.1 | 7.6 | 3.1x | 952 |
| 7 | `delatorre2018` | 2 | Yes | high | 657.6 | 658.6 | 9.4 | 1.0x | 315 |
| 1 | `mcelreath2025` | 2 | Yes | low | 47.2 | 48.1 | 0.6 | 14.5x | 4335 |
| 2 | `der2011` | 2 | Yes | low | 50.1 | 62.5 | 15.9 | 13.7x | 4163 |
| 3 | `izzo2015` | 2 | Yes | low | 60.3 | 61.2 | 0.9 | 11.4x | 3382 |
| 4 | `gooding1990` | 2 | Yes | low | 96.7 | 98.7 | 1.8 | 7.1x | 2108 |
| 5 | `negrete2024` | 2 | Yes | low | 194.1 | 196.2 | 2.5 | 3.5x | 1038 |
| 6 | `arora2013` | 2 | Yes | low | 226.7 | 229.3 | 6.9 | 3.0x | 902 |
| 7 | `delatorre2018` | 2 | Yes | low | 686.3 | 687.1 | 8.2 | 1.0x | 300 |

### One-revolution branch cases

#### Der article II

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 1 | No | high | 47.7 | 48.8 | 0.7 | 13.2x | 4238 |
| 2 | `der2011` | 1 | No | high | 49.8 | 50.7 | 0.7 | 12.7x | 4089 |
| 3 | `izzo2015` | 1 | No | high | 61.6 | 63.9 | 1.1 | 10.2x | 3332 |
| 4 | `gooding1990` | 1 | No | high | 97.1 | 99.8 | 2.2 | 6.5x | 2110 |
| 5 | `arora2013` | 1 | No | high | 147.5 | 150.3 | 4.1 | 4.3x | 1393 |
| 6 | `negrete2024` | 1 | No | high | 194.3 | 197.3 | 3.9 | 3.2x | 1039 |
| 7 | `delatorre2018` | 1 | No | high | 630.6 | 633.5 | 10.7 | 1.0x | 324 |
| 1 | `mcelreath2025` | 1 | No | low | 48.0 | 49.1 | 0.7 | 12.9x | 4248 |
| 2 | `der2011` | 1 | No | low | 50.1 | 51.0 | 0.7 | 12.3x | 4097 |
| 3 | `izzo2015` | 1 | No | low | 62.1 | 63.3 | 1.1 | 10.0x | 3305 |
| 4 | `gooding1990` | 1 | No | low | 97.0 | 99.5 | 2.3 | 6.4x | 2132 |
| 5 | `arora2013` | 1 | No | low | 149.3 | 152.1 | 4.5 | 4.1x | 1372 |
| 6 | `negrete2024` | 1 | No | low | 194.5 | 197.3 | 2.7 | 3.2x | 1036 |
| 7 | `delatorre2018` | 1 | No | low | 618.5 | 620.4 | 10.2 | 1.0x | 331 |
| 1 | `mcelreath2025` | 1 | Yes | high | 47.9 | 49.0 | 0.7 | 13.5x | 4249 |
| 2 | `der2011` | 1 | Yes | high | 49.7 | 50.5 | 0.7 | 13.0x | 4134 |
| 3 | `izzo2015` | 1 | Yes | high | 61.2 | 62.2 | 1.0 | 10.6x | 3369 |
| 4 | `gooding1990` | 1 | Yes | high | 96.4 | 99.7 | 2.6 | 6.7x | 2120 |
| 5 | `arora2013` | 1 | Yes | high | 157.2 | 160.0 | 4.2 | 4.1x | 1303 |
| 6 | `negrete2024` | 1 | Yes | high | 196.8 | 201.5 | 8.9 | 3.3x | 1026 |
| 7 | `delatorre2018` | 1 | Yes | high | 646.4 | 647.2 | 10.2 | 1.0x | 317 |
| 1 | `mcelreath2025` | 1 | Yes | low | 47.8 | 48.8 | 0.7 | 13.0x | 4285 |
| 2 | `der2011` | 1 | Yes | low | 49.9 | 50.8 | 0.7 | 12.5x | 4087 |
| 3 | `izzo2015` | 1 | Yes | low | 61.6 | 62.7 | 1.0 | 10.1x | 3310 |
| 4 | `gooding1990` | 1 | Yes | low | 97.7 | 100.1 | 2.2 | 6.4x | 2067 |
| 5 | `arora2013` | 1 | Yes | low | 162.2 | 165.5 | 4.7 | 3.8x | 1263 |
| 6 | `negrete2024` | 1 | Yes | low | 196.5 | 199.4 | 3.8 | 3.2x | 1025 |
| 7 | `delatorre2018` | 1 | Yes | low | 623.1 | 625.2 | 9.3 | 1.0x | 329 |

### Near-tangent robustness cases

#### Near tangent

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 0 | Yes | high | 47.4 | 48.5 | 0.7 | 6.7x | 4332 |
| 2 | `arora2013` | 0 | Yes | high | 54.6 | 55.4 | 0.7 | 5.8x | 3753 |
| 3 | `izzo2015` | 0 | Yes | high | 60.0 | 61.0 | 1.0 | 5.3x | 3438 |
| 4 | `der2011` | 0 | Yes | high | 65.7 | 66.9 | 1.0 | 4.8x | 3176 |
| 5 | `gooding1990` | 0 | Yes | high | 79.5 | 81.0 | 1.2 | 4.0x | 2591 |
| 6 | `negrete2024` | 0 | Yes | high | 128.4 | 130.2 | 1.2 | 2.5x | 1572 |
| 7 | `delatorre2018` | 0 | Yes | high | 316.1 | 318.9 | 8.8 | 1.0x | 644 |

### Zero-revolution general cases

#### Curtiss book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 33.1 | 33.6 | 0.3 | 990.6x | 6196 |
| 2 | `der2011` | 0 | Yes | low | 39.3 | 40.1 | 0.5 | 834.2x | 5166 |
| 3 | `jiang2016` | 0 | Yes | low | 40.9 | 41.7 | 0.5 | 802.9x | 5028 |
| 4 | `mcelreath2025` | 0 | Yes | low | 47.1 | 51.8 | 1.7 | 696.7x | 4342 |
| 5 | `arora2013` | 0 | Yes | low | 55.7 | 56.8 | 0.7 | 588.8x | 3672 |
| 6 | `izzo2015` | 0 | Yes | low | 60.1 | 61.0 | 1.0 | 546.3x | 3411 |
| 7 | `gooding1990` | 0 | Yes | low | 79.1 | 80.5 | 1.2 | 415.0x | 2581 |
| 8 | `thorne2004` | 0 | Yes | low | 97.5 | 99.2 | 1.9 | 336.5x | 2099 |
| 9 | `avanzini2008` | 0 | Yes | low | 118.5 | 121.5 | 2.7 | 276.9x | 1719 |
| 10 | `negrete2024` | 0 | Yes | low | 125.8 | 128.3 | 1.4 | 260.7x | 1602 |
| 11 | `delatorre2018` | 0 | Yes | low | 315.2 | 318.3 | 9.2 | 104.1x | 646 |
| 12 | `pan2016` | 0 | Yes | low | 32808.5 | 32679.4 | 483.2 | 1.0x | 7 |

#### Der article I

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | high | 35.3 | 35.9 | 0.4 | 885.7x | 5749 |
| 2 | `jiang2016` | 0 | No | high | 42.5 | 43.2 | 0.5 | 736.3x | 4803 |
| 3 | `mcelreath2025` | 0 | No | high | 47.2 | 48.0 | 0.6 | 662.1x | 4361 |
| 4 | `arora2013` | 0 | No | high | 56.6 | 57.6 | 0.7 | 552.5x | 3599 |
| 5 | `der2011` | 0 | No | high | 57.2 | 58.3 | 1.0 | 546.4x | 3614 |
| 6 | `izzo2015` | 0 | No | high | 60.8 | 61.7 | 0.9 | 514.3x | 3375 |
| 7 | `gooding1990` | 0 | No | high | 80.8 | 82.5 | 1.4 | 387.1x | 2500 |
| 8 | `thorne2004` | 0 | No | high | 101.2 | 103.1 | 2.1 | 309.0x | 2020 |
| 9 | `avanzini2008` | 0 | No | high | 124.8 | 127.5 | 2.8 | 250.5x | 1631 |
| 10 | `negrete2024` | 0 | No | high | 126.8 | 128.4 | 1.2 | 246.4x | 1588 |
| 11 | `delatorre2018` | 0 | No | high | 314.5 | 318.4 | 9.5 | 99.4x | 648 |
| 12 | `pan2016` | 0 | No | high | 31257.4 | 31297.2 | 134.9 | 1.0x | 7 |
| 1 | `battin1984` | 0 | Yes | low | 36.6 | 37.1 | 0.3 | 857.8x | 5568 |
| 2 | `der2011` | 0 | Yes | low | 39.9 | 40.6 | 0.5 | 785.5x | 5099 |
| 3 | `jiang2016` | 0 | Yes | low | 42.0 | 42.9 | 0.5 | 747.1x | 4864 |
| 4 | `mcelreath2025` | 0 | Yes | low | 45.5 | 46.2 | 0.6 | 689.1x | 4319 |
| 5 | `arora2013` | 0 | Yes | low | 53.6 | 54.4 | 0.6 | 585.6x | 3801 |
| 6 | `izzo2015` | 0 | Yes | low | 59.8 | 60.7 | 0.9 | 524.5x | 3428 |
| 7 | `gooding1990` | 0 | Yes | low | 79.5 | 82.4 | 1.3 | 394.5x | 2548 |
| 8 | `thorne2004` | 0 | Yes | low | 105.6 | 107.5 | 2.2 | 296.9x | 1952 |
| 9 | `avanzini2008` | 0 | Yes | low | 125.4 | 127.8 | 2.7 | 250.2x | 1635 |
| 10 | `negrete2024` | 0 | Yes | low | 127.4 | 128.8 | 1.1 | 246.2x | 1581 |
| 11 | `delatorre2018` | 0 | Yes | low | 288.7 | 289.2 | 9.3 | 108.6x | 721 |
| 12 | `pan2016` | 0 | Yes | low | 31364.5 | 31397.0 | 69.0 | 1.0x | 7 |

#### Der article II

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | high | 35.5 | 36.0 | 0.3 | 1016.3x | 5705 |
| 2 | `der2011` | 0 | No | high | 40.3 | 41.0 | 0.5 | 894.3x | 5042 |
| 3 | `jiang2016` | 0 | No | high | 40.7 | 41.4 | 0.5 | 885.5x | 5017 |
| 4 | `mcelreath2025` | 0 | No | high | 47.2 | 48.0 | 0.6 | 763.6x | 4324 |
| 5 | `izzo2015` | 0 | No | high | 60.2 | 61.2 | 1.0 | 599.0x | 3415 |
| 6 | `arora2013` | 0 | No | high | 63.3 | 64.5 | 0.8 | 569.4x | 3219 |
| 7 | `gooding1990` | 0 | No | high | 81.8 | 83.3 | 1.2 | 440.4x | 2486 |
| 8 | `thorne2004` | 0 | No | high | 90.2 | 92.2 | 2.0 | 399.5x | 2274 |
| 9 | `avanzini2008` | 0 | No | high | 121.3 | 125.0 | 3.1 | 297.0x | 1710 |
| 10 | `negrete2024` | 0 | No | high | 128.2 | 130.3 | 1.4 | 281.1x | 1572 |
| 11 | `delatorre2018` | 0 | No | high | 296.7 | 299.8 | 8.9 | 121.5x | 685 |
| 12 | `pan2016` | 0 | No | high | 36042.0 | 36043.5 | 146.6 | 1.0x | 6 |
| 1 | `battin1984` | 0 | Yes | high | 37.2 | 37.8 | 0.4 | 976.3x | 5472 |
| 2 | `der2011` | 0 | Yes | high | 39.7 | 40.4 | 0.5 | 912.6x | 5130 |
| 3 | `jiang2016` | 0 | Yes | high | 40.5 | 41.2 | 0.5 | 896.1x | 5024 |
| 4 | `mcelreath2025` | 0 | Yes | high | 47.2 | 48.2 | 0.6 | 768.2x | 4392 |
| 5 | `arora2013` | 0 | Yes | high | 53.9 | 54.9 | 0.6 | 672.5x | 3771 |
| 6 | `izzo2015` | 0 | Yes | high | 60.3 | 61.2 | 0.9 | 602.0x | 3415 |
| 7 | `gooding1990` | 0 | Yes | high | 81.1 | 82.7 | 1.3 | 447.6x | 2498 |
| 8 | `thorne2004` | 0 | Yes | high | 89.0 | 90.7 | 1.9 | 407.5x | 2291 |
| 9 | `avanzini2008` | 0 | Yes | high | 128.3 | 131.7 | 3.6 | 282.8x | 1588 |
| 10 | `negrete2024` | 0 | Yes | high | 128.3 | 129.8 | 1.2 | 282.8x | 1570 |
| 11 | `delatorre2018` | 0 | Yes | high | 306.9 | 321.9 | 9.7 | 118.2x | 659 |
| 12 | `pan2016` | 0 | Yes | high | 36276.9 | 36228.8 | 305.5 | 1.0x | 6 |

### Zero-revolution hyperbolic cases

#### GMAT hyperbolic

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | low | 35.4 | 35.9 | 0.3 | 1027.1x | 5726 |
| 2 | `der2011` | 0 | No | low | 39.7 | 40.4 | 0.4 | 917.5x | 5084 |
| 3 | `jiang2016` | 0 | No | low | 40.8 | 41.5 | 0.5 | 893.1x | 4998 |
| 4 | `mcelreath2025` | 0 | No | low | 48.5 | 49.3 | 0.6 | 751.1x | 4195 |
| 5 | `izzo2015` | 0 | No | low | 60.3 | 61.1 | 0.9 | 604.1x | 3392 |
| 6 | `arora2013` | 0 | No | low | 62.6 | 63.5 | 1.0 | 581.2x | 3284 |
| 7 | `gooding1990` | 0 | No | low | 80.8 | 82.2 | 1.0 | 450.3x | 2552 |
| 8 | `thorne2004` | 0 | No | low | 101.2 | 103.2 | 2.1 | 359.6x | 2030 |
| 9 | `avanzini2008` | 0 | No | low | 126.8 | 130.3 | 3.1 | 287.0x | 1613 |
| 10 | `negrete2024` | 0 | No | low | 132.5 | 134.1 | 1.2 | 274.6x | 1521 |
| 11 | `delatorre2018` | 0 | No | low | 329.7 | 331.9 | 8.9 | 110.4x | 617 |
| 12 | `pan2016` | 0 | No | low | 36395.9 | 36396.0 | 49.4 | 1.0x | 6 |
| 1 | `battin1984` | 0 | Yes | low | 34.5 | 35.0 | 0.3 | 1059.8x | 5902 |
| 2 | `der2011` | 0 | Yes | low | 39.7 | 40.5 | 0.5 | 920.5x | 5128 |
| 3 | `jiang2016` | 0 | Yes | low | 40.1 | 40.8 | 0.4 | 912.5x | 5079 |
| 4 | `mcelreath2025` | 0 | Yes | low | 48.0 | 48.8 | 0.6 | 761.9x | 4242 |
| 5 | `arora2013` | 0 | Yes | low | 53.4 | 54.2 | 0.9 | 685.1x | 3912 |
| 6 | `izzo2015` | 0 | Yes | low | 59.8 | 60.7 | 0.9 | 611.6x | 3426 |
| 7 | `gooding1990` | 0 | Yes | low | 80.7 | 82.0 | 1.1 | 453.2x | 2544 |
| 8 | `thorne2004` | 0 | Yes | low | 99.6 | 101.7 | 2.0 | 367.1x | 2065 |
| 9 | `avanzini2008` | 0 | Yes | low | 125.0 | 127.8 | 2.8 | 292.5x | 1636 |
| 10 | `negrete2024` | 0 | Yes | low | 132.7 | 135.4 | 1.5 | 275.5x | 1524 |
| 11 | `delatorre2018` | 0 | Yes | low | 335.0 | 334.9 | 13.2 | 109.2x | 619 |
| 12 | `pan2016` | 0 | Yes | low | 36563.6 | 36566.5 | 135.8 | 1.0x | 6 |

### Zero-revolution nominal cases

#### Battin book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `gauss1809` | 0 | Yes | low | 33.7 | 34.7 | 0.5 | 1059.5x | 6087 |
| 2 | `battin1984` | 0 | Yes | low | 37.0 | 38.0 | 0.4 | 965.9x | 5494 |
| 3 | `der2011` | 0 | Yes | low | 50.5 | 52.2 | 1.2 | 707.6x | 4030 |
| 4 | `jiang2016` | 0 | Yes | low | 51.0 | 52.7 | 1.1 | 700.2x | 4028 |
| 5 | `mcelreath2025` | 0 | Yes | low | 60.4 | 62.5 | 1.3 | 591.2x | 3436 |
| 6 | `arora2013` | 0 | Yes | low | 78.9 | 81.2 | 1.9 | 452.7x | 2634 |
| 7 | `izzo2015` | 0 | Yes | low | 83.0 | 85.5 | 1.7 | 430.4x | 2439 |
| 8 | `gooding1990` | 0 | Yes | low | 104.1 | 114.2 | 2.9 | 342.9x | 1980 |
| 9 | `thorne2004` | 0 | Yes | low | 126.8 | 130.7 | 5.1 | 281.7x | 1628 |
| 10 | `negrete2024` | 0 | Yes | low | 127.4 | 130.8 | 4.2 | 280.2x | 1583 |
| 11 | `vallado2013` | 0 | Yes | low | 131.6 | 134.5 | 4.4 | 271.4x | 1530 |
| 12 | `avanzini2008` | 0 | Yes | low | 156.1 | 161.6 | 5.1 | 228.8x | 1315 |
| 13 | `delatorre2018` | 0 | Yes | low | 309.0 | 312.3 | 14.9 | 115.6x | 666 |
| 14 | `pan2016` | 0 | Yes | low | 35707.4 | 35730.6 | 209.5 | 1.0x | 6 |

#### Vallado book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 39.1 | 40.6 | 0.6 | 8.8x | 5273 |
| 2 | `gauss1809` | 0 | Yes | low | 52.6 | 54.2 | 1.0 | 6.5x | 3860 |
| 3 | `jiang2016` | 0 | Yes | low | 56.1 | 58.2 | 1.4 | 6.1x | 3686 |
| 4 | `mcelreath2025` | 0 | Yes | low | 61.5 | 63.9 | 1.4 | 5.6x | 3402 |
| 5 | `vallado2013` | 0 | Yes | low | 67.6 | 70.3 | 1.0 | 5.1x | 2990 |
| 6 | `der2011` | 0 | Yes | low | 69.7 | 72.2 | 1.5 | 4.9x | 2919 |
| 7 | `arora2013` | 0 | Yes | low | 78.6 | 80.9 | 2.0 | 4.4x | 2617 |
| 8 | `izzo2015` | 0 | Yes | low | 82.3 | 84.8 | 1.7 | 4.2x | 2508 |
| 9 | `gooding1990` | 0 | Yes | low | 106.7 | 110.3 | 3.1 | 3.2x | 1894 |
| 10 | `pan2016` | 0 | Yes | low | 115.9 | 119.4 | 3.7 | 3.0x | 1775 |
| 11 | `thorne2004` | 0 | Yes | low | 125.1 | 129.0 | 5.0 | 2.7x | 1641 |
| 12 | `negrete2024` | 0 | Yes | low | 128.4 | 131.9 | 4.4 | 2.7x | 1569 |
| 13 | `avanzini2008` | 0 | Yes | low | 164.9 | 172.2 | 7.4 | 2.1x | 1246 |
| 14 | `delatorre2018` | 0 | Yes | low | 343.7 | 347.3 | 14.1 | 1.0x | 595 |

<!-- performance-comparison:end -->
