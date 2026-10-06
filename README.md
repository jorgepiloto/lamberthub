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

- Generated: 2026-10-06 14:08 UTC
- Commit: `e368bbb2f9bd`
- Environment: Linux, Python 3.12.14, AMD EPYC 7763 64-Core Processor

Times are in microseconds (lower is better).
**Speedup** is relative to the slowest solver for each case.
All solvers are JIT-warmed before timing begins.

### Near-minimum-energy multi-revolution cases

#### Near-minimum-energy

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 1 | Yes | high | 52.7 | 54.9 | 1.4 | 13.3x | 3974 |
| 2 | `der2011` | 1 | Yes | high | 55.5 | 57.2 | 1.0 | 12.7x | 3706 |
| 3 | `izzo2015` | 1 | Yes | high | 83.9 | 86.6 | 2.3 | 8.4x | 2447 |
| 4 | `gooding1990` | 1 | Yes | high | 117.7 | 122.6 | 4.4 | 6.0x | 1757 |
| 5 | `negrete2024` | 1 | Yes | high | 190.4 | 195.7 | 7.2 | 3.7x | 1059 |
| 6 | `arora2013` | 1 | Yes | high | 240.5 | 245.8 | 17.2 | 2.9x | 852 |
| 7 | `delatorre2018` | 1 | Yes | high | 702.6 | 704.9 | 12.4 | 1.0x | 296 |
| 1 | `mcelreath2025` | 1 | Yes | low | 53.1 | 55.1 | 1.2 | 14.0x | 3893 |
| 2 | `der2011` | 1 | Yes | low | 82.3 | 85.0 | 1.8 | 9.0x | 2512 |
| 3 | `izzo2015` | 1 | Yes | low | 83.2 | 86.0 | 2.0 | 8.9x | 2456 |
| 4 | `gooding1990` | 1 | Yes | low | 116.0 | 119.9 | 3.2 | 6.4x | 1781 |
| 5 | `negrete2024` | 1 | Yes | low | 190.2 | 194.5 | 5.1 | 3.9x | 1062 |
| 6 | `arora2013` | 1 | Yes | low | 262.4 | 268.3 | 18.0 | 2.8x | 784 |
| 7 | `delatorre2018` | 1 | Yes | low | 744.3 | 746.5 | 5.6 | 1.0x | 277 |
| 1 | `der2011` | 2 | Yes | high | 51.2 | 52.9 | 0.8 | 13.4x | 3994 |
| 2 | `mcelreath2025` | 2 | Yes | high | 53.3 | 55.3 | 1.1 | 12.8x | 3909 |
| 3 | `izzo2015` | 2 | Yes | high | 83.7 | 86.9 | 2.0 | 8.2x | 2450 |
| 4 | `gooding1990` | 2 | Yes | high | 118.2 | 122.4 | 3.4 | 5.8x | 1732 |
| 5 | `negrete2024` | 2 | Yes | high | 190.5 | 194.7 | 5.3 | 3.6x | 1059 |
| 6 | `arora2013` | 2 | Yes | high | 241.9 | 247.5 | 17.6 | 2.8x | 851 |
| 7 | `delatorre2018` | 2 | Yes | high | 685.0 | 685.7 | 11.7 | 1.0x | 299 |
| 1 | `der2011` | 2 | Yes | low | 51.0 | 52.6 | 0.8 | 14.1x | 4007 |
| 2 | `mcelreath2025` | 2 | Yes | low | 52.6 | 54.4 | 1.1 | 13.6x | 3890 |
| 3 | `izzo2015` | 2 | Yes | low | 84.7 | 88.1 | 3.1 | 8.5x | 2458 |
| 4 | `gooding1990` | 2 | Yes | low | 117.5 | 121.7 | 3.4 | 6.1x | 1731 |
| 5 | `negrete2024` | 2 | Yes | low | 190.3 | 194.5 | 5.1 | 3.8x | 1059 |
| 6 | `arora2013` | 2 | Yes | low | 257.4 | 262.9 | 18.5 | 2.8x | 799 |
| 7 | `delatorre2018` | 2 | Yes | low | 717.3 | 722.3 | 12.3 | 1.0x | 285 |

### One-revolution branch cases

#### Der article II

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 1 | No | high | 24.8 | 25.0 | 0.3 | 14.8x | 8420 |
| 2 | `der2011` | 1 | No | high | 27.5 | 27.7 | 0.6 | 13.3x | 7735 |
| 3 | `izzo2015` | 1 | No | high | 34.0 | 34.5 | 0.8 | 10.8x | 6108 |
| 4 | `gooding1990` | 1 | No | high | 49.2 | 50.0 | 1.1 | 7.4x | 4194 |
| 5 | `arora2013` | 1 | No | high | 85.4 | 86.1 | 2.0 | 4.3x | 2522 |
| 6 | `negrete2024` | 1 | No | high | 124.3 | 125.0 | 3.1 | 2.9x | 1627 |
| 7 | `delatorre2018` | 1 | No | high | 366.2 | 367.8 | 6.0 | 1.0x | 567 |
| 1 | `mcelreath2025` | 1 | No | low | 24.6 | 24.8 | 0.2 | 15.3x | 8253 |
| 2 | `der2011` | 1 | No | low | 26.9 | 27.1 | 0.5 | 14.0x | 7608 |
| 3 | `izzo2015` | 1 | No | low | 34.0 | 34.3 | 1.0 | 11.1x | 6449 |
| 4 | `gooding1990` | 1 | No | low | 48.5 | 49.0 | 1.0 | 7.8x | 4181 |
| 5 | `arora2013` | 1 | No | low | 81.9 | 82.5 | 2.2 | 4.6x | 2471 |
| 6 | `negrete2024` | 1 | No | low | 125.3 | 125.7 | 4.1 | 3.0x | 1646 |
| 7 | `delatorre2018` | 1 | No | low | 376.8 | 376.5 | 11.6 | 1.0x | 590 |
| 1 | `mcelreath2025` | 1 | Yes | high | 24.9 | 25.2 | 0.2 | 15.0x | 8175 |
| 2 | `der2011` | 1 | Yes | high | 27.0 | 27.3 | 0.8 | 13.8x | 7792 |
| 3 | `izzo2015` | 1 | Yes | high | 32.9 | 33.3 | 0.7 | 11.3x | 6409 |
| 4 | `gooding1990` | 1 | Yes | high | 50.8 | 51.8 | 1.6 | 7.3x | 4062 |
| 5 | `arora2013` | 1 | Yes | high | 91.2 | 91.8 | 3.0 | 4.1x | 2279 |
| 6 | `negrete2024` | 1 | Yes | high | 127.9 | 127.7 | 3.7 | 2.9x | 1666 |
| 7 | `delatorre2018` | 1 | Yes | high | 372.5 | 373.6 | 8.3 | 1.0x | 560 |
| 1 | `mcelreath2025` | 1 | Yes | low | 24.3 | 24.9 | 0.3 | 14.6x | 8495 |
| 2 | `der2011` | 1 | Yes | low | 27.0 | 27.2 | 0.7 | 13.2x | 7687 |
| 3 | `izzo2015` | 1 | Yes | low | 31.7 | 32.0 | 0.6 | 11.2x | 6314 |
| 4 | `gooding1990` | 1 | Yes | low | 51.2 | 51.7 | 0.9 | 6.9x | 4063 |
| 5 | `arora2013` | 1 | Yes | low | 94.3 | 95.2 | 3.2 | 3.8x | 2239 |
| 6 | `negrete2024` | 1 | Yes | low | 122.9 | 123.9 | 2.1 | 2.9x | 1709 |
| 7 | `delatorre2018` | 1 | Yes | low | 354.9 | 360.7 | 6.6 | 1.0x | 584 |

### Near-tangent robustness cases

#### Near tangent

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 0 | Yes | high | 26.1 | 26.3 | 0.4 | 7.1x | 7820 |
| 2 | `arora2013` | 0 | Yes | high | 30.5 | 30.8 | 0.2 | 6.1x | 6653 |
| 3 | `izzo2015` | 0 | Yes | high | 33.3 | 34.6 | 0.8 | 5.6x | 6141 |
| 4 | `der2011` | 0 | Yes | high | 35.2 | 35.9 | 0.5 | 5.3x | 5836 |
| 5 | `gooding1990` | 0 | Yes | high | 43.2 | 43.7 | 0.6 | 4.3x | 4739 |
| 6 | `negrete2024` | 0 | Yes | high | 81.8 | 81.8 | 3.6 | 2.3x | 2523 |
| 7 | `delatorre2018` | 0 | Yes | high | 186.1 | 187.5 | 4.5 | 1.0x | 1093 |

### Zero-revolution general cases

#### Curtiss book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 35.4 | 36.5 | 0.4 | 927.5x | 5750 |
| 2 | `der2011` | 0 | Yes | low | 42.3 | 43.7 | 0.7 | 774.8x | 4830 |
| 3 | `jiang2016` | 0 | Yes | low | 44.6 | 46.1 | 0.8 | 735.2x | 4586 |
| 4 | `mcelreath2025` | 0 | Yes | low | 54.0 | 55.9 | 1.1 | 607.5x | 3834 |
| 5 | `arora2013` | 0 | Yes | low | 63.2 | 64.8 | 1.2 | 518.5x | 3288 |
| 6 | `izzo2015` | 0 | Yes | low | 82.6 | 85.7 | 2.2 | 397.0x | 2534 |
| 7 | `gooding1990` | 0 | Yes | low | 92.6 | 95.7 | 1.7 | 354.0x | 2183 |
| 8 | `thorne2004` | 0 | Yes | low | 114.7 | 118.5 | 2.7 | 285.8x | 1782 |
| 9 | `negrete2024` | 0 | Yes | low | 123.9 | 127.6 | 3.5 | 264.7x | 1627 |
| 10 | `avanzini2008` | 0 | Yes | low | 149.2 | 154.2 | 4.6 | 219.8x | 1357 |
| 11 | `delatorre2018` | 0 | Yes | low | 336.5 | 339.7 | 13.4 | 97.5x | 607 |
| 12 | `pan2016` | 0 | Yes | low | 32794.8 | 32749.4 | 230.6 | 1.0x | 7 |

#### Der article I

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | high | 36.6 | 37.6 | 0.4 | 857.1x | 5507 |
| 2 | `jiang2016` | 0 | No | high | 46.5 | 48.0 | 0.9 | 675.1x | 4426 |
| 3 | `mcelreath2025` | 0 | No | high | 53.0 | 54.8 | 0.9 | 591.6x | 3887 |
| 4 | `der2011` | 0 | No | high | 63.2 | 65.7 | 1.7 | 496.2x | 3260 |
| 5 | `arora2013` | 0 | No | high | 63.3 | 65.3 | 1.1 | 495.6x | 3231 |
| 6 | `izzo2015` | 0 | No | high | 81.3 | 84.1 | 1.9 | 385.8x | 2531 |
| 7 | `gooding1990` | 0 | No | high | 96.3 | 99.5 | 2.0 | 325.8x | 2116 |
| 8 | `thorne2004` | 0 | No | high | 116.3 | 120.0 | 3.5 | 269.8x | 1739 |
| 9 | `negrete2024` | 0 | No | high | 125.1 | 128.6 | 3.6 | 250.8x | 1610 |
| 10 | `avanzini2008` | 0 | No | high | 152.1 | 157.5 | 5.2 | 206.3x | 1342 |
| 11 | `delatorre2018` | 0 | No | high | 335.6 | 339.0 | 13.3 | 93.5x | 610 |
| 12 | `pan2016` | 0 | No | high | 31377.0 | 31340.2 | 224.3 | 1.0x | 7 |
| 1 | `battin1984` | 0 | Yes | low | 38.1 | 39.0 | 0.4 | 830.8x | 5289 |
| 2 | `der2011` | 0 | Yes | low | 42.8 | 44.3 | 0.8 | 739.8x | 4812 |
| 3 | `jiang2016` | 0 | Yes | low | 45.9 | 47.3 | 0.8 | 689.7x | 4501 |
| 4 | `mcelreath2025` | 0 | Yes | low | 53.1 | 54.9 | 1.1 | 595.9x | 3872 |
| 5 | `arora2013` | 0 | Yes | low | 61.1 | 63.5 | 1.2 | 518.2x | 3360 |
| 6 | `izzo2015` | 0 | Yes | low | 82.8 | 85.6 | 1.8 | 382.2x | 2526 |
| 7 | `gooding1990` | 0 | Yes | low | 92.6 | 95.8 | 1.8 | 342.0x | 2185 |
| 8 | `thorne2004` | 0 | Yes | low | 121.9 | 125.7 | 3.8 | 259.8x | 1677 |
| 9 | `negrete2024` | 0 | Yes | low | 125.1 | 128.1 | 3.4 | 253.1x | 1609 |
| 10 | `avanzini2008` | 0 | Yes | low | 156.5 | 161.6 | 4.7 | 202.3x | 1311 |
| 11 | `delatorre2018` | 0 | Yes | low | 301.7 | 305.6 | 14.7 | 104.9x | 681 |
| 12 | `pan2016` | 0 | Yes | low | 31655.8 | 31629.1 | 202.1 | 1.0x | 7 |

#### Der article II

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | high | 36.9 | 37.9 | 0.4 | 982.2x | 5495 |
| 2 | `der2011` | 0 | No | high | 42.4 | 43.7 | 0.8 | 855.4x | 4824 |
| 3 | `jiang2016` | 0 | No | high | 44.0 | 45.4 | 0.7 | 824.2x | 4661 |
| 4 | `mcelreath2025` | 0 | No | high | 53.1 | 55.1 | 1.0 | 682.3x | 3878 |
| 5 | `arora2013` | 0 | No | high | 72.3 | 74.2 | 1.2 | 501.6x | 2843 |
| 6 | `izzo2015` | 0 | No | high | 82.8 | 85.5 | 1.9 | 437.5x | 2491 |
| 7 | `gooding1990` | 0 | No | high | 95.5 | 99.7 | 2.1 | 379.6x | 2139 |
| 8 | `thorne2004` | 0 | No | high | 106.2 | 109.6 | 2.8 | 341.4x | 1939 |
| 9 | `negrete2024` | 0 | No | high | 125.3 | 128.3 | 3.3 | 289.3x | 1606 |
| 10 | `avanzini2008` | 0 | No | high | 149.7 | 155.6 | 4.2 | 242.0x | 1367 |
| 11 | `delatorre2018` | 0 | No | high | 313.5 | 317.2 | 13.0 | 115.6x | 644 |
| 12 | `pan2016` | 0 | No | high | 36240.9 | 36231.3 | 450.8 | 1.0x | 6 |
| 1 | `battin1984` | 0 | Yes | high | 38.6 | 39.7 | 0.5 | 936.7x | 5254 |
| 2 | `der2011` | 0 | Yes | high | 42.5 | 44.0 | 0.8 | 851.2x | 4805 |
| 3 | `jiang2016` | 0 | Yes | high | 44.3 | 45.7 | 0.8 | 817.5x | 4613 |
| 4 | `mcelreath2025` | 0 | Yes | high | 52.9 | 54.8 | 1.0 | 684.2x | 3870 |
| 5 | `arora2013` | 0 | Yes | high | 60.3 | 62.8 | 1.1 | 599.9x | 3380 |
| 6 | `izzo2015` | 0 | Yes | high | 81.8 | 84.7 | 1.9 | 442.4x | 2491 |
| 7 | `gooding1990` | 0 | Yes | high | 95.1 | 98.3 | 1.7 | 380.3x | 2145 |
| 8 | `thorne2004` | 0 | Yes | high | 104.5 | 108.0 | 2.3 | 346.1x | 1956 |
| 9 | `negrete2024` | 0 | Yes | high | 125.1 | 127.9 | 3.3 | 289.2x | 1606 |
| 10 | `avanzini2008` | 0 | Yes | high | 155.6 | 160.9 | 5.9 | 232.5x | 1310 |
| 11 | `delatorre2018` | 0 | Yes | high | 328.5 | 333.7 | 14.1 | 110.1x | 617 |
| 12 | `pan2016` | 0 | Yes | high | 36177.6 | 36125.8 | 278.4 | 1.0x | 6 |

### Zero-revolution hyperbolic cases

#### GMAT hyperbolic

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | low | 36.9 | 38.1 | 0.5 | 1048.1x | 5503 |
| 2 | `der2011` | 0 | No | low | 42.9 | 44.5 | 1.0 | 902.5x | 4789 |
| 3 | `jiang2016` | 0 | No | low | 43.7 | 45.3 | 1.0 | 886.5x | 4737 |
| 4 | `mcelreath2025` | 0 | No | low | 54.8 | 56.6 | 1.2 | 706.7x | 3778 |
| 5 | `arora2013` | 0 | No | low | 70.0 | 71.9 | 1.5 | 552.8x | 3010 |
| 6 | `izzo2015` | 0 | No | low | 82.6 | 85.6 | 2.1 | 468.5x | 2485 |
| 7 | `gooding1990` | 0 | No | low | 96.4 | 99.6 | 2.3 | 401.6x | 2131 |
| 8 | `thorne2004` | 0 | No | low | 118.2 | 122.3 | 4.5 | 327.4x | 1691 |
| 9 | `negrete2024` | 0 | No | low | 129.9 | 133.1 | 3.8 | 297.9x | 1546 |
| 10 | `avanzini2008` | 0 | No | low | 157.6 | 162.5 | 5.6 | 245.7x | 1300 |
| 11 | `delatorre2018` | 0 | No | low | 360.4 | 363.1 | 13.6 | 107.4x | 569 |
| 12 | `pan2016` | 0 | No | low | 38716.1 | 38724.5 | 143.1 | 1.0x | 6 |
| 1 | `battin1984` | 0 | Yes | low | 37.3 | 38.5 | 0.6 | 999.5x | 5477 |
| 2 | `der2011` | 0 | Yes | low | 42.8 | 44.3 | 1.0 | 870.2x | 4843 |
| 3 | `jiang2016` | 0 | Yes | low | 43.5 | 45.4 | 1.0 | 855.8x | 4750 |
| 4 | `mcelreath2025` | 0 | Yes | low | 55.4 | 57.5 | 1.4 | 673.0x | 3752 |
| 5 | `arora2013` | 0 | Yes | low | 59.8 | 61.6 | 1.1 | 623.2x | 3487 |
| 6 | `izzo2015` | 0 | Yes | low | 81.8 | 84.7 | 1.9 | 455.3x | 2505 |
| 7 | `gooding1990` | 0 | Yes | low | 96.8 | 100.5 | 4.8 | 384.7x | 2107 |
| 8 | `thorne2004` | 0 | Yes | low | 117.4 | 121.2 | 4.3 | 317.4x | 1737 |
| 9 | `negrete2024` | 0 | Yes | low | 130.4 | 133.9 | 4.1 | 285.6x | 1545 |
| 10 | `avanzini2008` | 0 | Yes | low | 155.9 | 161.0 | 5.2 | 239.0x | 1311 |
| 11 | `delatorre2018` | 0 | Yes | low | 360.5 | 365.3 | 14.4 | 103.3x | 572 |
| 12 | `pan2016` | 0 | Yes | low | 37251.3 | 37398.5 | 258.7 | 1.0x | 6 |

### Zero-revolution nominal cases

#### Battin book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `gauss1809` | 0 | Yes | low | 15.5 | 15.6 | 0.3 | 1104.7x | 13173 |
| 2 | `battin1984` | 0 | Yes | low | 17.7 | 17.7 | 0.6 | 966.0x | 11796 |
| 3 | `der2011` | 0 | Yes | low | 19.0 | 20.5 | 0.1 | 900.3x | 10606 |
| 4 | `jiang2016` | 0 | Yes | low | 19.9 | 20.0 | 0.2 | 860.0x | 10538 |
| 5 | `mcelreath2025` | 0 | Yes | low | 23.6 | 24.1 | 0.2 | 723.7x | 8645 |
| 6 | `izzo2015` | 0 | Yes | low | 31.5 | 31.8 | 0.5 | 542.4x | 6550 |
| 7 | `arora2013` | 0 | Yes | low | 35.2 | 35.7 | 0.2 | 486.3x | 5780 |
| 8 | `gooding1990` | 0 | Yes | low | 40.2 | 40.5 | 0.3 | 425.8x | 4877 |
| 9 | `thorne2004` | 0 | Yes | low | 49.7 | 50.3 | 0.5 | 343.8x | 4088 |
| 10 | `avanzini2008` | 0 | Yes | low | 57.0 | 57.4 | 0.3 | 300.1x | 3541 |
| 11 | `vallado2013` | 0 | Yes | low | 68.3 | 68.9 | 0.5 | 250.6x | 2949 |
| 12 | `negrete2024` | 0 | Yes | low | 73.3 | 73.8 | 0.3 | 233.3x | 2739 |
| 13 | `delatorre2018` | 0 | Yes | low | 159.0 | 159.1 | 2.2 | 107.6x | 1340 |
| 14 | `pan2016` | 0 | Yes | low | 17104.8 | 17052.3 | 407.1 | 1.0x | 13 |

#### Vallado book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 18.7 | 18.8 | 0.1 | 9.4x | 10949 |
| 2 | `jiang2016` | 0 | Yes | low | 23.0 | 23.2 | 0.2 | 7.7x | 8833 |
| 3 | `mcelreath2025` | 0 | Yes | low | 23.9 | 24.3 | 0.1 | 7.4x | 8232 |
| 4 | `gauss1809` | 0 | Yes | low | 25.7 | 25.9 | 0.3 | 6.9x | 8175 |
| 5 | `der2011` | 0 | Yes | low | 25.8 | 26.1 | 0.6 | 6.8x | 8008 |
| 6 | `izzo2015` | 0 | Yes | low | 31.9 | 32.5 | 0.7 | 5.5x | 6449 |
| 7 | `vallado2013` | 0 | Yes | low | 33.2 | 33.5 | 0.4 | 5.3x | 6147 |
| 8 | `arora2013` | 0 | Yes | low | 35.3 | 35.7 | 0.6 | 5.0x | 5747 |
| 9 | `gooding1990` | 0 | Yes | low | 40.4 | 40.7 | 0.2 | 4.4x | 4997 |
| 10 | `pan2016` | 0 | Yes | low | 49.7 | 50.2 | 0.9 | 3.5x | 4111 |
| 11 | `thorne2004` | 0 | Yes | low | 52.7 | 53.5 | 1.2 | 3.3x | 4120 |
| 12 | `avanzini2008` | 0 | Yes | low | 62.0 | 62.7 | 1.0 | 2.8x | 3262 |
| 13 | `negrete2024` | 0 | Yes | low | 76.2 | 77.1 | 0.8 | 2.3x | 2719 |
| 14 | `delatorre2018` | 0 | Yes | low | 176.4 | 177.9 | 2.2 | 1.0x | 1193 |

<!-- performance-comparison:end -->
