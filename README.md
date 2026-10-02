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

- Generated: 2026-10-02 16:56 UTC
- Commit: `bd5b7d9bfce8`
- Environment: Linux, Python 3.12.14, AMD EPYC 9V74 80-Core Processor

Times are in microseconds (lower is better).
**Speedup** is relative to the slowest solver for each case.
All solvers are JIT-warmed before timing begins.

### Near-minimum-energy multi-revolution cases

#### Near-minimum-energy

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 1 | Yes | high | 43.8 | 44.5 | 0.6 | 14.7x | 4651 |
| 2 | `der2011` | 1 | Yes | high | 51.1 | 52.0 | 0.9 | 12.6x | 4050 |
| 3 | `izzo2015` | 1 | Yes | high | 59.3 | 60.3 | 0.9 | 10.8x | 3444 |
| 4 | `gooding1990` | 1 | Yes | high | 88.4 | 90.3 | 1.8 | 7.3x | 2297 |
| 5 | `negrete2024` | 1 | Yes | high | 190.0 | 192.2 | 2.2 | 3.4x | 1060 |
| 6 | `arora2013` | 1 | Yes | high | 207.0 | 210.3 | 7.4 | 3.1x | 988 |
| 7 | `delatorre2018` | 1 | Yes | high | 643.0 | 656.7 | 9.8 | 1.0x | 320 |
| 1 | `mcelreath2025` | 1 | Yes | low | 44.5 | 45.2 | 0.5 | 15.7x | 4591 |
| 2 | `izzo2015` | 1 | Yes | low | 59.6 | 60.7 | 0.9 | 11.7x | 3429 |
| 3 | `der2011` | 1 | Yes | low | 77.4 | 80.9 | 1.8 | 9.0x | 2678 |
| 4 | `gooding1990` | 1 | Yes | low | 87.8 | 89.9 | 1.5 | 7.9x | 2305 |
| 5 | `negrete2024` | 1 | Yes | low | 190.3 | 192.6 | 2.6 | 3.7x | 1059 |
| 6 | `arora2013` | 1 | Yes | low | 231.0 | 234.2 | 8.6 | 3.0x | 881 |
| 7 | `delatorre2018` | 1 | Yes | low | 697.4 | 697.3 | 9.0 | 1.0x | 292 |
| 1 | `mcelreath2025` | 2 | Yes | high | 44.3 | 45.0 | 0.5 | 14.6x | 4613 |
| 2 | `der2011` | 2 | Yes | high | 46.7 | 47.5 | 0.7 | 13.8x | 4405 |
| 3 | `izzo2015` | 2 | Yes | high | 59.7 | 60.7 | 0.9 | 10.8x | 3438 |
| 4 | `gooding1990` | 2 | Yes | high | 89.4 | 95.3 | 3.5 | 7.2x | 2275 |
| 5 | `negrete2024` | 2 | Yes | high | 191.4 | 194.1 | 2.8 | 3.4x | 1051 |
| 6 | `arora2013` | 2 | Yes | high | 211.4 | 214.1 | 5.4 | 3.1x | 964 |
| 7 | `delatorre2018` | 2 | Yes | high | 645.4 | 645.2 | 8.6 | 1.0x | 318 |
| 1 | `mcelreath2025` | 2 | Yes | low | 44.0 | 44.9 | 0.5 | 15.3x | 4609 |
| 2 | `der2011` | 2 | Yes | low | 47.5 | 48.4 | 0.8 | 14.2x | 4353 |
| 3 | `izzo2015` | 2 | Yes | low | 60.2 | 61.2 | 0.9 | 11.2x | 3392 |
| 4 | `gooding1990` | 2 | Yes | low | 90.1 | 91.9 | 1.5 | 7.5x | 2262 |
| 5 | `negrete2024` | 2 | Yes | low | 191.5 | 194.2 | 2.8 | 3.5x | 1052 |
| 6 | `arora2013` | 2 | Yes | low | 224.6 | 227.9 | 8.7 | 3.0x | 904 |
| 7 | `delatorre2018` | 2 | Yes | low | 675.1 | 678.5 | 11.2 | 1.0x | 304 |

### One-revolution branch cases

#### Der article II

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `der2011` | 1 | No | high | 50.2 | 52.4 | 1.1 | 13.3x | 4078 |
| 2 | `mcelreath2025` | 1 | No | high | 53.7 | 55.9 | 1.4 | 12.4x | 3890 |
| 3 | `izzo2015` | 1 | No | high | 81.9 | 84.9 | 2.2 | 8.1x | 2495 |
| 4 | `gooding1990` | 1 | No | high | 118.7 | 123.3 | 4.5 | 5.6x | 1751 |
| 5 | `arora2013` | 1 | No | high | 169.1 | 174.3 | 7.3 | 3.9x | 1209 |
| 6 | `negrete2024` | 1 | No | high | 189.9 | 194.4 | 7.1 | 3.5x | 1064 |
| 7 | `delatorre2018` | 1 | No | high | 666.2 | 669.6 | 14.4 | 1.0x | 310 |
| 1 | `der2011` | 1 | No | low | 50.3 | 52.0 | 1.0 | 12.8x | 4077 |
| 2 | `mcelreath2025` | 1 | No | low | 53.5 | 55.5 | 1.5 | 12.0x | 3883 |
| 3 | `izzo2015` | 1 | No | low | 81.0 | 84.1 | 2.2 | 7.9x | 2532 |
| 4 | `gooding1990` | 1 | No | low | 116.9 | 121.3 | 4.3 | 5.5x | 1742 |
| 5 | `arora2013` | 1 | No | low | 171.7 | 177.9 | 7.5 | 3.7x | 1188 |
| 6 | `negrete2024` | 1 | No | low | 189.5 | 194.0 | 7.2 | 3.4x | 1065 |
| 7 | `delatorre2018` | 1 | No | low | 643.3 | 643.5 | 14.8 | 1.0x | 319 |
| 1 | `der2011` | 1 | Yes | high | 50.6 | 53.4 | 1.3 | 13.5x | 4079 |
| 2 | `mcelreath2025` | 1 | Yes | high | 53.9 | 56.0 | 1.6 | 12.7x | 3855 |
| 3 | `izzo2015` | 1 | Yes | high | 83.0 | 86.1 | 2.7 | 8.2x | 2493 |
| 4 | `gooding1990` | 1 | Yes | high | 118.9 | 123.3 | 4.5 | 5.8x | 1713 |
| 5 | `arora2013` | 1 | Yes | high | 181.5 | 185.8 | 7.7 | 3.8x | 1135 |
| 6 | `negrete2024` | 1 | Yes | high | 191.0 | 196.4 | 9.3 | 3.6x | 1058 |
| 7 | `delatorre2018` | 1 | Yes | high | 684.9 | 686.6 | 18.1 | 1.0x | 306 |
| 1 | `der2011` | 1 | Yes | low | 50.3 | 54.5 | 1.5 | 12.9x | 4101 |
| 2 | `mcelreath2025` | 1 | Yes | low | 53.3 | 55.5 | 1.5 | 12.1x | 3908 |
| 3 | `izzo2015` | 1 | Yes | low | 82.0 | 86.6 | 2.7 | 7.9x | 2504 |
| 4 | `gooding1990` | 1 | Yes | low | 119.0 | 123.5 | 4.8 | 5.4x | 1724 |
| 5 | `arora2013` | 1 | Yes | low | 187.2 | 193.1 | 9.3 | 3.5x | 1095 |
| 6 | `negrete2024` | 1 | Yes | low | 190.5 | 194.8 | 5.5 | 3.4x | 1058 |
| 7 | `delatorre2018` | 1 | Yes | low | 648.1 | 649.0 | 14.6 | 1.0x | 318 |

### Near-tangent robustness cases

#### Near tangent

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `mcelreath2025` | 0 | Yes | high | 34.2 | 35.3 | 0.4 | 7.1x | 5949 |
| 2 | `arora2013` | 0 | Yes | high | 40.6 | 41.2 | 0.4 | 6.0x | 4992 |
| 3 | `izzo2015` | 0 | Yes | high | 46.0 | 46.9 | 0.6 | 5.3x | 4449 |
| 4 | `der2011` | 0 | Yes | high | 48.6 | 49.3 | 0.9 | 5.0x | 4263 |
| 5 | `gooding1990` | 0 | Yes | high | 59.3 | 61.4 | 1.0 | 4.1x | 3434 |
| 6 | `negrete2024` | 0 | Yes | high | 104.5 | 105.6 | 0.9 | 2.3x | 1937 |
| 7 | `delatorre2018` | 0 | Yes | high | 241.9 | 247.3 | 7.4 | 1.0x | 849 |

### Zero-revolution general cases

#### Curtiss book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 35.2 | 36.3 | 0.5 | 952.0x | 5760 |
| 2 | `der2011` | 0 | Yes | low | 42.8 | 44.5 | 1.1 | 783.2x | 4828 |
| 3 | `jiang2016` | 0 | Yes | low | 44.2 | 46.0 | 1.1 | 758.4x | 4649 |
| 4 | `mcelreath2025` | 0 | Yes | low | 52.9 | 54.8 | 1.3 | 634.7x | 3947 |
| 5 | `arora2013` | 0 | Yes | low | 63.5 | 65.2 | 1.3 | 528.5x | 3295 |
| 6 | `izzo2015` | 0 | Yes | low | 83.4 | 86.6 | 2.7 | 402.3x | 2473 |
| 7 | `gooding1990` | 0 | Yes | low | 93.2 | 96.6 | 2.1 | 359.9x | 2195 |
| 8 | `thorne2004` | 0 | Yes | low | 115.6 | 119.7 | 4.9 | 290.2x | 1771 |
| 9 | `negrete2024` | 0 | Yes | low | 123.0 | 126.3 | 3.7 | 272.8x | 1639 |
| 10 | `avanzini2008` | 0 | Yes | low | 148.3 | 153.5 | 7.1 | 226.2x | 1370 |
| 11 | `delatorre2018` | 0 | Yes | low | 334.0 | 337.4 | 13.8 | 100.4x | 618 |
| 12 | `pan2016` | 0 | Yes | low | 33547.3 | 33536.2 | 174.0 | 1.0x | 7 |

#### Der article I

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | high | 36.3 | 37.4 | 0.4 | 872.0x | 5572 |
| 2 | `jiang2016` | 0 | No | high | 45.8 | 47.6 | 1.1 | 691.4x | 4513 |
| 3 | `mcelreath2025` | 0 | No | high | 53.1 | 55.0 | 1.3 | 597.1x | 3886 |
| 4 | `arora2013` | 0 | No | high | 63.3 | 65.3 | 1.5 | 501.0x | 3280 |
| 5 | `der2011` | 0 | No | high | 64.4 | 67.1 | 1.7 | 492.0x | 3210 |
| 6 | `izzo2015` | 0 | No | high | 82.7 | 85.7 | 2.0 | 383.2x | 2487 |
| 7 | `gooding1990` | 0 | No | high | 95.5 | 98.9 | 3.4 | 331.9x | 2136 |
| 8 | `thorne2004` | 0 | No | high | 119.3 | 123.6 | 5.0 | 265.8x | 1710 |
| 9 | `negrete2024` | 0 | No | high | 124.5 | 127.6 | 3.8 | 254.5x | 1618 |
| 10 | `avanzini2008` | 0 | No | high | 154.2 | 160.0 | 7.6 | 205.5x | 1332 |
| 11 | `delatorre2018` | 0 | No | high | 333.7 | 337.1 | 13.0 | 95.0x | 615 |
| 12 | `pan2016` | 0 | No | high | 31697.2 | 31697.4 | 101.8 | 1.0x | 7 |
| 1 | `battin1984` | 0 | Yes | low | 38.1 | 39.1 | 0.4 | 834.3x | 5325 |
| 2 | `der2011` | 0 | Yes | low | 42.6 | 44.2 | 0.9 | 746.4x | 4835 |
| 3 | `jiang2016` | 0 | Yes | low | 45.2 | 46.9 | 1.0 | 703.3x | 4525 |
| 4 | `mcelreath2025` | 0 | Yes | low | 52.8 | 54.7 | 1.4 | 602.9x | 3912 |
| 5 | `arora2013` | 0 | Yes | low | 60.5 | 62.3 | 1.2 | 525.6x | 3397 |
| 6 | `izzo2015` | 0 | Yes | low | 82.6 | 85.2 | 2.1 | 385.4x | 2499 |
| 7 | `gooding1990` | 0 | Yes | low | 93.8 | 97.5 | 2.4 | 339.3x | 2167 |
| 8 | `thorne2004` | 0 | Yes | low | 123.7 | 127.5 | 4.9 | 257.2x | 1657 |
| 9 | `negrete2024` | 0 | Yes | low | 124.8 | 127.7 | 4.0 | 255.0x | 1615 |
| 10 | `avanzini2008` | 0 | Yes | low | 154.9 | 160.0 | 7.3 | 205.4x | 1324 |
| 11 | `delatorre2018` | 0 | Yes | low | 294.1 | 298.0 | 13.6 | 108.2x | 686 |
| 12 | `pan2016` | 0 | Yes | low | 31819.3 | 31858.6 | 309.7 | 1.0x | 7 |

#### Der article II

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | high | 36.2 | 37.3 | 0.5 | 1007.6x | 5533 |
| 2 | `der2011` | 0 | No | high | 42.9 | 44.4 | 1.0 | 850.4x | 4780 |
| 3 | `jiang2016` | 0 | No | high | 43.6 | 45.3 | 1.2 | 836.1x | 4701 |
| 4 | `mcelreath2025` | 0 | No | high | 52.4 | 55.2 | 1.4 | 696.8x | 3965 |
| 5 | `arora2013` | 0 | No | high | 72.9 | 76.2 | 2.1 | 500.5x | 2872 |
| 6 | `izzo2015` | 0 | No | high | 82.9 | 85.8 | 2.6 | 440.3x | 2498 |
| 7 | `gooding1990` | 0 | No | high | 95.3 | 98.7 | 2.5 | 382.7x | 2143 |
| 8 | `thorne2004` | 0 | No | high | 107.1 | 110.6 | 3.8 | 340.7x | 1905 |
| 9 | `negrete2024` | 0 | No | high | 124.6 | 127.6 | 3.8 | 292.7x | 1613 |
| 10 | `avanzini2008` | 0 | No | high | 151.8 | 157.2 | 6.6 | 240.4x | 1341 |
| 11 | `delatorre2018` | 0 | No | high | 311.9 | 315.0 | 13.7 | 117.0x | 658 |
| 12 | `pan2016` | 0 | No | high | 36481.3 | 36773.0 | 154.4 | 1.0x | 6 |
| 1 | `battin1984` | 0 | Yes | high | 38.4 | 39.4 | 0.5 | 953.6x | 5298 |
| 2 | `der2011` | 0 | Yes | high | 42.5 | 44.2 | 1.0 | 860.4x | 4799 |
| 3 | `jiang2016` | 0 | Yes | high | 44.0 | 45.8 | 1.0 | 831.0x | 4650 |
| 4 | `mcelreath2025` | 0 | Yes | high | 52.4 | 54.6 | 1.3 | 697.7x | 3950 |
| 5 | `arora2013` | 0 | Yes | high | 59.7 | 61.7 | 1.3 | 612.7x | 3432 |
| 6 | `izzo2015` | 0 | Yes | high | 82.0 | 84.8 | 2.0 | 446.2x | 2506 |
| 7 | `gooding1990` | 0 | Yes | high | 95.0 | 98.6 | 2.4 | 385.4x | 2120 |
| 8 | `thorne2004` | 0 | Yes | high | 105.9 | 109.4 | 3.6 | 345.6x | 1934 |
| 9 | `negrete2024` | 0 | Yes | high | 124.7 | 127.9 | 4.1 | 293.4x | 1616 |
| 10 | `avanzini2008` | 0 | Yes | high | 156.9 | 162.5 | 6.9 | 233.2x | 1304 |
| 11 | `delatorre2018` | 0 | Yes | high | 326.8 | 330.7 | 13.3 | 112.0x | 630 |
| 12 | `pan2016` | 0 | Yes | high | 36592.1 | 36605.5 | 157.8 | 1.0x | 6 |

### Zero-revolution hyperbolic cases

#### GMAT hyperbolic

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | No | low | 33.2 | 33.7 | 0.3 | 1036.5x | 6126 |
| 2 | `der2011` | 0 | No | low | 35.2 | 35.8 | 0.4 | 976.7x | 5804 |
| 3 | `jiang2016` | 0 | No | low | 35.8 | 36.4 | 0.4 | 962.2x | 5692 |
| 4 | `mcelreath2025` | 0 | No | low | 44.0 | 44.8 | 0.5 | 782.7x | 4640 |
| 5 | `izzo2015` | 0 | No | low | 59.4 | 60.7 | 0.9 | 579.0x | 3422 |
| 6 | `arora2013` | 0 | No | low | 60.1 | 61.0 | 0.7 | 572.1x | 3382 |
| 7 | `gooding1990` | 0 | No | low | 74.2 | 75.6 | 1.1 | 463.8x | 2734 |
| 8 | `thorne2004` | 0 | No | low | 94.5 | 96.3 | 1.7 | 363.9x | 2164 |
| 9 | `avanzini2008` | 0 | No | low | 118.7 | 121.3 | 2.9 | 289.9x | 1707 |
| 10 | `negrete2024` | 0 | No | low | 129.5 | 131.4 | 1.0 | 265.6x | 1556 |
| 11 | `delatorre2018` | 0 | No | low | 333.0 | 336.2 | 9.0 | 103.3x | 607 |
| 12 | `pan2016` | 0 | No | low | 34402.4 | 34411.5 | 89.2 | 1.0x | 6 |
| 1 | `battin1984` | 0 | Yes | low | 33.3 | 33.8 | 0.3 | 1033.2x | 6106 |
| 2 | `der2011` | 0 | Yes | low | 35.3 | 35.9 | 0.4 | 975.4x | 5742 |
| 3 | `jiang2016` | 0 | Yes | low | 35.8 | 37.1 | 0.4 | 961.2x | 5688 |
| 4 | `mcelreath2025` | 0 | Yes | low | 43.3 | 44.3 | 0.5 | 795.2x | 4735 |
| 5 | `arora2013` | 0 | Yes | low | 51.1 | 52.1 | 0.6 | 672.6x | 4004 |
| 6 | `izzo2015` | 0 | Yes | low | 58.9 | 60.1 | 0.9 | 584.0x | 3484 |
| 7 | `gooding1990` | 0 | Yes | low | 73.7 | 76.0 | 1.0 | 466.8x | 2724 |
| 8 | `thorne2004` | 0 | Yes | low | 95.0 | 97.2 | 1.8 | 362.2x | 2165 |
| 9 | `avanzini2008` | 0 | Yes | low | 118.5 | 122.0 | 3.1 | 290.3x | 1708 |
| 10 | `negrete2024` | 0 | Yes | low | 129.7 | 132.7 | 1.2 | 265.1x | 1556 |
| 11 | `delatorre2018` | 0 | Yes | low | 332.9 | 336.6 | 9.2 | 103.3x | 609 |
| 12 | `pan2016` | 0 | Yes | low | 34395.8 | 34414.8 | 49.0 | 1.0x | 6 |

### Zero-revolution nominal cases

#### Battin book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `gauss1809` | 0 | Yes | low | 29.5 | 30.4 | 0.4 | 1095.8x | 6840 |
| 2 | `battin1984` | 0 | Yes | low | 33.4 | 34.5 | 0.4 | 967.9x | 6048 |
| 3 | `der2011` | 0 | Yes | low | 42.6 | 44.3 | 1.0 | 758.1x | 4826 |
| 4 | `jiang2016` | 0 | Yes | low | 43.2 | 44.8 | 1.0 | 747.5x | 4751 |
| 5 | `mcelreath2025` | 0 | Yes | low | 53.5 | 55.6 | 1.4 | 604.5x | 3909 |
| 6 | `arora2013` | 0 | Yes | low | 73.4 | 75.8 | 1.9 | 440.2x | 2854 |
| 7 | `izzo2015` | 0 | Yes | low | 81.2 | 84.2 | 2.1 | 398.0x | 2532 |
| 8 | `gooding1990` | 0 | Yes | low | 93.9 | 97.5 | 3.0 | 344.4x | 2170 |
| 9 | `thorne2004` | 0 | Yes | low | 114.8 | 118.7 | 3.5 | 281.5x | 1776 |
| 10 | `negrete2024` | 0 | Yes | low | 123.2 | 126.4 | 3.7 | 262.4x | 1634 |
| 11 | `vallado2013` | 0 | Yes | low | 125.9 | 129.2 | 3.9 | 256.8x | 1591 |
| 12 | `avanzini2008` | 0 | Yes | low | 144.3 | 149.8 | 5.4 | 224.0x | 1420 |
| 13 | `delatorre2018` | 0 | Yes | low | 299.2 | 302.6 | 13.5 | 108.0x | 687 |
| 14 | `pan2016` | 0 | Yes | low | 32320.0 | 32287.5 | 150.4 | 1.0x | 7 |

#### Vallado book

| Rank | Solver | Revolutions | prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |
|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|
| 1 | `battin1984` | 0 | Yes | low | 35.1 | 36.3 | 0.5 | 9.5x | 5819 |
| 2 | `jiang2016` | 0 | Yes | low | 47.7 | 50.0 | 1.0 | 7.0x | 4324 |
| 3 | `gauss1809` | 0 | Yes | low | 47.9 | 49.5 | 1.1 | 7.0x | 4237 |
| 4 | `mcelreath2025` | 0 | Yes | low | 53.0 | 54.8 | 1.2 | 6.3x | 3881 |
| 5 | `der2011` | 0 | Yes | low | 60.7 | 63.2 | 1.9 | 5.5x | 3433 |
| 6 | `vallado2013` | 0 | Yes | low | 62.0 | 63.8 | 0.8 | 5.4x | 3242 |
| 7 | `arora2013` | 0 | Yes | low | 73.1 | 75.6 | 1.8 | 4.6x | 2852 |
| 8 | `izzo2015` | 0 | Yes | low | 81.2 | 84.2 | 1.9 | 4.1x | 2561 |
| 9 | `gooding1990` | 0 | Yes | low | 95.0 | 99.1 | 3.4 | 3.5x | 2139 |
| 10 | `pan2016` | 0 | Yes | low | 102.6 | 106.7 | 3.1 | 3.2x | 1991 |
| 11 | `thorne2004` | 0 | Yes | low | 114.0 | 118.3 | 4.2 | 2.9x | 1804 |
| 12 | `negrete2024` | 0 | Yes | low | 124.5 | 127.7 | 4.0 | 2.7x | 1617 |
| 13 | `avanzini2008` | 0 | Yes | low | 151.7 | 157.6 | 7.3 | 2.2x | 1366 |
| 14 | `delatorre2018` | 0 | Yes | low | 333.2 | 336.9 | 13.6 | 1.0x | 615 |

<!-- performance-comparison:end -->
