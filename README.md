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

| Group | What it tests |
|-------|---------------|
| `zero-rev-nominal` | Classic textbook cases, all solvers |
| `zero-rev-general` | General elliptic, zero-revolution solvers |
| `zero-rev-hyperbolic` | Hyperbolic transfers |
| `multi-rev-m1` | First multi-revolution branch (prograde and retrograde) |
| `multi-rev-hard` | Near-minimum-energy multi-revolution cases |
| `robust-zero-rev` | Near-tangent transfers that stress robustness |

<!-- performance-comparison:start -->

_This section is auto-generated from CI benchmark artifacts._

- Generated: 2026-05-25 14:38 UTC
- Commit: `ae2d794eb839`
- Environment: Linux, Python 3.12.13, AMD EPYC 9V74 80-Core Processor

Times are reported in microseconds. Lower values are better.

### multi-rev-hard

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `m1-prograde-high` | `mcelreath2025` | 48.2 | 49.1 | 0.6 | 4220 |
| `m1-prograde-high` | `izzo2015` | 60.6 | 61.7 | 1.0 | 3374 |
| `m1-prograde-high` | `gooding1990` | 103.7 | 105.8 | 2.2 | 1965 |
| `m1-prograde-high` | `negrete2024` | 195.0 | 197.4 | 3.2 | 1034 |
| `m1-prograde-high` | `arora2013` | 213.5 | 216.3 | 6.6 | 964 |
| `m1-prograde-high` | `der2011` | 321.4 | 324.5 | 11.0 | 640 |
| `m1-prograde-high` | `delatorre2018` | 659.3 | 662.7 | 10.2 | 315 |
| `m1-prograde-low` | `mcelreath2025` | 47.9 | 48.7 | 0.6 | 4240 |
| `m1-prograde-low` | `izzo2015` | 61.1 | 62.1 | 0.9 | 3343 |
| `m1-prograde-low` | `gooding1990` | 104.4 | 106.6 | 2.3 | 1974 |
| `m1-prograde-low` | `negrete2024` | 195.0 | 197.6 | 3.0 | 1034 |
| `m1-prograde-low` | `arora2013` | 234.9 | 238.2 | 9.9 | 874 |
| `m1-prograde-low` | `der2011` | 429.0 | 430.2 | 10.9 | 485 |
| `m1-prograde-low` | `delatorre2018` | 710.2 | 710.0 | 9.6 | 288 |
| `m2-prograde-high` | `mcelreath2025` | 47.9 | 48.9 | 0.6 | 4227 |
| `m2-prograde-high` | `izzo2015` | 61.6 | 63.8 | 1.1 | 3298 |
| `m2-prograde-high` | `gooding1990` | 106.2 | 108.4 | 2.3 | 1937 |
| `m2-prograde-high` | `negrete2024` | 196.1 | 198.4 | 2.8 | 1026 |
| `m2-prograde-high` | `arora2013` | 216.9 | 220.4 | 8.3 | 950 |
| `m2-prograde-high` | `der2011` | 268.1 | 271.3 | 9.5 | 754 |
| `m2-prograde-high` | `delatorre2018` | 655.7 | 656.2 | 10.3 | 313 |
| `m2-prograde-low` | `mcelreath2025` | 47.9 | 48.8 | 0.6 | 4242 |
| `m2-prograde-low` | `izzo2015` | 61.6 | 62.8 | 1.0 | 3288 |
| `m2-prograde-low` | `gooding1990` | 105.6 | 107.5 | 2.9 | 1923 |
| `m2-prograde-low` | `negrete2024` | 196.0 | 198.4 | 2.8 | 1027 |
| `m2-prograde-low` | `arora2013` | 229.7 | 237.6 | 10.6 | 895 |
| `m2-prograde-low` | `der2011` | 268.6 | 271.4 | 9.7 | 761 |
| `m2-prograde-low` | `delatorre2018` | 684.7 | 685.8 | 8.6 | 299 |

### multi-rev-m1

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `prograde-high` | `mcelreath2025` | 47.8 | 48.8 | 0.7 | 4299 |
| `prograde-high` | `izzo2015` | 61.7 | 62.8 | 0.9 | 3311 |
| `prograde-high` | `gooding1990` | 103.9 | 105.9 | 2.2 | 1952 |
| `prograde-high` | `arora2013` | 159.2 | 161.9 | 4.5 | 1294 |
| `prograde-high` | `negrete2024` | 197.6 | 200.3 | 2.9 | 1019 |
| `prograde-high` | `der2011` | 274.2 | 277.0 | 9.7 | 761 |
| `prograde-high` | `delatorre2018` | 653.2 | 654.6 | 9.2 | 313 |
| `prograde-low` | `mcelreath2025` | 47.5 | 48.5 | 0.6 | 4316 |
| `prograde-low` | `izzo2015` | 62.0 | 63.1 | 0.9 | 3274 |
| `prograde-low` | `gooding1990` | 104.2 | 106.6 | 2.2 | 1949 |
| `prograde-low` | `arora2013` | 163.5 | 166.3 | 4.3 | 1253 |
| `prograde-low` | `negrete2024` | 197.7 | 200.7 | 3.9 | 1018 |
| `prograde-low` | `der2011` | 253.1 | 257.1 | 10.0 | 806 |
| `prograde-low` | `delatorre2018` | 628.4 | 630.0 | 9.8 | 324 |
| `retrograde-high` | `mcelreath2025` | 47.5 | 48.3 | 0.5 | 4277 |
| `retrograde-high` | `izzo2015` | 62.1 | 63.1 | 0.9 | 3292 |
| `retrograde-high` | `gooding1990` | 103.0 | 105.2 | 2.2 | 1979 |
| `retrograde-high` | `arora2013` | 148.0 | 150.7 | 3.9 | 1391 |
| `retrograde-high` | `negrete2024` | 195.3 | 197.7 | 2.6 | 1031 |
| `retrograde-high` | `der2011` | 309.7 | 313.0 | 10.4 | 662 |
| `retrograde-high` | `delatorre2018` | 631.1 | 632.6 | 10.2 | 319 |
| `retrograde-low` | `mcelreath2025` | 47.7 | 48.6 | 0.6 | 4272 |
| `retrograde-low` | `izzo2015` | 62.7 | 64.7 | 1.0 | 3283 |
| `retrograde-low` | `gooding1990` | 102.4 | 104.5 | 2.2 | 1983 |
| `retrograde-low` | `arora2013` | 149.5 | 151.9 | 3.7 | 1368 |
| `retrograde-low` | `negrete2024` | 195.0 | 197.4 | 2.1 | 1033 |
| `retrograde-low` | `der2011` | 288.8 | 291.8 | 9.6 | 709 |
| `retrograde-low` | `delatorre2018` | 624.8 | 626.9 | 9.4 | 326 |

### robust-zero-rev

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `near-tangent-prograde-high` | `mcelreath2025` | 62.3 | 64.7 | 1.4 | 3250 |
| `near-tangent-prograde-high` | `arora2013` | 67.5 | 69.6 | 1.6 | 3097 |
| `near-tangent-prograde-high` | `izzo2015` | 87.1 | 90.0 | 2.1 | 2355 |
| `near-tangent-prograde-high` | `gooding1990` | 112.1 | 116.3 | 3.8 | 1811 |
| `near-tangent-prograde-high` | `negrete2024` | 131.4 | 134.6 | 4.9 | 1536 |
| `near-tangent-prograde-high` | `der2011` | 300.4 | 306.1 | 22.3 | 685 |
| `near-tangent-prograde-high` | `delatorre2018` | 354.8 | 359.0 | 15.0 | 579 |

### zero-rev-general

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `curtiss-book` | `battin1984` | 33.5 | 34.0 | 0.3 | 6120 |
| `curtiss-book` | `mcelreath2025` | 47.2 | 48.0 | 0.6 | 4330 |
| `curtiss-book` | `arora2013` | 56.6 | 57.4 | 0.6 | 3614 |
| `curtiss-book` | `izzo2015` | 61.2 | 62.5 | 1.1 | 3344 |
| `curtiss-book` | `gooding1990` | 84.0 | 85.6 | 1.3 | 2398 |
| `curtiss-book` | `thorne2004` | 99.3 | 101.3 | 2.3 | 2075 |
| `curtiss-book` | `der2011` | 103.3 | 105.0 | 2.2 | 1970 |
| `curtiss-book` | `avanzini2008` | 119.9 | 122.6 | 2.9 | 1699 |
| `curtiss-book` | `jiang2016` | 120.0 | 122.0 | 2.6 | 1683 |
| `curtiss-book` | `negrete2024` | 127.6 | 129.3 | 1.1 | 1583 |
| `curtiss-book` | `delatorre2018` | 321.1 | 324.2 | 8.7 | 630 |
| `curtiss-book` | `pan2016` | 32902.9 | 32937.5 | 244.7 | 7 |
| `der-article-i-prograde-low` | `battin1984` | 36.6 | 37.2 | 0.4 | 5563 |
| `der-article-i-prograde-low` | `mcelreath2025` | 47.4 | 48.3 | 0.6 | 4300 |
| `der-article-i-prograde-low` | `arora2013` | 53.7 | 54.6 | 0.6 | 3806 |
| `der-article-i-prograde-low` | `izzo2015` | 60.6 | 61.7 | 1.0 | 3393 |
| `der-article-i-prograde-low` | `gooding1990` | 81.7 | 83.2 | 1.2 | 2489 |
| `der-article-i-prograde-low` | `thorne2004` | 106.7 | 108.7 | 2.6 | 1910 |
| `der-article-i-prograde-low` | `avanzini2008` | 125.8 | 128.9 | 3.5 | 1610 |
| `der-article-i-prograde-low` | `negrete2024` | 129.6 | 131.4 | 1.2 | 1556 |
| `der-article-i-prograde-low` | `der2011` | 139.1 | 141.4 | 2.9 | 1468 |
| `der-article-i-prograde-low` | `jiang2016` | 199.1 | 215.8 | 36.1 | 1038 |
| `der-article-i-prograde-low` | `delatorre2018` | 285.5 | 289.4 | 9.1 | 712 |
| `der-article-i-prograde-low` | `pan2016` | 31046.7 | 31061.6 | 129.6 | 7 |
| `der-article-i-retrograde-high` | `battin1984` | 35.2 | 35.9 | 0.3 | 5754 |
| `der-article-i-retrograde-high` | `mcelreath2025` | 47.4 | 48.4 | 0.6 | 4323 |
| `der-article-i-retrograde-high` | `arora2013` | 56.7 | 57.6 | 0.6 | 3588 |
| `der-article-i-retrograde-high` | `izzo2015` | 61.6 | 62.7 | 0.9 | 3310 |
| `der-article-i-retrograde-high` | `gooding1990` | 83.0 | 84.8 | 1.2 | 2432 |
| `der-article-i-retrograde-high` | `thorne2004` | 101.9 | 104.0 | 2.1 | 2002 |
| `der-article-i-retrograde-high` | `avanzini2008` | 124.1 | 127.0 | 3.3 | 1629 |
| `der-article-i-retrograde-high` | `negrete2024` | 129.0 | 130.9 | 1.2 | 1565 |
| `der-article-i-retrograde-high` | `der2011` | 180.9 | 184.2 | 5.9 | 1130 |
| `der-article-i-retrograde-high` | `jiang2016` | 255.4 | 256.4 | 8.5 | 819 |
| `der-article-i-retrograde-high` | `delatorre2018` | 318.8 | 322.2 | 8.9 | 634 |
| `der-article-i-retrograde-high` | `pan2016` | 31338.1 | 31326.3 | 53.8 | 7 |
| `der-article-ii-prograde-high` | `battin1984` | 37.2 | 37.8 | 0.3 | 5472 |
| `der-article-ii-prograde-high` | `mcelreath2025` | 47.3 | 48.2 | 0.6 | 4315 |
| `der-article-ii-prograde-high` | `arora2013` | 53.8 | 54.7 | 0.5 | 3766 |
| `der-article-ii-prograde-high` | `izzo2015` | 62.0 | 63.2 | 1.0 | 3307 |
| `der-article-ii-prograde-high` | `gooding1990` | 85.0 | 86.9 | 1.3 | 2384 |
| `der-article-ii-prograde-high` | `thorne2004` | 88.6 | 91.7 | 3.0 | 2323 |
| `der-article-ii-prograde-high` | `jiang2016` | 108.3 | 110.3 | 2.0 | 1886 |
| `der-article-ii-prograde-high` | `der2011` | 119.7 | 121.9 | 2.5 | 1704 |
| `der-article-ii-prograde-high` | `avanzini2008` | 128.2 | 133.9 | 3.8 | 1584 |
| `der-article-ii-prograde-high` | `negrete2024` | 130.0 | 132.2 | 1.2 | 1548 |
| `der-article-ii-prograde-high` | `delatorre2018` | 311.7 | 315.3 | 9.1 | 647 |
| `der-article-ii-prograde-high` | `pan2016` | 36317.0 | 36430.4 | 149.6 | 6 |
| `der-article-ii-retrograde-high` | `battin1984` | 35.4 | 35.9 | 0.3 | 5747 |
| `der-article-ii-retrograde-high` | `mcelreath2025` | 47.1 | 48.0 | 0.6 | 4365 |
| `der-article-ii-retrograde-high` | `izzo2015` | 61.1 | 62.3 | 1.1 | 3321 |
| `der-article-ii-retrograde-high` | `arora2013` | 63.3 | 64.2 | 0.8 | 3203 |
| `der-article-ii-retrograde-high` | `gooding1990` | 85.2 | 86.9 | 1.4 | 2380 |
| `der-article-ii-retrograde-high` | `thorne2004` | 89.6 | 91.6 | 1.8 | 2290 |
| `der-article-ii-retrograde-high` | `jiang2016` | 109.4 | 111.1 | 2.6 | 1881 |
| `der-article-ii-retrograde-high` | `avanzini2008` | 122.2 | 125.4 | 3.5 | 1661 |
| `der-article-ii-retrograde-high` | `negrete2024` | 129.6 | 131.5 | 1.3 | 1553 |
| `der-article-ii-retrograde-high` | `der2011` | 139.5 | 141.8 | 3.5 | 1466 |
| `der-article-ii-retrograde-high` | `delatorre2018` | 301.7 | 304.5 | 9.0 | 677 |
| `der-article-ii-retrograde-high` | `pan2016` | 36091.1 | 36143.1 | 237.0 | 6 |

### zero-rev-hyperbolic

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `gmat-prograde` | `battin1984` | 39.6 | 40.8 | 0.5 | 5140 |
| `gmat-prograde` | `mcelreath2025` | 63.1 | 65.4 | 1.3 | 3257 |
| `gmat-prograde` | `arora2013` | 63.4 | 65.5 | 1.4 | 3297 |
| `gmat-prograde` | `izzo2015` | 84.1 | 87.3 | 2.1 | 2438 |
| `gmat-prograde` | `gooding1990` | 108.3 | 115.4 | 5.2 | 1881 |
| `gmat-prograde` | `thorne2004` | 124.3 | 128.8 | 5.7 | 1654 |
| `gmat-prograde` | `negrete2024` | 135.8 | 139.2 | 5.1 | 1486 |
| `gmat-prograde` | `jiang2016` | 159.3 | 164.2 | 7.2 | 1291 |
| `gmat-prograde` | `avanzini2008` | 163.2 | 168.3 | 6.6 | 1259 |
| `gmat-prograde` | `der2011` | 177.8 | 199.2 | 31.2 | 1187 |
| `gmat-prograde` | `delatorre2018` | 365.1 | 370.9 | 15.5 | 559 |
| `gmat-prograde` | `pan2016` | 42659.1 | 42597.4 | 187.8 | 5 |
| `gmat-retrograde` | `battin1984` | 39.2 | 40.4 | 0.5 | 5155 |
| `gmat-retrograde` | `mcelreath2025` | 63.1 | 65.5 | 1.5 | 3262 |
| `gmat-retrograde` | `arora2013` | 74.2 | 76.4 | 2.0 | 2849 |
| `gmat-retrograde` | `izzo2015` | 83.9 | 96.9 | 3.7 | 2457 |
| `gmat-retrograde` | `gooding1990` | 107.6 | 112.0 | 5.0 | 1890 |
| `gmat-retrograde` | `thorne2004` | 127.4 | 131.6 | 5.4 | 1609 |
| `gmat-retrograde` | `negrete2024` | 135.6 | 139.3 | 5.5 | 1491 |
| `gmat-retrograde` | `der2011` | 141.2 | 145.7 | 6.8 | 1457 |
| `gmat-retrograde` | `jiang2016` | 157.5 | 162.1 | 6.8 | 1288 |
| `gmat-retrograde` | `avanzini2008` | 164.7 | 170.1 | 7.2 | 1256 |
| `gmat-retrograde` | `delatorre2018` | 365.1 | 370.3 | 15.7 | 560 |
| `gmat-retrograde` | `pan2016` | 43008.7 | 42973.1 | 197.4 | 5 |

### zero-rev-nominal

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `battin-book` | `gauss1809` | 32.7 | 33.7 | 0.5 | 6192 |
| `battin-book` | `battin1984` | 35.7 | 36.7 | 0.4 | 5659 |
| `battin-book` | `mcelreath2025` | 60.0 | 61.9 | 1.3 | 3431 |
| `battin-book` | `arora2013` | 76.0 | 78.1 | 1.4 | 2747 |
| `battin-book` | `izzo2015` | 83.1 | 85.5 | 1.5 | 2470 |
| `battin-book` | `jiang2016` | 102.4 | 105.2 | 2.4 | 2014 |
| `battin-book` | `gooding1990` | 107.7 | 111.3 | 2.7 | 1893 |
| `battin-book` | `thorne2004` | 120.6 | 125.1 | 4.9 | 1707 |
| `battin-book` | `negrete2024` | 128.6 | 131.8 | 4.1 | 1562 |
| `battin-book` | `vallado2013` | 130.0 | 132.9 | 3.9 | 1548 |
| `battin-book` | `der2011` | 150.0 | 154.5 | 6.2 | 1364 |
| `battin-book` | `avanzini2008` | 153.7 | 160.0 | 5.1 | 1330 |
| `battin-book` | `delatorre2018` | 308.8 | 312.5 | 13.8 | 659 |
| `battin-book` | `pan2016` | 35382.9 | 35364.0 | 101.6 | 6 |
| `vallado-book` | `battin1984` | 37.5 | 38.5 | 0.5 | 5406 |
| `vallado-book` | `gauss1809` | 51.2 | 52.5 | 0.8 | 3951 |
| `vallado-book` | `mcelreath2025` | 60.8 | 63.2 | 1.3 | 3377 |
| `vallado-book` | `vallado2013` | 65.4 | 66.8 | 0.7 | 3076 |
| `vallado-book` | `arora2013` | 76.7 | 78.8 | 1.3 | 2658 |
| `vallado-book` | `izzo2015` | 85.0 | 88.1 | 2.1 | 2420 |
| `vallado-book` | `gooding1990` | 106.9 | 110.3 | 2.9 | 1893 |
| `vallado-book` | `pan2016` | 113.1 | 116.4 | 2.7 | 1791 |
| `vallado-book` | `thorne2004` | 118.5 | 122.5 | 4.3 | 1732 |
| `vallado-book` | `negrete2024` | 130.4 | 133.8 | 4.4 | 1548 |
| `vallado-book` | `avanzini2008` | 162.1 | 167.4 | 5.6 | 1276 |
| `vallado-book` | `der2011` | 187.8 | 192.9 | 8.6 | 1098 |
| `vallado-book` | `delatorre2018` | 344.7 | 348.3 | 13.8 | 592 |
| `vallado-book` | `jiang2016` | 394.6 | 399.3 | 15.0 | 512 |

<!-- performance-comparison:end -->
