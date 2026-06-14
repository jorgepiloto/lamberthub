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

<!-- performance-comparison:start -->

_This section is generated from CI benchmark artifacts._

- Generated: 2026-06-14 06:39 UTC
- Commit: `5fff8c99d908`
- Environment: Linux, Python 3.12.13, AMD EPYC 9V74 80-Core Processor

Times are reported in microseconds. Lower values are better.

### multi-rev-hard

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `m1-prograde-high` | `mcelreath2025` | 46.8 | 47.8 | 0.5 | 4332 |
| `m1-prograde-high` | `izzo2015` | 61.0 | 62.2 | 1.0 | 3351 |
| `m1-prograde-high` | `gooding1990` | 102.4 | 104.7 | 2.0 | 1992 |
| `m1-prograde-high` | `negrete2024` | 194.7 | 197.4 | 3.8 | 1034 |
| `m1-prograde-high` | `arora2013` | 213.2 | 215.6 | 7.1 | 969 |
| `m1-prograde-high` | `der2011` | 314.0 | 320.9 | 10.1 | 649 |
| `m1-prograde-high` | `delatorre2018` | 651.5 | 651.9 | 9.9 | 315 |
| `m1-prograde-low` | `mcelreath2025` | 47.3 | 48.1 | 0.5 | 4292 |
| `m1-prograde-low` | `izzo2015` | 61.5 | 62.6 | 0.9 | 3319 |
| `m1-prograde-low` | `gooding1990` | 102.5 | 104.6 | 1.9 | 1981 |
| `m1-prograde-low` | `negrete2024` | 194.9 | 197.1 | 2.8 | 1033 |
| `m1-prograde-low` | `arora2013` | 233.1 | 238.1 | 9.9 | 872 |
| `m1-prograde-low` | `der2011` | 421.0 | 423.1 | 10.5 | 486 |
| `m1-prograde-low` | `delatorre2018` | 704.6 | 707.3 | 7.7 | 291 |
| `m2-prograde-high` | `mcelreath2025` | 47.5 | 48.4 | 0.5 | 4300 |
| `m2-prograde-high` | `izzo2015` | 61.4 | 62.8 | 1.0 | 3343 |
| `m2-prograde-high` | `gooding1990` | 104.0 | 106.3 | 2.0 | 1954 |
| `m2-prograde-high` | `negrete2024` | 196.1 | 202.7 | 8.4 | 1027 |
| `m2-prograde-high` | `arora2013` | 215.3 | 221.5 | 10.9 | 957 |
| `m2-prograde-high` | `der2011` | 267.9 | 271.3 | 9.9 | 758 |
| `m2-prograde-high` | `delatorre2018` | 650.2 | 651.0 | 9.4 | 316 |
| `m2-prograde-low` | `mcelreath2025` | 47.2 | 48.4 | 0.5 | 4315 |
| `m2-prograde-low` | `izzo2015` | 61.9 | 63.0 | 0.9 | 3303 |
| `m2-prograde-low` | `gooding1990` | 103.8 | 106.5 | 2.0 | 1960 |
| `m2-prograde-low` | `negrete2024` | 196.3 | 201.4 | 7.2 | 1027 |
| `m2-prograde-low` | `arora2013` | 229.3 | 236.1 | 10.9 | 895 |
| `m2-prograde-low` | `der2011` | 267.9 | 271.4 | 9.6 | 764 |
| `m2-prograde-low` | `delatorre2018` | 680.2 | 681.3 | 9.4 | 299 |

### multi-rev-m1

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `prograde-high` | `mcelreath2025` | 47.5 | 48.7 | 0.7 | 4324 |
| `prograde-high` | `izzo2015` | 61.4 | 62.7 | 1.0 | 3319 |
| `prograde-high` | `gooding1990` | 104.3 | 106.7 | 2.3 | 1950 |
| `prograde-high` | `arora2013` | 158.7 | 161.3 | 6.2 | 1298 |
| `prograde-high` | `negrete2024` | 197.8 | 200.3 | 2.4 | 1018 |
| `prograde-high` | `der2011` | 273.8 | 277.8 | 10.3 | 747 |
| `prograde-high` | `delatorre2018` | 643.2 | 646.1 | 10.8 | 320 |
| `prograde-low` | `mcelreath2025` | 47.5 | 48.4 | 0.6 | 4269 |
| `prograde-low` | `izzo2015` | 61.9 | 63.0 | 0.9 | 3309 |
| `prograde-low` | `gooding1990` | 104.1 | 106.4 | 2.2 | 1940 |
| `prograde-low` | `arora2013` | 160.5 | 163.9 | 5.1 | 1273 |
| `prograde-low` | `negrete2024` | 197.8 | 200.4 | 3.0 | 1019 |
| `prograde-low` | `der2011` | 251.6 | 271.4 | 14.3 | 812 |
| `prograde-low` | `delatorre2018` | 618.8 | 621.7 | 11.4 | 331 |
| `retrograde-high` | `mcelreath2025` | 47.0 | 47.8 | 0.7 | 4291 |
| `retrograde-high` | `izzo2015` | 61.9 | 63.1 | 1.0 | 3282 |
| `retrograde-high` | `gooding1990` | 104.0 | 106.4 | 2.2 | 1960 |
| `retrograde-high` | `arora2013` | 146.6 | 149.5 | 3.8 | 1398 |
| `retrograde-high` | `negrete2024` | 195.3 | 197.6 | 2.9 | 1031 |
| `retrograde-high` | `der2011` | 310.6 | 313.7 | 10.1 | 658 |
| `retrograde-high` | `delatorre2018` | 630.6 | 632.3 | 9.9 | 324 |
| `retrograde-low` | `mcelreath2025` | 47.7 | 48.6 | 0.6 | 4287 |
| `retrograde-low` | `izzo2015` | 61.9 | 62.9 | 1.0 | 3284 |
| `retrograde-low` | `gooding1990` | 103.4 | 105.7 | 2.3 | 1969 |
| `retrograde-low` | `arora2013` | 149.3 | 152.0 | 3.8 | 1375 |
| `retrograde-low` | `negrete2024` | 195.4 | 197.7 | 2.6 | 1030 |
| `retrograde-low` | `der2011` | 290.8 | 294.4 | 10.4 | 706 |
| `retrograde-low` | `delatorre2018` | 616.2 | 617.8 | 9.5 | 331 |

### robust-zero-rev

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `near-tangent-prograde-high` | `mcelreath2025` | 48.5 | 49.6 | 0.6 | 4241 |
| `near-tangent-prograde-high` | `arora2013` | 54.1 | 55.0 | 0.8 | 3784 |
| `near-tangent-prograde-high` | `izzo2015` | 61.4 | 62.5 | 1.0 | 3334 |
| `near-tangent-prograde-high` | `gooding1990` | 85.2 | 87.2 | 1.8 | 2383 |
| `near-tangent-prograde-high` | `negrete2024` | 129.3 | 131.0 | 1.3 | 1558 |
| `near-tangent-prograde-high` | `der2011` | 252.3 | 255.1 | 10.1 | 816 |
| `near-tangent-prograde-high` | `delatorre2018` | 321.1 | 324.2 | 9.2 | 634 |

### zero-rev-general

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `curtiss-book` | `battin1984` | 33.4 | 34.1 | 0.4 | 6123 |
| `curtiss-book` | `mcelreath2025` | 46.9 | 47.8 | 0.6 | 4337 |
| `curtiss-book` | `arora2013` | 56.5 | 57.4 | 0.8 | 3627 |
| `curtiss-book` | `izzo2015` | 60.0 | 61.0 | 1.0 | 3422 |
| `curtiss-book` | `gooding1990` | 83.6 | 85.2 | 1.5 | 2422 |
| `curtiss-book` | `thorne2004` | 97.6 | 99.7 | 2.2 | 2086 |
| `curtiss-book` | `der2011` | 101.1 | 103.2 | 2.3 | 2016 |
| `curtiss-book` | `jiang2016` | 119.5 | 121.8 | 2.7 | 1698 |
| `curtiss-book` | `avanzini2008` | 120.0 | 123.5 | 3.3 | 1702 |
| `curtiss-book` | `negrete2024` | 127.4 | 128.9 | 1.2 | 1587 |
| `curtiss-book` | `delatorre2018` | 322.3 | 325.8 | 9.7 | 636 |
| `curtiss-book` | `pan2016` | 32582.3 | 32569.6 | 130.7 | 7 |
| `der-article-i-prograde-low` | `battin1984` | 36.4 | 37.1 | 0.4 | 5602 |
| `der-article-i-prograde-low` | `mcelreath2025` | 46.7 | 47.9 | 0.6 | 4373 |
| `der-article-i-prograde-low` | `arora2013` | 53.7 | 54.7 | 0.7 | 3789 |
| `der-article-i-prograde-low` | `izzo2015` | 60.1 | 61.2 | 1.0 | 3415 |
| `der-article-i-prograde-low` | `gooding1990` | 82.6 | 84.2 | 1.5 | 2474 |
| `der-article-i-prograde-low` | `thorne2004` | 106.4 | 108.4 | 2.4 | 1929 |
| `der-article-i-prograde-low` | `avanzini2008` | 126.0 | 129.1 | 3.5 | 1621 |
| `der-article-i-prograde-low` | `negrete2024` | 129.0 | 130.8 | 1.4 | 1563 |
| `der-article-i-prograde-low` | `der2011` | 137.5 | 139.5 | 3.6 | 1488 |
| `der-article-i-prograde-low` | `jiang2016` | 197.2 | 200.2 | 6.6 | 1038 |
| `der-article-i-prograde-low` | `delatorre2018` | 284.9 | 288.4 | 9.5 | 714 |
| `der-article-i-prograde-low` | `pan2016` | 31106.2 | 31144.7 | 154.0 | 7 |
| `der-article-i-retrograde-high` | `battin1984` | 35.3 | 35.9 | 0.4 | 5759 |
| `der-article-i-retrograde-high` | `mcelreath2025` | 46.8 | 47.8 | 0.6 | 4358 |
| `der-article-i-retrograde-high` | `arora2013` | 56.1 | 57.1 | 0.7 | 3616 |
| `der-article-i-retrograde-high` | `izzo2015` | 60.1 | 61.2 | 1.0 | 3374 |
| `der-article-i-retrograde-high` | `gooding1990` | 82.9 | 84.5 | 1.5 | 2433 |
| `der-article-i-retrograde-high` | `thorne2004` | 100.5 | 102.3 | 2.2 | 2038 |
| `der-article-i-retrograde-high` | `avanzini2008` | 122.9 | 125.7 | 3.1 | 1656 |
| `der-article-i-retrograde-high` | `negrete2024` | 128.2 | 130.0 | 1.4 | 1570 |
| `der-article-i-retrograde-high` | `der2011` | 176.1 | 178.9 | 5.8 | 1161 |
| `der-article-i-retrograde-high` | `jiang2016` | 248.1 | 250.9 | 8.9 | 826 |
| `der-article-i-retrograde-high` | `delatorre2018` | 316.6 | 319.4 | 9.1 | 637 |
| `der-article-i-retrograde-high` | `pan2016` | 30979.6 | 30990.9 | 106.4 | 7 |
| `der-article-ii-prograde-high` | `battin1984` | 36.9 | 37.5 | 0.4 | 5490 |
| `der-article-ii-prograde-high` | `mcelreath2025` | 47.1 | 48.1 | 0.6 | 4355 |
| `der-article-ii-prograde-high` | `arora2013` | 53.1 | 54.1 | 0.7 | 3804 |
| `der-article-ii-prograde-high` | `izzo2015` | 60.4 | 61.3 | 1.0 | 3371 |
| `der-article-ii-prograde-high` | `gooding1990` | 84.0 | 85.6 | 1.4 | 2415 |
| `der-article-ii-prograde-high` | `thorne2004` | 88.3 | 90.1 | 2.7 | 2316 |
| `der-article-ii-prograde-high` | `jiang2016` | 107.6 | 109.3 | 2.1 | 1889 |
| `der-article-ii-prograde-high` | `der2011` | 118.5 | 120.6 | 2.8 | 1721 |
| `der-article-ii-prograde-high` | `avanzini2008` | 126.4 | 129.8 | 3.5 | 1610 |
| `der-article-ii-prograde-high` | `negrete2024` | 129.7 | 133.0 | 1.9 | 1556 |
| `der-article-ii-prograde-high` | `delatorre2018` | 313.1 | 320.7 | 9.6 | 649 |
| `der-article-ii-prograde-high` | `pan2016` | 35787.9 | 35899.0 | 176.0 | 6 |
| `der-article-ii-retrograde-high` | `battin1984` | 35.4 | 35.9 | 0.4 | 5719 |
| `der-article-ii-retrograde-high` | `mcelreath2025` | 46.7 | 47.7 | 0.6 | 4350 |
| `der-article-ii-retrograde-high` | `izzo2015` | 60.6 | 61.6 | 0.9 | 3356 |
| `der-article-ii-retrograde-high` | `arora2013` | 63.3 | 64.3 | 0.8 | 3207 |
| `der-article-ii-retrograde-high` | `gooding1990` | 85.1 | 86.6 | 1.4 | 2370 |
| `der-article-ii-retrograde-high` | `thorne2004` | 89.8 | 92.0 | 1.9 | 2286 |
| `der-article-ii-retrograde-high` | `jiang2016` | 106.9 | 108.8 | 2.1 | 1909 |
| `der-article-ii-retrograde-high` | `avanzini2008` | 120.8 | 123.6 | 3.1 | 1682 |
| `der-article-ii-retrograde-high` | `negrete2024` | 129.5 | 131.3 | 1.4 | 1558 |
| `der-article-ii-retrograde-high` | `der2011` | 138.2 | 140.3 | 3.3 | 1469 |
| `der-article-ii-retrograde-high` | `delatorre2018` | 299.4 | 302.1 | 9.1 | 677 |
| `der-article-ii-retrograde-high` | `pan2016` | 36280.6 | 36334.5 | 206.3 | 6 |

### zero-rev-hyperbolic

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `gmat-prograde` | `battin1984` | 34.8 | 35.3 | 0.4 | 5850 |
| `gmat-prograde` | `mcelreath2025` | 47.9 | 49.0 | 0.7 | 4233 |
| `gmat-prograde` | `arora2013` | 52.7 | 53.6 | 0.7 | 3886 |
| `gmat-prograde` | `izzo2015` | 61.9 | 62.9 | 1.1 | 3312 |
| `gmat-prograde` | `gooding1990` | 85.2 | 86.6 | 1.8 | 2366 |
| `gmat-prograde` | `thorne2004` | 104.4 | 106.6 | 2.5 | 1990 |
| `gmat-prograde` | `avanzini2008` | 126.4 | 129.3 | 3.6 | 1619 |
| `gmat-prograde` | `jiang2016` | 131.6 | 133.9 | 3.0 | 1544 |
| `gmat-prograde` | `negrete2024` | 133.5 | 135.3 | 1.3 | 1513 |
| `gmat-prograde` | `der2011` | 141.8 | 143.7 | 3.2 | 1434 |
| `gmat-prograde` | `delatorre2018` | 334.4 | 336.8 | 9.4 | 616 |
| `gmat-prograde` | `pan2016` | 37101.9 | 37101.5 | 229.7 | 6 |
| `gmat-retrograde` | `battin1984` | 35.2 | 35.8 | 0.4 | 5743 |
| `gmat-retrograde` | `mcelreath2025` | 47.6 | 51.5 | 1.0 | 4238 |
| `gmat-retrograde` | `arora2013` | 62.0 | 62.8 | 0.8 | 3277 |
| `gmat-retrograde` | `izzo2015` | 63.0 | 64.0 | 1.0 | 3233 |
| `gmat-retrograde` | `gooding1990` | 84.3 | 86.1 | 1.4 | 2402 |
| `gmat-retrograde` | `thorne2004` | 106.7 | 108.5 | 2.4 | 1953 |
| `gmat-retrograde` | `der2011` | 115.6 | 117.9 | 3.1 | 1785 |
| `gmat-retrograde` | `avanzini2008` | 128.0 | 131.4 | 3.4 | 1593 |
| `gmat-retrograde` | `jiang2016` | 132.5 | 134.6 | 2.8 | 1541 |
| `gmat-retrograde` | `negrete2024` | 133.4 | 135.1 | 1.4 | 1509 |
| `gmat-retrograde` | `delatorre2018` | 329.2 | 333.4 | 9.7 | 617 |
| `gmat-retrograde` | `pan2016` | 37044.4 | 37066.7 | 413.7 | 6 |

### zero-rev-nominal

| Case | Solver | Median | Mean | IQR | Rounds |
|------|--------|-------:|-----:|----:|-------:|
| `battin-book` | `gauss1809` | 33.3 | 34.5 | 0.6 | 6074 |
| `battin-book` | `battin1984` | 36.6 | 37.8 | 0.5 | 5600 |
| `battin-book` | `mcelreath2025` | 61.2 | 63.7 | 1.2 | 3344 |
| `battin-book` | `arora2013` | 78.1 | 80.7 | 1.8 | 2661 |
| `battin-book` | `izzo2015` | 83.7 | 86.9 | 1.9 | 2492 |
| `battin-book` | `jiang2016` | 105.3 | 109.1 | 3.8 | 1960 |
| `battin-book` | `gooding1990` | 109.1 | 113.2 | 3.8 | 1854 |
| `battin-book` | `thorne2004` | 124.1 | 129.1 | 6.4 | 1644 |
| `battin-book` | `negrete2024` | 129.1 | 133.2 | 5.2 | 1567 |
| `battin-book` | `vallado2013` | 131.4 | 134.8 | 4.6 | 1531 |
| `battin-book` | `der2011` | 156.2 | 161.5 | 8.3 | 1317 |
| `battin-book` | `avanzini2008` | 156.9 | 162.5 | 6.5 | 1300 |
| `battin-book` | `delatorre2018` | 310.9 | 315.6 | 15.4 | 658 |
| `battin-book` | `pan2016` | 36477.1 | 36498.6 | 196.2 | 6 |
| `vallado-book` | `battin1984` | 38.2 | 39.5 | 0.5 | 5341 |
| `vallado-book` | `gauss1809` | 51.9 | 53.5 | 0.8 | 3915 |
| `vallado-book` | `mcelreath2025` | 62.1 | 66.3 | 1.7 | 3330 |
| `vallado-book` | `vallado2013` | 67.1 | 69.3 | 0.9 | 2988 |
| `vallado-book` | `arora2013` | 78.3 | 80.9 | 1.9 | 2610 |
| `vallado-book` | `izzo2015` | 84.0 | 87.2 | 2.1 | 2493 |
| `vallado-book` | `gooding1990` | 109.8 | 113.4 | 3.6 | 1854 |
| `vallado-book` | `pan2016` | 115.2 | 119.6 | 3.7 | 1777 |
| `vallado-book` | `thorne2004` | 122.8 | 127.2 | 5.0 | 1675 |
| `vallado-book` | `negrete2024` | 130.5 | 134.0 | 5.2 | 1545 |
| `vallado-book` | `avanzini2008` | 167.0 | 172.9 | 7.2 | 1236 |
| `vallado-book` | `der2011` | 193.4 | 198.6 | 10.6 | 1061 |
| `vallado-book` | `delatorre2018` | 356.2 | 360.9 | 15.9 | 583 |
| `vallado-book` | `jiang2016` | 401.0 | 406.7 | 16.2 | 516 |

<!-- performance-comparison:end -->
