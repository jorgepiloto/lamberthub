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
| `avanzini2008`| G. Avanzini, “A simple lambert algorithm,” *Journal of Guidance, Control, and Dynamics*, vol. 31, no. 6, pp. 1587–1594, 2008.                          |
| `arora2013`   | N. Arora and R. P. Russell, “A fast and robust multiple revolution lambert algorithm using a cosine transformation,” Paper AAS, vol. 13, p. 728, 2013.  |
| `vallado2013` | D. A. Vallado, *Fundamentals of astrodynamics and applications*. Springer Science & Business Media, 2013, vol. 12.                                       |
| `izzo2015`    | D. Izzo, “Revisiting lambert’s problem,” *Celestial Mechanics and Dynamical Astronomy*, vol. 121, no. 1, pp. 1–15, 2015.                                |

<!-- vale on -->

## Using a solver

Any Lambert's problem algorithm implemented in `lamberthub` is a Python function
which accepts the following parameters:

```python
from lamberthub import authorYYYY


v1, v2 = authorYYYY(
    mu, r1, r2, tof, M=0, prograde=True, low_path=True,  # Type of solution
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
| `prograde`    | `bool`    | Controls the inclination of the final orbit. If `True`, inclination between 0 and 90 degrees. If `False`, inclination between 90 and 180 degrees. |
| `low_path`    | `bool`    | Selects the type of path when more than two solutions are available. No specific advantage unless there are mission constraints. |
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
possible solution, so the `low_path` flag does not play any role in this
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

v1, v2 = gooding1990(mu_sun, r1, r2, tof, M=0, prograde=True)
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

