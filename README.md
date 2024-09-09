# lamberthub: a hub of Lambert's problem solvers

<img align="left" width=350px src="https://github.com/jorgepiloto/lamberthub/raw/main/doc/source/_static/lamberts_problem_geometry.png"/>

A Python library designed to provide solutions to Lambert's problem, a
classical problem in astrodynamics that involves determining the orbit of a
spacecraft given two points in space and the time of flight between them. The
problem is essential for trajectory planning, particularly for interplanetary
missions.

Lamberthub implements multiple algorithms, each named after its author and
publication year, for solving different variations of Lambert's problem. These
algorithms can handle different types of orbits, including multi-revolution
paths and direct transfers.

<br>

[![Python](https://img.shields.io/pypi/pyversions/lamberthub?logo=pypi)](https://pypi.org/project/lamberthub/)
[![PyPI](https://img.shields.io/pypi/v/lamberthub.svg?logo=python&logoColor=white)](https://pypi.org/project/lamberthub/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GH-CI](https://github.com/jorgepiloto/lamberthub/actions/workflows/ci_cd.yml/badge.svg)](https://github.com/jorgepiloto/lamberthub/actions/workflows/ci_cd.yml)
[![codecov](https://codecov.io/gh/jorgepiloto/lamberthub/branch/main/graph/badge.svg?token=3BY2J5AB8D)](https://codecov.io/gh/jorgepiloto/lamberthub)
[![DOI](https://zenodo.org/badge/364482782.svg)](https://zenodo.org/badge/latestdoi/364482782)


## Installation

Install `lamberthub` by running:

```console
python -m pip install lamberthub
```

## Available solvers

| Algorithm     | Reference                                                                                                                                               |
|---------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|
| `gauss1809`   | C. F. Gauss, *Theoria motus corporum coelestium in sectionibus conicis solem ambientium*. 1809.                                                         |
| `battin1984`  | R. H. Battin and R. M. Vaughan, “An elegant lambert algorithm,” *Journal of Guidance, Control, and Dynamics*, vol. 7, no. 6, pp. 662–670, 1984.         |
| `gooding1990` | R. Gooding, “A procedure for the solution of lambert’s orbital boundary-value problem,” *Celestial Mechanics and Dynamical Astronomy*, vol. 48, no. 2, pp. 145–165, 1990. |
| `avanzini2008`| G. Avanzini, “A simple lambert algorithm,” *Journal of Guidance, Control, and Dynamics*, vol. 31, no. 6, pp. 1587–1594, 2008.                          |
| `arora2013`   | N. Arora and R. P. Russell, “A fast and robust multiple revolution lambert algorithm using a cosine transformation,” Paper AAS, vol. 13, p. 728, 2013.  |
| `vallado2013` | D. A. Vallado, *Fundamentals of astrodynamics and applications*. Springer Science & Business Media, 2013, vol. 12.                                       |
| `izzo2015`    | D. Izzo, “Revisiting lambert’s problem,” *Celestial Mechanics and Dynamical Astronomy*, vol. 121, no. 1, pp. 1–15, 2015.                                |

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
the year of publication. Any of the solvers hosted by the `ALL_SOLVERS` macro
can be used.

### Parameters and Returns

| Parameters    | Description |
|---------------|-------------|
| `mu`          | The gravitational parameter, i.e., mass of the attracting body times the gravitational constant. Equivalent to gravitational constant times the mass of the attractor body. |
| `r1`          | Initial position vector. |
| `r2`          | Final position vector. |
| `tof`         | Time of flight between initial and final vectors. |
| `M`           | The number of revolutions. If zero (default), direct transfer is assumed. |
| `prograde`    | Controls the inclination of the final orbit. If `True`, inclination between 0 and 90 degrees. If `False`, inclination between 90 and 180 degrees. |
| `low_path`    | Selects the type of path when more than two solutions are available. No specific advantage unless there are mission constraints. |
| `maxiter`     | Maximum number of iterations allowed when computing the solution. |
| `atol`        | Absolute tolerance for the iterative method. |
| `rtol`        | Relative tolerance for the iterative method. |
| `full_output` | If `True`, returns additional information such as the number of iterations. |

| Returns       | Description |
|---------------|-------------|
| `v1`          | Initial velocity vector. |
| `v2`          | Final velocity vector. |
| `numiter`     | Number of iterations (only if `full_output` is `True`). |
| `tpi`         | Time per iteration (only if `full_output` is `True`). |
