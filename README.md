# lamberthub: a hub of Lambert's problem solvers

<img align="left" width=350px src="https://github.com/jorgepiloto/lamberthub/raw/main/doc/source/_static/lamberts_problem_geometry.png"/>

[![Python](https://img.shields.io/pypi/pyversions/lamberthub?logo=pypi)](https://pypi.org/project/lamberthub/)
[![PyPI](https://img.shields.io/pypi/v/lamberthub.svg?logo=python&logoColor=white)](https://pypi.org/project/lamberthub/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GH-CI](https://github.com/jorgepiloto/lamberthub/actions/workflows/ci_cd.yml/badge.svg)](https://github.com/jorgepiloto/lamberthub/actions/workflows/ci_cd.yml)
[![codecov](https://codecov.io/gh/jorgepiloto/lamberthub/branch/main/graph/badge.svg?token=3BY2J5AB8D)](https://codecov.io/gh/jorgepiloto/lamberthub)
[![DOI](https://zenodo.org/badge/364482782.svg)](https://zenodo.org/badge/latestdoi/364482782)

A collection of Lambert's problem solvers implemented using modern Python.

**Install the latest stable release by running:**

```bash
python -m pip install lamberthub
```

## Available solvers

Verify available solvers by running:

```python
from lamberthub import ALL_SOLVERS


print([solver.__name__ for solver in ALL_SOLVERS])
```

```pycon
>>> [
    'gauss1809', 'battin1984', 'gooding1990', 'avanzini2008',  'arora2013', 
    'vallado2013', 'izzo2015'
]
```

## Using a solver

Any Lambert's problem algorithm implemented in `lamberthub` is a Python function
which accepts the following parameters:

```python
from lamberthub import authorYYYY


v1, v2 = authorYYYY(
    mu, r1, r2, tof, M, prograde=True, low_path=True, maxiter=35, 
    atol=1e-5, rtol=1e-7, full_output=False
)
```

where `author` is the name of the author which developed the solver and `YYYY`
the year of publication. Any of the solvers hosted by the `ALL_SOLVERS` macro
can be used.

**Parameters**
- `mu`: the gravitational parameter, that is the mass of the attracting body
- times the gravitational constant.
- `r1`: initial position vector.
- `r2`: final position vector.
- `tof`: time of flight between initial and final vectors.

**Additional parameters**
- `M`: the number of revolutions. If zero (default), direct transfer is assumed.
- `prograde`: this parameter controls the inclination of the final orbit. If set
- to `True`, the transfer will have an inclination between 0 and 90 degrees
- while if `False` inclinations between 90 and 180 are provided.
- `low_path`: selects the type of path when more than two solutions are available.
- There is no actual advantage on one or another solution, unless you have
- particular constrains on your mission.
- `maxiter`: maximum number of iterations allowed when computing the solution.
- `atol`: absolute tolerance for the iterative method.
- `rtol`: relative tolerance for the iterative method.
- `full_output`: if `True`, it returns additional information such us the number
  of iterations. 

**Returns**
* `v1`: initial velocity vector.
* `v2`: final velocity vector.

**Additional returns**
* `numiter`: number of iterations. Only if `full_output` has been set to `True`.
* `tpi`: time per iteration. Only if `full_output` has been set to `True`.
