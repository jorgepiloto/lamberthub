# lamberthub: a hub of Lambert's problem solvers

<img align="left" width=350px src="https://github.com/jorgepiloto/lamberthub/blob/main/docs/source/_static/lamberts_problem_geometry.png"/>

![Python shield](https://img.shields.io/badge/%F0%9F%90%8D%20Python-%3E%3D%203.8-blue)
![CI shield](https://github.com/jorgepiloto/lamberthub/actions/workflows/ci_actions.yml/badge.svg?branch=main)
[![Documentation Status](https://readthedocs.org/projects/lamberthub/badge/?version=latest)](https://lamberthub.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/jorgepiloto/lamberthub/branch/main/graph/badge.svg?token=3BY2J5AB8D)](https://codecov.io/gh/jorgepiloto/lamberthub)
[![DOI](https://zenodo.org/badge/364482782.svg)](https://zenodo.org/badge/latestdoi/364482782)
[![PyPI version](https://badge.fury.io/py/lamberthub.svg)](https://badge.fury.io/py/lamberthub)

A collection of Lambert's problem solvers implemented using modern Python.

**Install the latest stable release by running:**

```bash
pip install lamberthub
```

Just in case you are interested on knowing what the problem is about, how to
solve it or which applications it has, please check the [official documentation
of the project](https://lamberthub.readthedocs.io/en/latest/). For further
information about software internals, refer to [API
documentation](https://lamberthub.readthedocs.io/en/latest/autoapi/lamberthub/index.html).


## Which solvers are available?

Once installed, you can start by checking which solvers `lamberthub` ships with
by running:

```python
from lamberthub import ALL_SOLVERS
print([solver.__name__ for solver in ALL_SOLVERS])
```

At this moment, the following algorithms are available:

```bash
>>> ['gauss1809', 'battin1984', 'gooding1990', 'avanzini2008', 'arora2013', 'vallado2013', 'izzo2015']
```

## How can I use a solver?

Any Lambert's problem algorithm implemented in `lamberthub` is a Python function
which accepts the following parameters:

```python
# Import a solver of your choice from the ones listed above
from lamberthub import authorYYYY
v1, v2 = authorYYYY(mu, r1, r2, tof, prograde=True, low_path=True, maxiter=35, atol=1e-5, rtol=1e-7, full_output=False)
```

where `author` is the name of the author which developed the solver and `YYYY`
the year of publication. Any of the solvers hosted by the `ALL_SOLVERS` macro
can be used.

**Parameters**
* `mu`: the gravitational parameter, that is the mass of the attracting body
  times the gravitational constant.
* `r1`: initial position vector.
* `r2`: final position vector.
* `tof`: time of flight between initial and final vectors.

**Additional parameters**
* `prograde`: this parameter controls the inclination of the final orbit. If set
  to `True`, the transfer will have an inclination between 0 and 90 degrees
  while if `False` inclinations between 90 and 180 are provided.
* `low_path`: selects the type of path when more than two solutions are available.
  There is no actual advantage on one or another solution, unless you have
  particular constrains on your mission.
* `maxiter`: maximum number of iterations allowed when computing the solution.
* `atol`: absolute tolerance for the iterative method.
* `rtol`: relative tolerance for the iterative method.
* `full_output`: if `True`, it returns additional information such us the number
  of iterations. 

**Returns**
* `v1`: initial velocity vector.
* `v2`: final velocity vector.

**Additional returns**
* `numiter`: number of iterations. Only if `full_output` has been set to `True`.
* `tpi`: time per iteration. Only if `full_output` has been set to `True`.

## Documentation and performance comparison tools

The [official lamberthub
documentation](https://lamberthub.readthedocs.io/en/latest/) contains different
how-to guides, explanations, tutorials and references related to the package.

If you are interested in the performance comparison tools provided, please refer
to [this documentation
chapter](https://lamberthub.readthedocs.io/en/latest/explanations/performance_comparison.html),
were you can find a brief tutorials on how to use those tools.
