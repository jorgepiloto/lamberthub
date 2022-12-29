# lamberthub: a hub of Lambert's problem solvers

<img align="left" width=350px src="https://github.com/jorgepiloto/lamberthub/raw/main/doc/source/_static/lamberts_problem_geometry.png"/>

![Python shield](https://img.shields.io/badge/Python-%3E%3D%203.8-blue)
[![GH-CI](https://github.com/jorgepiloto/lamberthub/actions/workflows/ci_cd.yml/badge.svg)](https://github.com/jorgepiloto/lamberthub/actions/workflows/ci_cd.yml)
[![codecov](https://codecov.io/gh/jorgepiloto/lamberthub/branch/main/graph/badge.svg?token=3BY2J5AB8D)](https://codecov.io/gh/jorgepiloto/lamberthub)
[![DOI](https://zenodo.org/badge/364482782.svg)](https://zenodo.org/badge/latestdoi/364482782)
[![PyPI](https://img.shields.io/pypi/v/lamberthub.svg?logo=python&logoColor=white)](https://pypi.org/project/lamberthub/)
[![Apache 2.0](https://img.shields.io/badge/License-Apache-lightgray.svg)](https://opensource.org/licenses/Apache-2.0)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat)](https://github.com/psf/black)

A collection of Lambert's problem solvers implemented using modern Python.

**Install the latest stable release by running:**

```bash
python -m pip install lamberthub
```

Just in case you are interested on knowing what the problem is about, how to
solve it or which applications it has, please check the [official lamberthub
documentation]. For further information about software internals, refer to [API
documentation].


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

## Input and output parameters

All solvers accept and return the same amount of parameters. The tables below
explain which input and output parameters are required and which others are
optional.

**Input parameters**
These parameters are always required to be passed to the desired solver.

| Input parameters | Description                                                                                            |
|:----------------:|:------------------------------------------------------------------------------------------------------:|
| ``mu``           | The gravitational parameter, that is the mass of the attracting body times the gravitational constant. |
| ``r1``           | The initial position vector.                                                                           |
| ``r2``           | The final position vector.                                                                             |
| ``tof``          | The time of flight between initial and final vectors.                                                  |

**Optional parameters**
These parameters are optional. If no value is provided for them, a default one is used, see their description.

| Optional input parameters |                                                                            Description                                                                            |
|:-------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|           ``M``           |                                                          The number of revolutions. Default value is 0.                                                           |
|       ``prograde``        | Default to ``True``, which assumes the transfer to have an orbit inclination between 0 and 90 degrees. If ``False``, inclination lies between 90 and 180 degrees. |
|       ``low_path``        |               Default to ``True``, which selects the low path arc when more than two solutions exist. If ``False``, the high path arc is selected.                |
|        ``maxiter``        |                                                  The maximum allowed number of iterations. Default value is 35.                                                   |
|         ``atol``          |                                              The absolute tolerance of the numerical method. Default value is 1e-7.                                               |
|         ``rtol``          |                                              The relative tolerance of the numerical method. Default value is 1e-5.                                               |
|      ``full_output``      |                                         If ``True``, it returns additional information such us the number of iterations.                                          |

**Output parameters**
These paremeters are always returned by any solver.

| Output parameters | Description                  |
|:-----------------:|:----------------------------:|
| ``v1``            | The initial velocity vector. |
| ``v2``            | The final velocity vector.   |

**Optional output parameters**
These parameters are only returned if the ``full_output`` input parameter has been set to ``True``.

| Optional output parameters |        Description        |
|:--------------------------:|:-------------------------:|
|        ``numiter``         | The number of iterations. |
|          ``tpi``           |  The time per iteration.  |

## Documentation and performance comparison tools

The [official lamberthub documentation] contains different how-to guides,
explanations, tutorials and references related to the package.

If you are interested in the performance comparison tools provided, please refer
to [performance comparison section]. Here you can find a brief tutorials on how
to use those tools.


<!-- Links and References -->

[official lamberthub documentation]: https://jorgemartinez.space/projects/lamberthub/index.html
[API documentation]: https://jorgemartinez.space/projects/lamberthub/autoapi/lamberthub/index.html
[performance comparison section]: https://jorgemartinez.space/projects/lamberthub/explanations/performance_comparison.html
