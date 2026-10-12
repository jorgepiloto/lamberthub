# lamberthub

<p align="center">
  <a href="https://github.com/jorgepiloto/lamberthub">
    <img width="350" src="https://raw.githubusercontent.com/jorgepiloto/lamberthub/main/doc/source/_static/logo.png" alt="lamberthub">
  </a>
</p>

<!-- vale off -->
[![Python](https://img.shields.io/pypi/pyversions/lamberthub?logo=pypi)](https://pypi.org/project/lamberthub/)
[![PyPI](https://img.shields.io/pypi/v/lamberthub.svg?logo=python&logoColor=white)](https://pypi.org/project/lamberthub/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CI](https://github.com/jorgepiloto/lamberthub/actions/workflows/pull-request.yml/badge.svg?branch=main)](https://github.com/jorgepiloto/lamberthub/actions/workflows/pull-request.yml)
[![Coverage](https://codecov.io/gh/jorgepiloto/lamberthub/branch/main/graph/badge.svg?token=3BY2J5AB8D)](https://codecov.io/gh/jorgepiloto/lamberthub)
[![DOI](https://zenodo.org/badge/364482782.svg)](https://zenodo.org/badge/latestdoi/364482782)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
<!-- vale on -->

## About

`lamberthub` provides a common Python interface to Lambert's problem
algorithms. These algorithms determine an orbit from two position vectors and
a time of flight. The package supports direct and multi-revolution transfers
through solvers named after their authors and publication years.

## Installation

Review the [prerequisites] before installing `lamberthub`. Install the latest
release from PyPI:

```console
python -m pip install lamberthub
```

The [installation guidelines] also describe how to install the development
version from GitHub.

## Documentation

The [official documentation] contains the following sections:

- [Getting Started] explains what the package provides, its prerequisites, and
  how to install it.
- [User Guide] documents the shared solver interface, parameters, and return
  values.
- [API Reference] is generated from the package source and NumPy-style
  documentation strings.
- [Examples] provides a complete direct-transfer calculation.
- [Changelog] records user-visible changes.

## Troubleshooting

For troubleshooting or bug reports, [open an issue]. Include the `lamberthub`
version, Python version, operating system, a minimal reproducer, and the full
error message.

## License

`lamberthub` is licensed under the GNU General Public License v3.0. See the
[LICENSE] file for the full text.

## Changelog

The [project changelog] records notable changes for each release.

[prerequisites]: https://jorgemartinez.space/projects/lamberthub/version/stable/getting-started/prerequisites.html
[installation guidelines]: https://jorgemartinez.space/projects/lamberthub/version/stable/getting-started/installation.html
[official documentation]: https://jorgemartinez.space/projects/lamberthub/version/stable/
[Getting Started]: https://jorgemartinez.space/projects/lamberthub/version/stable/getting-started.html
[User Guide]: https://jorgemartinez.space/projects/lamberthub/version/stable/user-guide.html
[API Reference]: https://jorgemartinez.space/projects/lamberthub/version/stable/api-reference.html
[Examples]: https://jorgemartinez.space/projects/lamberthub/version/stable/examples.html
[Changelog]: https://jorgemartinez.space/projects/lamberthub/version/stable/changelog.html
[open an issue]: https://github.com/jorgepiloto/lamberthub/issues/new
[LICENSE]: https://github.com/jorgepiloto/lamberthub/blob/main/LICENSE
[project changelog]: https://jorgemartinez.space/projects/lamberthub/version/stable/changelog.html
