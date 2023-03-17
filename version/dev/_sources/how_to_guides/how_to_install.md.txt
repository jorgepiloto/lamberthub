---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# How to install

This page contains the guide for installing the `lamberthub` package. If you
experience any kind of problem during one of the steps shown in the following
lines, please open an new issue (or select similar ones) in the [issues
board](https://github.com/jorgepiloto/lamberthub/issues).

## Install from PyPI using pip

The installation process is similar to other python packages, meaning that you
only need to run:

```bash
pip install lamberthub
```

previous command installs the latest stable version of the library. Once
done, you can open the Python terminal and import the package and verify its
version by running:

```{code-cell} ipython3
import lamberthub
print(f"Current lamberthub version is {lamberthub.__version__}")
```
