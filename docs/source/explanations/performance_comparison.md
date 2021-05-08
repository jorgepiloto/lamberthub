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

# Solvers performance comparison

One of the goals of `lamberthub` is to provide different performance comparison
tools for new authors addressing the Lambert's problem. At the moment, the
following utilities are available:

* **Iterations performance plotter:** generates a contour map showing the number of
  iterations for a particular combination of the transfer angle and the
  non-dimensional time of flight.

In this page, all these tools are presented together with the results output by
them.


## Iterations contour maps

A first performance comparison between the current implemented solvers can be
performed via the `IterationsPlotter` class. This class provides all the
necessary methods to output a figure showing the relation between the transfer
angle, the non-dimensional time of flight and the number of iterations.

The following code shows how to create a collection of these figures for a given
set of solvers. Taking advantage of the `ALL_SOLVERS` macro provided by
`lamberthub`, it is possible to build a graphical representation for all the
algoritms currently available.

```{code-cell}
import matplotlib.pyplot as plt

from lamberthub import ALL_SOLVERS
from lamberthub.plotting import IterationsPlotter

def plot_iterations_performance(SOLVERS_LIST):
    """ Plots iterations performance for each one of the solvers from given set """

    for solver in SOLVERS_LIST:
        fig, ax = plt.subplots()
        iter_plotter = IterationsPlotter(ax=ax, fig=fig)
        iter_plotter.plot_performance(solver, maxiter=10, cmap=plt.get_cmap("YlOrRd"))
    plt.show()

plot_iterations_performance(ALL_SOLVERS)
```
