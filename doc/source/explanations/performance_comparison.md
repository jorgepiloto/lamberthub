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

* **Time required per iteration:** shows a contour map in which the time per
  iteration is shown against the transfer angle and the non-dimensional time. It
  also computes the required average time. 

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
algorithms currently available.

```{code-cell}
import matplotlib.pyplot as plt

from lamberthub import ALL_SOLVERS
from lamberthub.plotting import IterationsPlotter

def plot_iterations_performance(SOLVERS_LIST):
    """ Plots iterations performance for each one of the solvers from given set """

    for solver in SOLVERS_LIST:
        fig, ax = plt.subplots()
        iter_plotter = IterationsPlotter(ax=ax, fig=fig)
        iter_plotter.plot_performance(solver, maxiter=10)
    plt.show()

plot_iterations_performance(ALL_SOLVERS)
```

## Time required per iteration

The number of iterations does not provide a full insight of the problem, as an
algorithm with low iterations number might require a lot of time per iteration
it finds a solution. Therefore, a plotter in which the time per iteration is
shown against a combination of transfer angle and non-dimensional time of
flight is required.

```{note}
Notice this performance tool is directly linked to machine specifications, the
number of processes being run in it, the current implementation and the number
of samples.
```

The following code snippet runs the time performance plotter for all the
available solvers. Results might vary in your local machine as these ones are
executed with a low number of samples in order to avoid `RuntimeErrors` in the
server computer which builds the documentation. 

```{code-cell}
import matplotlib.pyplot as plt

from lamberthub import ALL_SOLVERS
from lamberthub.plotting import TPIPlotter

def plot_time_performance(SOLVERS_LIST):
    """ Plots the contour maps and a bar chart for a given set of solvers """

    # Plots iterations performance for each one of the solvers from given set
    for solver in SOLVERS_LIST:
        fig_ctime, ax_ctime = plt.subplots()
        iter_plotter = TPIPlotter(ax=ax_ctime, fig=fig_ctime)
        iter_plotter.plot_performance(solver, N_samples=5)
    plt.show()

plot_time_performance(ALL_SOLVERS)
```

## Total time required

```{code-cell}
import matplotlib.pyplot as plt

from lamberthub import ALL_SOLVERS
from lamberthub.plotting import TTCPlotter

def plot_time_performance(SOLVERS_LIST):
    """ Plots the contour maps and a bar chart for a given set of solvers """

    # Plots iterations performance for each one of the solvers from given set
    for solver in SOLVERS_LIST:
        fig_ctime, ax_ctime = plt.subplots()
        iter_plotter = TTCPlotter(ax=ax_ctime, fig=fig_ctime)
        iter_plotter.plot_performance(solver, N_samples=5)
    plt.show()

plot_time_performance(ALL_SOLVERS)
```
