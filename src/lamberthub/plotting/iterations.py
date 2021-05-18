""" Holds plotting utilities related with required number of iterations """

import numpy as np
from cmaps import sunshine_9lev
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

from lamberthub.plotting._base import TauThetaPlotter, _vec_measure_performance
from lamberthub.utils.misc import get_solver_name


class IterationsPlotter(TauThetaPlotter):
    """
    A class for plotting solvers number of iterations as functions of the
    transfer angle and the non-dimensional time of flight.
    """

    def __init__(self, ax=None, fig=None):
        """
        Initializes any instance of the plotter.

        Parameters
        ----------
        ax: matplotlib.Axes
            The axes in which the lines will be drawn.
        fig: matplotlib.Figure
            The figure instance for the plot.

        """

        # Initialize the mother class
        super().__init__(ax, fig)

    def _get_iterations(self, solver, theta_span, tau_span):
        """Returns a meshgrid holding the number of iterations for each one of
        the different evaluated points."""

        # Compute the meshgrid holding the number of required iterations
        NN_ITER, _, _ = _vec_measure_performance(
            solver, theta_span[np.newaxis, :], tau_span[:, np.newaxis]
        )
        return NN_ITER.astype("int")

    def plot_performance(self, solver, maxiter=10, step=1, cmap=sunshine_9lev):
        """
        Returns a graphical representation on the iteration performance for a
        particular solver.

        Parameters
        ----------
        solver: function
            The solver who's performance is to be assessed.
        maxiter: int
            The maximum number of iterations.
        step: float
            Step for drawing the colorbar ticks.
        cmap: matplotlib.colors.Colormap
            The map for colouring the grid.

        Notes
        -----
        The method does not accept neither a transfer angle nor time of flight
        span vectors. This is to force a common basis when addressing the
        performance from the number of iterations point of view. Only figure
        customization options are valid as input parameters.

        """

        # Compute the span and the associated meshgrid
        theta_span, tau_span = self._get_spans()
        THETA, TAU = np.meshgrid(theta_span, tau_span)

        # Solve for the number of iterations
        NN_ITER = self._get_iterations(solver, theta_span, tau_span)
        MEAN_NN_ITER = np.mean(NN_ITER)

        # Prepare the levels for the contour
        levels = MaxNLocator(nbins=11).tick_values(0, maxiter)
        bd_norm = BoundaryNorm(levels, ncolors=cmap.N)

        # Generate meshgrid
        self.collection = self.ax.pcolor(
            THETA,
            TAU,
            NN_ITER[:-1, :-1],  # For pcolor, the last rows need to be removed
            cmap=cmap,
            vmin=1,
            edgecolors="k",
            linewidths=1,
            norm=bd_norm,
        )
        self.collection.cmap.set_under("black")

        # Draw the colorbar with proper diemensions and label
        label = f"Number of iterations\nAverage {MEAN_NN_ITER:.2f} iter"
        self._draw_colorbar(maxiter, step, label, cmap, "white")

        # Finally, draw the ticks, labels and title
        self._draw_ticks()
        self._draw_labels()
        self.ax.set_title(get_solver_name(solver))

        return self.ax
