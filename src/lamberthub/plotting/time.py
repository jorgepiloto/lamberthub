""" Holds plotting utilities related with required number of iterations """


import numpy as np
from cmaps import sunshine_9lev
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

from lamberthub.plotting._base import TauThetaPlotter, _vec_measure_performance
from lamberthub.utils.misc import get_solver_name


class TPIPlotter(TauThetaPlotter):
    """ Time per iteration """

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

    def _get_tpi(self, solver, theta_span, tau_span):
        """ Computes the time per iteration for the whole meshgrid """

        # Compute the meshgrid holding the number of required iterations
        _, TPI, _ = _vec_measure_performance(
            solver, theta_span[np.newaxis, :], tau_span[:, np.newaxis]
        )
        return TPI.astype("float")

    def plot_performance(
        self, solver, N_samples=10, maxtpi=200, step=20, cmap=sunshine_9lev
    ):
        """
        Returns a graphical representation on the time per iteration performance
        for a particular solver.

        Parameters
        ----------
        solver: function
            The solver who's performance is to be assessed.
        N_samples: int
            Number of samples to be computed. The higher, the less spurious
            values in the result.
        maxtpi: float
            The maximum value for the time per iterations (in microseconds).
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

        # Compute the time per iteration in the form of meshgrid
        TPI_SUM = np.zeros(THETA.shape)

        for _ in range(N_samples):

            # Compute the meshgrid holding the number of required iterations
            TPI = self._get_tpi(solver, theta_span, tau_span)
            # Add the result to the total amount
            TPI_SUM += TPI

        # The time per iteration is computed and stored in microseconds per iteration
        TPI = TPI_SUM / N_samples * 1e6

        # Store the mean time in microseconds per iteration
        MEAN_TPI = np.mean(TPI)

        # Prepare the levels for the contour
        levels = MaxNLocator(nbins=11).tick_values(0, maxtpi)
        bd_norm = BoundaryNorm(levels, ncolors=cmap.N)

        # Generate discrete contour map
        self.collection = self.ax.pcolor(
            THETA,
            TAU,
            TPI[:-1, :-1],  # For pcolor, the last rows need to be removed
            cmap=cmap,
            vmin=1e-12,
            edgecolors="k",
            linewidths=1,
            norm=bd_norm,
        )
        self.collection.cmap.set_under("black")

        # Draw the colorbar with proper diemensions and label
        label = f"Time per iteration\nAverage {MEAN_TPI:.2f} " + r"$\mu$" + "s / iter"
        self._draw_colorbar(maxtpi, step, label, cmap, "black")

        # Finally, draw the ticks, labels and title
        self._draw_ticks()
        self._draw_labels()
        self.ax.set_title(get_solver_name(solver))

        return self.ax


class TTCPlotter(TauThetaPlotter):
    """ Total time computation """

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

    def _get_ttc(self, solver, theta_span, tau_span):
        """ Computes the total time of computation for the whole meshgrid """

        # Compute the meshgrid holding the number of required iterations
        _, _, TTC = _vec_measure_performance(
            solver, theta_span[np.newaxis, :], tau_span[:, np.newaxis]
        )

        return TTC.astype("float")

    def plot_performance(
        self, solver, N_samples=10, maxttc=1000, step=100, cmap=sunshine_9lev
    ):
        """
        Returns a graphical representation on the time per iteration performance
        for a particular solver.

        Parameters
        ----------
        solver: function
            The solver who's performance is to be assessed.
        N_samples: int
            Number of samples to be computed. The higher, the less spurious
            values in the result.
        maxttc: float
            The maximum value for the total time of computation (in microseconds).
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

        # Compute the time per iteration in the form of meshgrid
        TTC_SUM = np.zeros(THETA.shape)

        for _ in range(N_samples):

            # Compute the meshgrid holding the number of required iterations
            TTC = self._get_ttc(solver, theta_span, tau_span)
            # Add the result to the total amount
            TTC_SUM += TTC

        # The time per iteration is computed and stored in microseconds per iteration
        TTC = TTC_SUM / N_samples * 1e6

        # Store the mean time in microseconds per iteration
        MEAN_TTC = np.mean(TTC)

        # Prepare the levels for the contour
        levels = MaxNLocator(nbins=11).tick_values(0, maxttc)
        bd_norm = BoundaryNorm(levels, ncolors=cmap.N)

        # Generate discrete contour map
        self.collection = self.ax.pcolor(
            THETA,
            TAU,
            TTC[:-1, :-1],  # For pcolor, the last rows need to be removed
            cmap=cmap,
            vmin=1,
            edgecolors="k",
            linewidths=1,
            norm=bd_norm,
        )
        self.collection.cmap.set_under("black")

        # Draw the colorbar with proper diemensions and label
        label = f"Total computation time\nAverage {MEAN_TTC:.2f} " + r"$\mu$" + "s"
        self._draw_colorbar(maxttc, step, label, cmap, "black")

        # Finally, draw the ticks, labels and title
        self._draw_ticks()
        self._draw_labels()
        self.ax.set_title(get_solver_name(solver))

        return self.ax
