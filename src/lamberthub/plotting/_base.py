""" A collection of common routines for plotting ones """

import time

import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin


class TauThetaPlotter:
    """ A class for modelling a discrete grid contour plotter. """

    def __init__(self, ax=None, fig=None, Nres=50):
        """
        Initializes any instance of the plotter.

        Parameters
        ----------
        func: function
            The core function which provides the results to be shown.
        ax: matplotlib.Axes
            The axes in which the lines will be drawn.
        fig: matplotlib.Figure
            The figure instance for the plot.
        Nres: int
            Number of total elements

        """

        # Check if axes are available
        if ax is None:
            _, ax = plt.subplots()

        # Check if figure available
        if fig is None:
            fig, _ = plt.subplots()

        # Assign previous figure and axes. Impose equal aspect ratio.
        self.ax, self.fig = ax, fig
        self.ax.set_aspect("equal")

        # Assign the number of points per row and column
        self.Nres = Nres

    def _get_spans(self, p=0.999):
        """
        Returns a lineal span for transfer angle and non-dimensional time of flight.

        Parameters
        ----------
        p: float
            Percentage of the final value. This is required due to singularities
            in some of the solvers at transfer angles of 2pi.

        Returns
        -------
        theta_span: np.array
            An array of linearly spaced transfer angles.
        tau_span: np.array
            An array of linearly spaced non-dimensional transfer times.

        """

        # Generate a meshgrid holding any combination of transfer angle and
        # non-dimensional time of flight. The 2 * pi value is avoided by
        # computing an approximate one. Nevertheless, this last value will not
        # be used due to the way `pcolor` function operates.
        theta_span, tau_span = [
            np.linspace(0, p * 2 * np.pi, self.Nres) for _ in range(2)
        ]
        return theta_span, tau_span

    def _draw_colorbar(self, maxval, step, label, cmap, color_vmin):
        """Draws the colorbar for the figure.

        Parameters
        ----------
        maxval: float
            The maximum value of the figure.
        step: float
            The step for drawing each of the colorbar ticks.
        label: str
            The title of the colorbar.
        cmap: matplotlib.cmap
            The colormap used in the contour plot.

        """

        # Generate the colorbar
        self.cbar = self.fig.colorbar(self.collection)
        self.cbar.ax.get_yaxis().set_ticks([])

        # Append the label and make its position
        self.cbar.set_label(label)
        self.cbar.ax.get_yaxis().labelpad = 15

        # Properly size the aspect ratio of the colorbar
        digits = int(np.log10(maxval)) + 1
        cbar_title = r"$\times 10^" + f"{digits-2}$" if digits > 2 else None
        self.cbar.ax.set_title(cbar_title)

        # Compute the step which separates two different levels
        step = maxval / cmap.N

        # Draw a beautiful colorbar with the legend for the number of iterations
        # in the middle
        for n in range(int(cmap.N)):

            # Select suitable font-color
            fontcolor = "black" if n != 0 else color_vmin

            # Draw the number within the scale
            self.cbar.ax.text(
                0.5 * maxval,
                step / 2 + step * n,
                str(int(step * n / 10 ** (digits - 2))),
                ha="center",
                va="center",
                color=fontcolor,
            )

    def _draw_ticks(self):
        """ Draws the ticks within the axes """

        # Set the X-ticks
        self.ax.set_xticks(np.array([0, 0.5, 1, 1.5, 2]) * np.pi)
        self.ax.set_xticklabels(
            [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{2\pi}{3}$", r"$2\pi$"]
        )

        # Set the Y-ticks
        self.ax.set_yticks(np.array([0, 0.5, 1, 1.5, 2]) * np.pi)
        self.ax.set_yticklabels(
            [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{2\pi}{3}$", r"$2\pi$"]
        )

    def _draw_labels(self):
        """ Draws axes labels """

        # Set axes labels and title
        self.ax.set_xlabel(r"$\Delta \theta$")
        self.ax.set_ylabel(r"$\Delta \tau$")


def _measure_performance(solver, theta, tau):
    """
    Computes the number of iterations from a particular value of theta and the
    transfer angle.

    Parameters
    ----------
    solver: function
        The Lambert's problem solver function.
    theta: float
        The transfer angle.
    tau: float
        The non-dimensional time of flight.

    Returns
    -------

    Notes
    -----
    The customization is null to prevent users from shooting themselves and
    creating performance comparisons under different boundary conditions.

    """

    # Generate final position vector by rotating initial one a given theta
    r1_vec = np.array([1, 0])
    r2_vec = (
        2.0 * np.array([[cos(theta), sin(theta)], [sin(theta), cos(theta)]]) @ r1_vec
    )

    # Make vectors three-dimensional
    r1_vec = np.append(r1_vec, np.array([0]), axis=0)
    r2_vec = np.append(r2_vec, np.array([0]), axis=0)

    # Compute the norms, the chord and semi-perimeter
    r1, r2 = [np.linalg.norm(rr) for rr in [r1_vec, r2_vec]]
    c = (r1 ** 2 + r2 ** 2 - 2 * r1 * r2 * cos(theta)) ** 0.5
    s = (r1 + r2 + c) / 2

    # Compute the dimensional time from the non-dimensional one using
    # Lancaster's expression. This relation is more intuitive as it relates
    # revolution cases with multiples of pi.
    mu = 1.00
    tof = tau / (8 * mu / s ** 3) ** 0.5

    # Filter non-valid input: null value is returned if no iterations were run
    if tof == 0 or theta == 0:
        return 0, 0, 0

    # Solve the problem but only collect the number of iterations
    tic = time.perf_counter()
    _, _, numiter, tpi = solver(
        mu,
        r1_vec,
        r2_vec,
        tof,
        M=0,
        prograde=True,
        maxiter=35,
        atol=1e-5,
        rtol=1e-7,
        full_output=True,
    )
    tac = time.perf_counter()

    return numiter, tpi, (tac - tic)


# Vectorize the solver
_vec_measure_performance = np.vectorize(
    _measure_performance, otypes=[np.ndarray, np.ndarray, np.ndarray], excluded=[0]
)
