""" Holds plotting utilities related with required number of iterations """

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from numpy import cos, pi, sin


def _iterations_from_theta_and_T(solver, theta, T):
    """
    Computes the number of iterations from a particular value of theta and the
    transfer angle.

    Parameters
    ----------
    solver: function
        The Lambert's problem solver function.
    theta: float
        The transfer angle.
    T: float
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
    tof = T / (8 * mu / s ** 3) ** 0.5

    # Filter non-valid input: null value is returned if no iterations were run
    if tof == 0 or theta == 0:
        return 0

    # Solve the problem but only collect the number of iterations
    _, _, numiter = solver(
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

    return numiter


# Vectorize the solver
_iterations_from_theta_and_T_vec = np.vectorize(
    _iterations_from_theta_and_T, otypes=[np.ndarray], excluded=[0]
)


class IterationsPlotter:
    """
    A class for plotting solvers number of iterations as functions of the
    transfer angle and the non-dimensional time of flight.
    """

    def __init__(self, ax, fig):
        """
        Initializes any instance of the plotter.

        Parameters
        ----------
        ax: matplotlib.Axes
            The axes in which the lines will be drawn.
        fig: matplotlib.Figure
            The figure instance for the plot.

        """

        # Link Axes and Figures to the instance
        self.fig, self.ax = fig, ax

        # Force the aspect ratio of the figure
        self.ax.set_aspect("equal")

    def _get_solver_name(self, solver):
        """
        Returns solver's name.

        Parameters
        ----------
        solver: function
            Solver function.

        Returns
        -------
        name: str
            String representation for the solver.

        """

        # Get the function name and its length
        raw_name = str(solver.__name__.capitalize())
        len_name = len(raw_name)

        # Break into author name and implementation year
        author_name = raw_name[0 : (len_name - 4)]
        year_name = raw_name[(len_name - 4) :]

        # Assemble and return the solver's name
        name = author_name + " " + year_name
        return name

    def plot_performance(self, solver, maxiter=10, cmap=None):
        """
        Returns a graphical representation on the iteration performance for a
        particular solver.

        Parameters
        ----------
        solver: function
            The solver who's performance is to be assessed.

        Notes
        -----
        The method does not accept neither a transfer angle nor time of flight
        span vectors. This is to force a common basis when addressing the
        performance from the number of iterations point of view. Only figure
        customization options are valid as input parameters.

        """

        # Generate a meshgrid holding any combination of transfer angle and
        # non-dimensional time of flight. The 2 * pi value is avoided by
        # computing an approximate one. Nevertheless, this last value will not
        # be used due to the way `pcolor` function operates.
        theta_span, T_span = [np.linspace(0, 2 * pi * 0.999, 50) for _ in range(2)]
        THETA, TOF = np.meshgrid(theta_span, T_span)

        # Compute the meshgrid holding the number of required iterations
        NN_ITER = _iterations_from_theta_and_T_vec(
            solver, theta_span[np.newaxis, :], T_span[:, np.newaxis]
        ).astype("int")

        # Prepare the levels for the contour
        levels = MaxNLocator(nbins=maxiter + 1).tick_values(0, maxiter)
        cmap = plt.get_cmap("YlOrRd") if cmap is None else cmap
        bd_norm = BoundaryNorm(levels, ncolors=cmap.N)

        # Generate meshgrid
        c = self.ax.pcolor(
            THETA,
            TOF,
            NN_ITER[:-1, :-1],  # For pcolor, the last rows need to be removed
            cmap=cmap,
            vmin=1,
            edgecolors="k",
            linewidths=1,
            norm=bd_norm,
        )
        c.cmap.set_under("k")
        cbar = self.fig.colorbar(c)
        cbar.ax.get_yaxis().set_ticks([])

        # Draw a beautiful colorbar with the legend for the number of iterations
        # in the middle
        for n in range(maxiter):
            color = "white" if n == 0 else "black"
            cbar.ax.text(
                4.5,
                (maxiter * n + 5) / maxiter,
                str(n),
                ha="center",
                va="center",
                color=color,
            )
            cbar.ax.get_yaxis().labelpad = maxiter * 1.5
        cbar.set_label("Number of iterations")

        # Set the ticks
        self.ax.set_xticks(np.array([0, 0.5, 1, 1.5, 2]) * pi)
        self.ax.set_xticklabels(
            [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{2\pi}{3}$", r"$2\pi$"]
        )
        self.ax.set_yticks(np.array([0, 0.5, 1, 1.5, 2]) * pi)
        self.ax.set_yticklabels(
            [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{2\pi}{3}$", r"$2\pi$"]
        )

        # Set axes labels and title
        self.ax.set_xlabel(r"$\Delta \theta$")
        self.ax.set_ylabel(r"$\Delta \tau$")
        self.ax.set_title(self._get_solver_name(solver))

        return self.ax
