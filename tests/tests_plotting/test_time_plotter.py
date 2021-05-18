import matplotlib.pyplot as plt
import pytest

from lamberthub import gooding1990
from lamberthub.plotting import TPIPlotter, TTCPlotter


@pytest.mark.mpl_image_compare
def test_time_per_iteration_performance_plotter():
    fig, ax = plt.subplots()
    ipp = TPIPlotter(ax=ax, fig=fig)
    ipp.plot_performance(gooding1990, N_samples=3)
    return fig


@pytest.mark.mpl_image_compare
def test_total_time_performance_plotter():
    fig, ax = plt.subplots()
    ipp = TTCPlotter(ax=ax, fig=fig)
    ipp.plot_performance(gooding1990, N_samples=3)
    return fig
