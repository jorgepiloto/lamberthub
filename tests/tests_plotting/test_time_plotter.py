import matplotlib.pyplot as plt
import pytest

from lamberthub import gooding1990
from lamberthub.plotting import TimePlotter


@pytest.mark.mpl_image_compare
def test_time_per_iteration_performance_plotter():
    fig, ax = plt.subplots()
    ipp = TimePlotter(ax=ax, fig=fig)
    ipp.plot_performance(gooding1990)
    return fig
