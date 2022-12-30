import matplotlib.pyplot as plt

from lamberthub.plotting import IterationsPlotter
from lamberthub.universal_solvers.gooding import gooding1990



plotter = IterationsPlotter()
plotter.plot_performance(gooding1990)
plt.show()
