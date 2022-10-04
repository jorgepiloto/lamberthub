"""A module containing various plotting utilities."""

import matplotlib.pyplot as plt
import matplotlib as mpl

# CHECK: https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar

def generate_discrete_cmap(base_cmap_name=None, n_levels=10, desired_cmap_name=None):
    """Generate a discrete contour map with desired number of levels.

    Parameters
    ----------
    base_cmap_name : str, optional
        Name of the ``LinearSegmentedColormap`` to be used as base map.
    n_levels : int
        Desired number of levels for the map.

    Returns
    -------
    cmap : ~matplotlib.LinearSegmentedColormap
        Discrete contour map.

    """
    # If no base name was provided, use the yellow + orange + red
    if base_cmap_name is None:
        base_cmap_name = "YlOrRd"

    # If no output name is provided, name it "lamberthub"
    if desired_cmap_name is None:
        desired_cmap_name = "lamberthub_cmap"

    # Get the desired colormap 
    base_cmap = mpl.colormaps[base_cmap_name]

    # Retireve all its colors
    cmaplist = [base_cmap(i) for i in range(base_cmap.N)]

    # Force the first color to be black
    cmaplist[0] = (0, 0, 0, 1.0)

    return mpl.colors.LinearSegmentedColormap.from_list(desired_cmap_name, cmaplist, cmap.N)
