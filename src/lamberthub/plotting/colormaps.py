"""A module containing various colormaps definitions."""

from matplotlib.colors import LinearSegmentedColormap


def generate_cmap_from_colors_list(colors_list, cmap_name, n_levels):
    return LinearSegmentedColormap.from_list(
        cmap_name, colors_list, N=n_levels
    )


sunshine_9lev = generate_cmap_from_colors_list(
    [
        (255, 255, 255),
        (255, 245, 204),
        (255, 230, 112),
        (255, 204, 51),
        (255, 175, 51),
        (255, 153, 51),
        (255, 111, 51),
        (255, 85, 0),
        (230, 40, 30),
        (200, 30, 20),
    ],
    "sunshine_9lev",
    10,
)
