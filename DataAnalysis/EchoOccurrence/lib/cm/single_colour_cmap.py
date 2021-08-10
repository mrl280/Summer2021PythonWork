import numpy as np
import matplotlib.cm as cm

from matplotlib import colors
from matplotlib.colors import ListedColormap


def single_colour_cmap(named_color):
    """

    Produce a colourmap that is just all one solid colour.
    This was developed for the purpose of shading areas boolean areas.

    :param named_color: str:
            Matplotlib named colour.
    :return new_cmap: matplotlib.colors.Colormap:
            A colourmap, all of which is just named_color
    """
    cmap_template = cm.get_cmap('jet', 256)
    new_cmap = cmap_template(np.linspace(0, 1, 256))
    new_cmap[:] = colors.to_rgba(named_color)

    return ListedColormap(new_cmap)
