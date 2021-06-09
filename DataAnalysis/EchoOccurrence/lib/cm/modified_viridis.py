import numpy as np
import matplotlib.cm as cm

from matplotlib.colors import ListedColormap


def modified_viridis_1():
    """
    :return: A modified viridis colormap
    """
    jet = cm.get_cmap('viridis_r', 256)
    newcolours = jet(np.linspace(0, 1, 256))
    white = np.array([0, 0, 0, 0])  # RGBA colours
    newcolours[:10, :] = white  # White out the bottom of the bar
    newcmp = ListedColormap(newcolours)
    return newcmp


def modified_viridis_2():
    """
    :return: A modified viridis colormap
    """
    jet = cm.get_cmap('viridis_r', 256)
    newcolours = jet(np.linspace(0, 1, 256))
    white = np.array([0, 0, 0, 0])  # RGBA colours
    newcolours[:50, :] = white  # White out the bottom of the bar
    newcmp = ListedColormap(newcolours)
    return newcmp