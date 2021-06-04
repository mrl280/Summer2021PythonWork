import numpy as np
import matplotlib.cm as cm

from matplotlib.colors import ListedColormap


def modified_jet():
    """
    :return: A modified jet colormap
    """
    jet = cm.get_cmap('jet', 256)
    newcolours = jet(np.linspace(0, 1, 256))
    purple = np.array([155/256, 53/256, 161/256, 1])  # RGBA colours
    white = np.array([0, 0, 0, 0])
    newcolours[:35, :] = white  # Make the first few colours white
    newcolours[220:255, :] = purple  # Make the last few colours purple
    newcmp = ListedColormap(newcolours)
    return newcmp