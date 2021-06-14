import numpy as np
import matplotlib.cm as cm

from matplotlib.colors import ListedColormap


def modified_jet(levels=6):
    """
    Take the jet colourmap, and white-out the bottom level
    :param levels: int: The number of levels in the colour bar
    :return: The jet colourmap, but with the bottom bar whited out
    """
    jet = cm.get_cmap('jet', 256)
    newcolours = jet(np.linspace(0, 1, 256))
    white = np.array([0, 0, 0, 0])  # RGBA colours

    bottom_bar = int(len(newcolours) / levels)
    newcolours[:bottom_bar, :] = white  # White out the bottom of the bar
    newcmp = ListedColormap(newcolours)

    return newcmp
