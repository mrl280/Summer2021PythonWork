import numpy as np
import matplotlib.cm as cm

from matplotlib.colors import ListedColormap


def modified_viridis(levels=6):
    """
    Take the viridis colourmap, and white-out the bottom level
    :param levels: int: The number of levels in the colour bar
    :return: The viridis colourmap, but with the bottom bar whited out
    """
    jet = cm.get_cmap('viridis_r', 256)
    newcolours = jet(np.linspace(0, 1, 256))
    white = np.array([0, 0, 0, 0])  # RGBA colours

    bottom_bar = int(len(newcolours) / levels)
    newcolours[:bottom_bar, :] = white  # White out the bottom of the bar
    newcmp = ListedColormap(newcolours)

    return newcmp


if __name__ == "__main__":
    """ Testing """

    modified_viridis()
