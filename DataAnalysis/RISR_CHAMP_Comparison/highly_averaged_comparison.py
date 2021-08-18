import bz2
import datetime
import glob
import os
import pathlib

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
import matplotlib.patches as patches

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch
from DataAnalysis.RISR_CHAMP_Comparison.read_in_flyover_events import read_in_flyover_events


def highly_averaged_comparison(in_file):
    """

    This program compares highly averaged CHAMP and RISR data.

    This programs runs on all the events listed in in_file.
     It you want to investigate a single event, see highly_averaged_comparison_single_event()

    Instructions from Koustov:
        - Consider 3 min RISR data. Average over beams with elevation >50 deg.
        - Consider CHAMP averaged density over geographic latitudes 75 + 76 deg. Longitudinally +/- 20 deg would be
            acceptable, but +/- 10 deg is better.  Updated instruction: It looks as you have only 1 or 2 beams to
            average for 1 degree of RISR coverage in latitude. This does not look good, perhaps we have to go back to
            2 degrees of coverage.
        - Consider those CHAMP points that are within 3 min from the center time for RISR measurement.
        - You will get ONE point for one pass over Resolute. I wonder what you would get

    :param in_file: str:
            The name of the .csv file containing the list of flyovers.
            Should contain the following columns:

                Time:               datetime of flyover
                Distance:           great circle distance between CHAMP and RISR at time of flyover
                champ_data:         is there CHAMP data for this event?  YES or NO
                risr_data:          is there RISR data for this event? YES or NO
                risr_start_day:     Often a single RISR file will contain data for several days. So, this is the day
                                     used in the RISR file name.

    :return fig:  matplotlib.pyplot.figure:
            The figure. It can then be modified, added to, printed out, or saved in whatever format is desired.
    """

    lat_range = (75, 77)  # deg
    lon_width = 10  # deg

    risr_lon_lat = (-94.91, 74.73)  # RISR location in geographical coordinates

    title_fontsize = 30
    risr_data_resolution = 3  # minutes
    density_range = (0, 30)
    bisector_colour = "red"

    event_df = read_in_flyover_events(in_file)
    print("Here is the event dataframe:")
    print(event_df.head())

    print("Preparing the figure...")
    fig = plt.figure(figsize=[20, 12], constrained_layout=True, dpi=300)
    axes = add_axes(fig)
    apply_subplot_formatting(axes=axes, density_range=density_range, lat_range=lat_range, lon_width=lon_width,
                             risr_lon_lat=risr_lon_lat, bisector_colour=bisector_colour)

    fig.suptitle("RISR CHAMP Density Comparison (Highly Averaged Approach)"
                 "\nVarious Events from Sept-Dec 2009.  " +
                 "Produced by " + str(os.path.basename(__file__)), fontsize=title_fontsize)

    risr_scatter_color = "red"
    champ_scatter_color = "blue"

    # starting_datetime = datetime.datetime.utcfromtimestamp(start_epoch)  # Note we are grabbing UTC datetime objs
    # ending_datetime = datetime.datetime.utcfromtimestamp(end_epoch)

    # TODO: Add data to the plots

    return fig


def apply_subplot_formatting(axes, density_range, lat_range, lon_width, risr_lon_lat, bisector_colour):
    """"

    Format subplots.

    This includes plotting the RISR radar location itself and shading a valid areas of spatial overlap.

    :param axes: dictionary of matplotlib.axes:
            A dictionary containing all the axes items.
    :param density_range: (float, float):
            The x and y limits for the density comparison scatter plots.
    :param lat_range:
    :param lon_width: float:
            The width of longitude around RISR we accept
    :param risr_lon_lat: (float, float):
            The location of RISR radar in geographical coordinates (longitude, latitude)
    :param bisector_colour: str:
            The colour of the bisector, any matplotlib named colour.
    """

    risr_lon = risr_lon_lat[0]
    risr_lat = risr_lon_lat[1]

    for axis_key, axis_item in axes.items():

        if axis_key == "map":
            print("Formatting map subplots...")

            for event_index, ax in axis_item.items():
                # ax.set_extent([-180, 180, 90, lat_lower_extent], crs=ccrs.PlateCarree())
                ax.set_extent([-135, -45, 80, 70], crs=ccrs.PlateCarree())
                ax.gridlines()
                # ax.add_feature(cfeature.OCEAN)
                ax.add_feature(cfeature.LAND)

                # Create a rectangular patch on the map to outline the valid area of overlap
                # It looks better if we use two smaller patches - gives a bit of a curve.
                # We have to recreate these rectangles every time because you can't reuse an artist on multiple axes
                rect1 = patches.Rectangle(xy=(risr_lon - lon_width, lat_range[0]), width=lon_width,
                                          height=(lat_range[1] - lat_range[0]), fill=True, alpha=0.2,
                                          transform=ccrs.Geodetic(),
                                          label="Spatial Overlap Considered")
                rect2 = patches.Rectangle(xy=(risr_lon, lat_range[0]), width=lon_width,
                                          height=(lat_range[1] - lat_range[0]),
                                          fill=True, alpha=0.2, transform=ccrs.Geodetic())
                ax.add_patch(rect1)
                ax.add_patch(rect2)

                # Plot RISR as a triangle
                ax.plot(risr_lon, risr_lat, "*", color="green", markersize=5, transform=ccrs.Geodetic(),
                        label="RISR Radar")

        elif axis_key == "data":
            print("Formatting the large scatter plot...")
            ax = axis_item
            ax.set_ylim(density_range)
            ax.set_xlim(density_range)

            ax.set_xlabel("Density (RISR) $10^{10}$ [m$^{-3}$]", fontsize=24)
            ax.set_ylabel("Density (CHAMP) $10^{10}$ [m$^{-3}$]", fontsize=24)

            ax.yaxis.set_major_locator(MultipleLocator(5))
            ax.xaxis.set_major_locator(MultipleLocator(5))
            ax.yaxis.set_minor_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(1))

            ax.set_title("CHAMP/RISR Density Comparison", fontsize=24)

            ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
            ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
            ax.tick_params(axis='both', which='major', labelsize=20)
            ax.plot([ax.get_ylim()[0], ax.get_ylim()[1]], [ax.get_xlim()[0], ax.get_xlim()[1]],
                    linestyle='--', linewidth=2, color=bisector_colour)

        else:
            raise Exception("Error in apply_subplot_formatting(): axis_key of " + str(axis_key) + " not recognized.")


def add_axes(fig):
    """
    :param fig: matplotlib.pyplot.figure:
            The figure to draw the axes on

    :return axes: dictionary of dictionary of axes:
            A dictionary containing all the axes items.
    """

    gs = fig.add_gridspec(ncols=15, nrows=3)

    axes = dict()  # Remember that splices don't include last index
    proj = ccrs.NorthPolarStereo()  # RISR is in the Northern Hemisphere

    # We will have a small map for each flyover event
    axes['map'] = {0: fig.add_subplot(gs[0, 0], projection=proj),
                   1: fig.add_subplot(gs[1, 0], projection=proj),
                   2: fig.add_subplot(gs[2, 0], projection=proj),

                   3: fig.add_subplot(gs[0, 1], projection=proj),
                   4: fig.add_subplot(gs[1, 1], projection=proj),
                   5: fig.add_subplot(gs[2, 1], projection=proj),

                   6: fig.add_subplot(gs[0, 2], projection=proj),
                   7: fig.add_subplot(gs[1, 2], projection=proj),
                   8: fig.add_subplot(gs[2, 2], projection=proj),

                   9: fig.add_subplot(gs[0, 3], projection=proj),
                   10: fig.add_subplot(gs[1, 3], projection=proj),
                   11: fig.add_subplot(gs[2, 3], projection=proj),

                   12: fig.add_subplot(gs[0, 4], projection=proj),
                   13: fig.add_subplot(gs[1, 4], projection=proj),
                   14: fig.add_subplot(gs[2, 4], projection=proj),

                   15: fig.add_subplot(gs[0, 5], projection=proj),
                   16: fig.add_subplot(gs[1, 5], projection=proj),
                   17: fig.add_subplot(gs[2, 5], projection=proj),

                   18: fig.add_subplot(gs[0, 6], projection=proj),
                   19: fig.add_subplot(gs[1, 6], projection=proj),
                   20: fig.add_subplot(gs[2, 6], projection=proj),

                   21: fig.add_subplot(gs[0, 7], projection=proj),
                   22: fig.add_subplot(gs[1, 7], projection=proj),
                   23: fig.add_subplot(gs[2, 7], projection=proj),
                   }

    # We will have one large plot for the scatter
    axes['data'] = fig.add_subplot(gs[0:3, 8:16])

    return axes


if __name__ == '__main__':
    """ Testing """

    loc_root = str(pathlib.Path().parent.absolute())
    out_dir = loc_root + "/out"

    in_file = "risr_champ_500km_conjunctions.csv"

    fig = highly_averaged_comparison(in_file=in_file)

    plt.show()

    out_fig = out_dir + "/risr_champ_highly_averaged_comparisons-all_events"

    print("Saving plot as " + out_fig)
    # fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
