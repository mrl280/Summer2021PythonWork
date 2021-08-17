import bz2
import datetime
import glob
import os
import pathlib

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt

from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch


def risr_champ_conjunction_map(year, month, day, risr_start_day, start_epoch, end_epoch):
    """

    Scatter all data for a RISR/CHAMP conjunction on a map.  Hopefully this will give us a better idea of what is going
     on.

    Program also prints out some helpful information.

    :param year: int:
            The year of the event to plot.
    :param month: int:
            The month of the event to plot.
    :param day: int
            The day of the event to plot.
    :param risr_start_day: int:
            Often a single RISR file will contain data for several days.
            So, this is the day used in the RISR file name.
    :param start_epoch: int:
            Event starting epoch.
    :param end_epoch: int:
            Event ending epoch.

    :return fig:  matplotlib.pyplot.figure:
            The figure. It can then be modified, added to, printed out, or saved in whatever format is desired.
    """

    lat_extreme = 70  # deg
    risr_lon, risr_lat = -94.91, 74.73  # RISR
    title_fontsize = 14

    risr_color = "red"
    champ_color = "blue"

    starting_datetime = datetime.datetime.utcfromtimestamp(start_epoch)  # Note we are grabbing UTC datetime objs
    ending_datetime = datetime.datetime.utcfromtimestamp(end_epoch)

    # RISR is in the northern hemisphere
    proj = ccrs.NorthPolarStereo()

    fig = plt.figure(figsize=[7, 8], constrained_layout=True, dpi=300)
    gs = fig.add_gridspec(ncols=1, nrows=1)
    ax = fig.add_subplot(gs[0, 0], projection=proj)

    fig.suptitle("RISR CHAMP Conjunction. "
                 "\nScatter from " + str(starting_datetime) + " to " + str(ending_datetime) + " (UTC)." +
                 "\nProduced by " + str(os.path.basename(__file__)), fontsize=title_fontsize)

    ax.set_extent([-180, 180, 90, lat_extreme], crs=ccrs.PlateCarree())
    ax.gridlines()
    # ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)

    # Plot RISR as a red dot
    ax.plot([risr_lon, risr_lon], [risr_lat, risr_lat], marker="^", color=risr_color, markersize=5,
            transform=ccrs.Geodetic(), label="RISR")

    print("\nAdding RISR scatter to the plot...")
    add_risr_points(ax=ax, year=year, month=month, day=risr_start_day, start_epoch=start_epoch, end_epoch=end_epoch,
                    color=risr_color)

    print("\nAdding champ scatter to the plot...")
    add_champ_points(ax=ax, year=year, month=month, day=day, start_epoch=start_epoch, end_epoch=end_epoch,
                     color=champ_color)

    plt.legend(loc='upper right', prop={'size': 20})

    return fig


def add_risr_points(ax, year, month, day, start_epoch, end_epoch, color):
    """

    Add a scatter point for each RISR measurement - this should show us our different beam options.

    :param ax: matplotlib.axes:
            The axes to draw on.
    :param year: int:
            The year of the event to plot.
    :param month: int:
            The month of the event to plot.
    :param day: int
            The day of the event to plot.
    :param start_epoch: int:
            Event starting epoch.
    :param end_epoch: int:
            Event ending epoch.
    :param color: str:
            matplotlib named colour to use for the scatter.
    """

    station = "ran"
    date = str(year) + str(month) + str(day)
    minute_res = 1

    loc_root = str((pathlib.Path().parent.absolute()).parent.absolute())
    in_dir = loc_root + "/DataReading/RISR/data/" + station + "/" + station + str(year) + str(month) + str(day)

    in_file = in_dir + "/" + station + date + "." + str(minute_res) + "min.pbz2"
    print(in_file)

    data_stream = bz2.BZ2File(in_file, "rb")
    df = pd.read_pickle(data_stream)

    df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]

    print("Here is RISR's max altitude: " + str(np.max(df['gdalt'])))
    print("Here is RISR's min altitude: " + str(np.min(df['gdalt'])))

    print("Here are RISR's elevations:")
    print(df['elv'].unique())

    df = df.loc[df['elv'] >= 49]

    print("Restricting to RISR altitude...")
    df = df.loc[(df['gdalt'] >= 300) & (df['gdalt'] <= 340)]

    print(df.head())

    print("Here are the RISR beams with measurements at this altitude: ")
    print(df['wdBmnum'].unique())

    ax.scatter(df['gdlon'], df['gdlat'], marker='o', color=color, s=5, transform=ccrs.Geodetic(),
               label="RISR Measurements")


def add_champ_points(ax, year, month, day, start_epoch, end_epoch, color):
    """

    Add a scatter point for each CHAMP measurement - this should tell us where about it flew

    :param ax: matplotlib.axes:
            The axes to draw on.
    :param year: int:
            The year of the event to plot.
    :param month: int:
            The month of the event to plot.
    :param day: int
            The day of the event to plot.
    :param start_epoch: int:
            Event starting epoch.
    :param end_epoch: int:
            Event ending epoch.
    :param color: str:
            matplotlib named colour to use for the scatter.
    """

    in_dir = "data/champ"
    for in_file in glob.iglob(in_dir + "/*PLPT*.pbz2"):
        file_name = str(os.path.basename(in_file))

        # Check to see if the file is the one we are looking for
        if str(year) in file_name and str(month) in file_name and str(day) in file_name:
            print("We are going to use the following .dat file: " + file_name)

            data_stream = bz2.BZ2File(in_file, "rb")
            df = pd.read_pickle(data_stream)

            df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]

            print(df.head())

            ax.scatter(df['gdlon'], df['gdlat'], marker='o', color=color, s=5, transform=ccrs.Geodetic(),
                       label="CHAMP Measurements")

            print("Here is CHAMP's max altitude: " + str(np.max(df['radius_km']) - 6358))
            print("Here is CHAMP's min altitude: " + str(np.min(df['radius_km']) - 6358))


if __name__ == '__main__':
    """ Testing """

    year = 2009
    month = 10
    day = 16

    # Sometimes RISR events cover multiple days - this the date on the RISR file
    risr_start_day = 15

    # On Oct 16, 2009 - CHAMP is over RISR at 22:51:56
    # _, start_epoch = build_datetime_epoch(year=year, month=month, day=day, hour=22, minute=47, second=00)
    # _, end_epoch = build_datetime_epoch(year=year, month=month, day=day, hour=22, minute=57, second=00)
    _, start_epoch = build_datetime_epoch(year=year, month=month, day=day, hour=22, minute=50, second=00)
    _, end_epoch = build_datetime_epoch(year=year, month=month, day=day, hour=22, minute=54, second=00)

    loc_root = str(pathlib.Path().parent.absolute())
    out_dir = loc_root + "/out"

    fig = risr_champ_conjunction_map(year=year, month=month, day=day, risr_start_day=risr_start_day,
                                     start_epoch=start_epoch, end_epoch=end_epoch)

    plt.show()

    out_fig = out_dir + "/risr_champ_conj-" + str(year) + str(month) + str(day)

    print("Saving plot as " + out_fig)
    fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
