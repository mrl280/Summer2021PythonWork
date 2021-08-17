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


def highly_averaged_comparison(year, month, day, risr_start_day, start_epoch, end_epoch):
    """

    This program compares highly averaged CHAMP and RISR data.

    Instructions from Koustov:
        - Consider 3 min RISR data. Average over beams with elevation >50 deg.
        - Consider CHAMP averaged density over geographic latitudes 75 + 76 deg. Longitudinally +/- 20 deg would be
            acceptable, but +/- 10 deg is better
        - Consider those CHAMP points that are within 3 min from the center time for RISR measurement.
        - You will get ONE point for one pass over Resolute. I wonder what you would get

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
    risr_data_resolution = 3  # minutes

    risr_color = "red"
    champ_color = "blue"

    starting_datetime = datetime.datetime.utcfromtimestamp(start_epoch)  # Note we are grabbing UTC datetime objs
    ending_datetime = datetime.datetime.utcfromtimestamp(end_epoch)

    proj = ccrs.NorthPolarStereo()  # RISR is in the northern hemisphere

    print("Preparing the figure...")
    fig = plt.figure(figsize=[8, 5], constrained_layout=True, dpi=300)
    gs = fig.add_gridspec(ncols=2, nrows=1)
    map_axis = fig.add_subplot(gs[0, 0], projection=proj)
    data_axis = fig.add_subplot(gs[0, 1])

    fig.suptitle("RISR CHAMP Density Comparison (Highly Averaged Approach)"
                 "\nScatter from " + str(starting_datetime) + " to " + str(ending_datetime) + " (UTC)" +
                 "\nProduced by " + str(os.path.basename(__file__)), fontsize=title_fontsize)

    print("Formatting the map...")
    map_axis.set_extent([-180, 180, 90, lat_extreme], crs=ccrs.PlateCarree())
    map_axis.gridlines()
    # map_axis.add_feature(cfeature.OCEAN)
    map_axis.add_feature(cfeature.LAND)

    # Plot RISR as a triangle
    map_axis.plot([risr_lon, risr_lon], [risr_lat, risr_lat], marker="^", color=risr_color, markersize=5,
                  transform=ccrs.Geodetic(), label="RISR")

    print("\nAdding RISR scatter to the map plot...")
    station = "ran"
    date = str(year) + str(month) + str(risr_start_day)
    loc_root = str((pathlib.Path().parent.absolute()).parent.absolute())
    in_dir = loc_root + "/DataReading/RISR/data/" + station + "/" + station + date
    in_file = in_dir + "/" + station + date + "." + str(risr_data_resolution) + "min.pbz2"
    print("     RISR data obtained from: + " + in_file)

    data_stream = bz2.BZ2File(in_file, "rb")
    risr_df = pd.read_pickle(data_stream)

    # Restrict RISR data
    risr_df = risr_df.loc[(risr_df['epoch'] >= start_epoch) & (risr_df['epoch'] <= end_epoch)]
    risr_df = risr_df.loc[risr_df['elv'] >= 49]
    risr_df = risr_df.loc[(risr_df['gdalt'] >= 300) & (risr_df['gdalt'] <= 340)]
    risr_df.reset_index(drop=True, inplace=True)

    map_axis.scatter(risr_df['gdlon'], risr_df['gdlat'], marker='o', color=risr_color, s=5, transform=ccrs.Geodetic(),
                     label="RISR Measurements")

    print("\nAdding champ scatter to the plot...")
    champ_data_found = False
    in_dir = "data/champ"
    for in_file in glob.iglob(in_dir + "/*PLPT*.pbz2"):
        if champ_data_found:
            continue  # Nothing more to do here

        file_name = str(os.path.basename(in_file))

        # Check to see if the file is the one we are looking for
        if str(year) in file_name and str(month) in file_name and str(day) in file_name:
            champ_data_found = True
            print("     CHAMP data obtained from: + " + file_name)

            data_stream = bz2.BZ2File(in_file, "rb")
            champ_df = pd.read_pickle(data_stream)

    # Restrict CHAMP data
    champ_df = champ_df.loc[(champ_df['epoch'] >= start_epoch) & (champ_df['epoch'] <= end_epoch)]
    champ_df.reset_index(drop=True, inplace=True)

    map_axis.scatter(champ_df['gdlon'], champ_df['gdlat'], marker='o', color=champ_color, s=5,
                     transform=ccrs.Geodetic(),
                     label="CHAMP Scatter")

    print("     Here is CHAMP's max altitude: " + str(np.max(champ_df['radius_km']) - 6358))
    print("     Here is CHAMP's min altitude: " + str(np.min(champ_df['radius_km']) - 6358))

    # TODO: Draw valid area of overlap on the map and compute the ONE point for this overlap

    return fig


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

    fig = highly_averaged_comparison(year=year, month=month, day=day, risr_start_day=risr_start_day,
                                     start_epoch=start_epoch, end_epoch=end_epoch)

    plt.show()

    out_fig = out_dir + "/risr_champ_highly_averaged_comparison-" + str(year) + str(month) + str(day)

    # print("Saving plot as " + out_fig)
    # fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
