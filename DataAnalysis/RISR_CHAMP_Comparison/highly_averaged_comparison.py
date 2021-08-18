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
    density_range = (0, 15)
    bisector_colour = "red"

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
    map_axis.plot(risr_lon, risr_lat, "^", color=risr_color, markersize=5, transform=ccrs.Geodetic(),
                  label="RISR Radar")

    # Create a Rectangle patch outlining the valid area of overlap
    lat_min = 75
    lat_max = 77
    lon_width = 10
    rect1 = patches.Rectangle(xy=(risr_lon - lon_width, lat_min), width=lon_width, height=(lat_max - lat_min),
                              fill=True, alpha=0.2, transform=ccrs.Geodetic(), label="Spatial Overlap Considered")
    rect2 = patches.Rectangle(xy=(risr_lon, lat_min), width=lon_width, height=(lat_max - lat_min), fill=True, alpha=0.2,
                              transform=ccrs.Geodetic())
    map_axis.add_patch(rect1)
    map_axis.add_patch(rect2)
    map_axis.set_title("Lats: " + str(lat_min) + " to " + str(lat_max) + ";" +
                       " Lons: " + str(risr_lon - lon_width) + " to " + str(risr_lon + lon_width))

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
    risr_df = risr_df.loc[(risr_df['gdlat'] >= lat_min) & (risr_df['gdlat'] <= lat_max)]
    risr_df = risr_df.loc[(risr_df['gdlon'] >= risr_lon - lon_width) & (risr_df['gdlon'] <= risr_lon + lon_width)]
    risr_df = risr_df.loc[risr_df['Ne'].notna()]
    risr_df.reset_index(drop=True, inplace=True)

    map_axis.scatter(risr_df['gdlon'], risr_df['gdlat'], marker='o', color=risr_color, s=5, transform=ccrs.Geodetic(),
                     label="RISR Measurements Considered")

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
    champ_df = champ_df.loc[(champ_df['gdlat'] >= lat_min) & (champ_df['gdlat'] <= lat_max)]
    champ_df = champ_df.loc[(champ_df['gdlon'] >= risr_lon - lon_width) & (champ_df['gdlon'] <= risr_lon + lon_width)]

    champ_df.reset_index(drop=True, inplace=True)

    map_axis.scatter(champ_df['gdlon'], champ_df['gdlat'], marker='o', color=champ_color, s=5,
                     transform=ccrs.Geodetic(),
                     label="CHAMP Measurements Considered")

    print("     Here is CHAMP's max altitude: " + str(np.max(champ_df['radius_km']) - 6358))
    print("     Here is CHAMP's min altitude: " + str(np.min(champ_df['radius_km']) - 6358))

    print("Formatting the data axis...")
    data_axis.set_ylim(density_range)
    data_axis.set_xlim(density_range)

    data_axis.set_xlabel("Density (RISR) $10^{10}$ [m$^{-3}$]")
    data_axis.set_ylabel("Density (CHAMP) $10^{10}$ [m$^{-3}$]")

    data_axis.yaxis.set_major_locator(MultipleLocator(5))
    data_axis.xaxis.set_major_locator(MultipleLocator(5))
    data_axis.yaxis.set_minor_locator(MultipleLocator(1))
    data_axis.xaxis.set_minor_locator(MultipleLocator(1))

    data_axis.set_title("CHAMP/RISR Density Comparison")

    data_axis.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
    data_axis.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
    data_axis.plot([data_axis.get_ylim()[0], data_axis.get_ylim()[1]],
                   [data_axis.get_xlim()[0], data_axis.get_xlim()[1]],
                   linestyle='--', linewidth=2, color=bisector_colour)

    # Compute density averages for each instrument
    print("\nHere are the RISR densities in m^-3: ")
    print(risr_df['Ne'])
    mean_density_risr = np.mean(risr_df['Ne']) / 1e10
    median_density_risr = np.median(risr_df['Ne']) / 1e10
    print("The mean RISR density is: " + str(mean_density_risr) + " m^3")
    print("The median RISR density is: " + str(median_density_risr) + " m^3")

    champ_df['Ne'] = champ_df['Ne_cm-3'] * 1e6  # Put all densities in per m^3
    print("\nHere are the CHAMP densities in m^-3: ")
    print(champ_df['Ne'])
    mean_density_champ = np.mean(champ_df['Ne']) / 1e10
    median_density_champ = np.median(champ_df['Ne']) / 1e10
    print("The mean CHAMP density is: " + str(mean_density_champ) + " 10^10 m^3")
    print("The median CHAMP density is: " + str(median_density_champ) + " 10^10 m^3")

    data_axis.plot(mean_density_risr, mean_density_champ, 'o', color="green", label="mean")
    data_axis.plot(median_density_risr, median_density_champ, 'o', color="purple", label="median")

    map_axis.legend(loc='upper left')
    data_axis.legend(loc='upper left')

    return fig


if __name__ == '__main__':
    """ Testing """

    # On Oct 16, 2009 - CHAMP is over RISR at 22:51:56 - first event Koustov suggested
    year = 2009
    month = 10
    day = 16
    risr_start_day = 15  # Sometimes RISR events cover multiple days - this the date on the RISR file
    _, start_epoch = build_datetime_epoch(year=year, month=month, day=day, hour=22, minute=49, second=00)
    _, end_epoch = build_datetime_epoch(year=year, month=month, day=day, hour=22, minute=55, second=00)

    # On Oct 16, 2009 - CHAMP is over RISR at  09:22:15 - low density
    # year = 2009
    # month = 10
    # day = 16
    # risr_start_day = 15  # Sometimes RISR events cover multiple days - this the date on the RISR file
    # _, start_epoch = build_datetime_epoch(year=year, month=month, day=day, hour=9, minute=20, second=00)
    # _, end_epoch = build_datetime_epoch(year=year, month=month, day=day, hour=9, minute=24, second=00)

    # # On Oct 15, 2009 - CHAMP is over RISR at  22:39:01 - probably too far away
    # year = 2009
    # month = 10
    # day = 15
    # risr_start_day = 15  # Sometimes RISR events cover multiple days - this the date on the RISR file
    # _, start_epoch = build_datetime_epoch(year=year, month=month, day=day, hour=22, minute=37, second=00)
    # _, end_epoch = build_datetime_epoch(year=year, month=month, day=day, hour=22, minute=41, second=00)

    loc_root = str(pathlib.Path().parent.absolute())
    out_dir = loc_root + "/out"

    fig = highly_averaged_comparison(year=year, month=month, day=day, risr_start_day=risr_start_day,
                                     start_epoch=start_epoch, end_epoch=end_epoch)

    plt.show()

    out_fig = out_dir + "/risr_champ_highly_averaged_comparison-" + \
              str(year) + str(month) + str(day) + "-" + str(start_epoch)

    print("Saving plot as " + out_fig)
    fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
