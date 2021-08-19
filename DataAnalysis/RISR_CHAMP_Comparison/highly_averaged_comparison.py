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

from DataAnalysis.RISR_CHAMP_Comparison.read_in_flyover_events import read_in_flyover_events


def highly_averaged_comparison(in_file, lat_range, lon_width):
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
    :param lat_range: (float, float):
            The allowed latitude range.
    :param lon_width: float:
            How many degrees of longitude to use on either side of RISR radar
            Allowed longitudes will be [risr_lon - lon_width, risr_lon + lon_width]

    :return fig:  matplotlib.pyplot.figure:
            The figure. It can then be modified, added to, printed out, or saved in whatever format is desired.
    """

    altitude_range = (300, 340)  # km

    risr_lon_lat = (-94.91, 74.73)  # RISR location in geographical coordinates
    lon_range = (risr_lon_lat[0] - lon_width, risr_lon_lat[0] + lon_width)

    lat_lon_string = "Lats: " + str(round(lat_range[0], 1)) + " to " + str(round(lat_range[1], 1)) + "; " + \
                     "Lons: " + str(round(lon_range[0], 1)) + " to " + str(round(lon_range[1], 1))

    title_fontsize = 26
    density_range = (0, 25)
    bisector_colour = "red"

    event_df = read_in_flyover_events(in_file)
    print("\nHere is the event dataframe:")
    print(event_df.head())
    print("Total number of events: " + str(len(event_df)))

    print("Preparing the figure...")
    fig = plt.figure(figsize=[20, 12], constrained_layout=True, dpi=300)
    axes = add_axes(fig)
    apply_subplot_formatting(axes=axes, density_range=density_range, lat_range=lat_range, lon_width=lon_width,
                             risr_lon_lat=risr_lon_lat, bisector_colour=bisector_colour)

    fig.suptitle("RISR CHAMP Density Comparison (Highly Averaged Approach).  " + lat_lon_string +
                 "\nVarious Events from Sept-Dec 2009.  " +
                 "Produced by " + str(os.path.basename(__file__)), fontsize=title_fontsize)

    risr_scatter_color = "red"
    champ_scatter_color = "blue"

    print("Computing datetime objects and epochs...")
    pattern = '%Y-%m-%d %H:%M'
    epoch_reference_datetime = datetime.datetime(year=1970, month=1, day=1, hour=0, second=0)
    datetimes, epoch = [], []
    for i in range(len(event_df)):
        datetime_here = datetime.datetime.strptime(event_df['Time'].iat[i], pattern)
        epoch.append((datetime_here - epoch_reference_datetime).total_seconds())
        datetimes.append(datetime_here)
    event_df['datetime'] = datetimes
    event_df['epoch'] = epoch

    # We need to consider a little bit of time on either side of the flyover
    time_delta_s = 180
    event_df['starting_epoch'] = event_df['epoch'] - time_delta_s
    event_df['ending_epoch'] = event_df['epoch'] + time_delta_s

    risr_matched_medians, champ_matched_medians = [], []

    print("Looping through the events and reading in the data...")
    for i in range(len(event_df)):

        ax = axes['map'][i]
        datetime_here = event_df['datetime'].iat[i]
        start_epoch = event_df['starting_epoch'].iat[i]
        end_epoch = event_df['ending_epoch'].iat[i]

        # Grab and restrict CHAMP data
        champ_df = get_champ_df(year=datetime_here.year, month=datetime_here.month, day=datetime_here.day)
        champ_df = champ_df.loc[(champ_df['epoch'] >= start_epoch) & (champ_df['epoch'] <= end_epoch)]

        # Plot the CHAMP flyover line.  This lets us see where the satellite went event if we aren't considering any
        #  points from the flyover
        ax.plot(champ_df['gdlon'], champ_df['gdlat'], linestyle='--', linewidth=0.75, color=champ_scatter_color,
                transform=ccrs.Geodetic(), label="CHAMP Trajectory")

        champ_df['Ne'] = champ_df['Ne_cm-3'] * 1e6  # Put densities in per m^3 so they are in the same units as RISR
        champ_df = champ_df.loc[champ_df['Ne'].notna()]
        champ_df = champ_df.loc[(champ_df['gdlat'] >= lat_range[0]) & (champ_df['gdlat'] <= lat_range[1])]
        champ_df = champ_df.loc[(champ_df['gdlon'] >= lon_range[0]) & (champ_df['gdlon'] <= lon_range[1])]
        champ_df.reset_index(drop=True, inplace=True)

        if len(champ_df) > 0:
            champ_matched_medians.append(np.median(champ_df['Ne']))
        else:
            champ_matched_medians.append(np.nan)

        # Plot CHAMP measurments considered
        ax.scatter(champ_df['gdlon'], champ_df['gdlat'], marker='o', color=champ_scatter_color, s=5,
                   transform=ccrs.Geodetic(), label="CHAMP Points Considered", zorder=4)

        # Grab and restrict RISR data
        risr_df = get_risr_df(year=datetime_here.year, month=datetime_here.month,
                              file_start_day=event_df['risr_start_day'].iat[i])
        risr_df = risr_df.loc[(risr_df['epoch'] >= start_epoch) & (risr_df['epoch'] <= end_epoch)]
        risr_df = risr_df.loc[risr_df['elv'] >= 49]
        risr_df = risr_df.loc[(risr_df['gdalt'] >= altitude_range[0]) & (risr_df['gdalt'] <= altitude_range[1])]
        risr_df = risr_df.loc[(risr_df['gdlat'] >= lat_range[0]) & (risr_df['gdlat'] <= lat_range[1])]
        risr_df = risr_df.loc[(risr_df['gdlon'] >= lon_range[0]) & (risr_df['gdlon'] <= lon_range[1])]
        risr_df = risr_df.loc[risr_df['Ne'].notna()]
        risr_df.reset_index(drop=True, inplace=True)

        ax.scatter(risr_df['gdlon'], risr_df['gdlat'], marker='o', color=risr_scatter_color, s=5,
                   transform=ccrs.Geodetic(), label="RISR Points Considered")

        if len(risr_df) > 0:
            risr_matched_medians.append(np.median(risr_df['Ne']))
        else:
            risr_matched_medians.append(np.nan)

        # Title plot
        ax.set_title(str(event_df['datetime'].iat[i].date()) + "\n" + str(event_df['datetime'].iat[i].time()))
        ax.legend(loc='upper center', prop={'size': 4})

    risr_matched_medians = np.asarray(risr_matched_medians) / 1e10
    champ_matched_medians = np.asarray(champ_matched_medians) / 1e10

    # Put the matched data into a dataframe so we can view it easier
    median_matched_df = pd.DataFrame({'risr': risr_matched_medians,
                                      'champ': champ_matched_medians,
                                      'champ/risr ratio': champ_matched_medians / risr_matched_medians})

    median_matched_df = median_matched_df.loc[(median_matched_df['champ/risr ratio'].notna())]

    print("\nHere is the median matched data:")
    print(median_matched_df)

    # Complete the large scatter plot
    axes['data'].plot(median_matched_df['risr'], median_matched_df['champ'], 'o', color="purple",
                      label="Matched Medians")

    # # Do not consider any points less than 1
    # median_matched_df = median_matched_df.loc[(median_matched_df['risr'] >= 1) & (median_matched_df['champ'] >= 1)]

    # Find the best fit line, and print it out onto the plot
    x_data = median_matched_df['risr']
    y_data = median_matched_df['champ']
    A = np.vstack([x_data, np.ones(len(x_data))]).T
    m, b = np.linalg.lstsq(A, y_data, rcond=None)[0]
    axes['data'].plot(x_data, m * x_data + b, 'm',
                      label="Best fit line (y=" + str(round(m, 3)) + "x + " + str(round(b, 3)) + ")")

    axes['data'].legend(loc='upper left', prop={'size': 20})

    return fig


def get_risr_df(year, month, file_start_day):
    """

    Read in RISR Long Pulse data from a pickled file.

    :param year: int:
            The year to consider.
    :param month: int:
            The month to consider.
    :param file_start_day: int:
            The day to consider.  Often a single RISR file will contain data for several days.
            So, this is the day used in the RISR file name.  The day the event in the file starts.

    :return risr_df: pandas.DataFrame:
            Dataframe containing RISR data.  Data is unfiltered.
    """

    year = str(year)
    if month < 10:
        month = "0" + str(month)
    else:
        month = str(month)
    if file_start_day < 10:
        file_start_day = "0" + str(file_start_day)
    else:
        file_start_day = str(file_start_day)

    station = "ran"
    date = year + month + file_start_day
    loc_root = str((pathlib.Path().parent.absolute()).parent.absolute())
    in_dir = loc_root + "/DataReading/RISR/data/" + station + "/" + station + date

    try:
        in_file = in_dir + "/" + station + date + ".1min.pbz2"  # Try for 1 minute resolution data
        data_stream = bz2.BZ2File(in_file, "rb")
    except FileNotFoundError:
        in_file = in_dir + "/" + station + date + ".3min.pbz2"  # Use 3 minute data
        data_stream = bz2.BZ2File(in_file, "rb")

    print("     RISR data obtained from: + " + in_file)
    risr_df = pd.read_pickle(data_stream)

    return risr_df


def get_champ_df(year, month, day):
    """

    Read in CHAMP PLPT data from a pickled file.

    :param year: int:
            The year to consider.
    :param month: int:
            The month to consider.
    :param day: int:
            The day to consider.

    :return champ_df: pandas.DataFrame:
            Dataframe containing CHAMP data.  Data is unfiltered
    """

    in_dir = "data/champ"
    for in_file in glob.iglob(in_dir + "/*PLPT*.pbz2"):
        file_name = str(os.path.basename(in_file))
        file_year = file_name[13:17]
        file_month = file_name[18:20]
        file_day = file_name[21:23]

        # Check to see if the file is the one we are looking for
        if str(year) == file_year and str(month) in file_month and str(day) in file_day:
            print("     CHAMP data obtained from: + " + file_name)

            data_stream = bz2.BZ2File(in_file, "rb")
            return pd.read_pickle(data_stream)

    raise Exception("Error in get_champ_df(): no data found for " +
                    "year=" + str(year) + " month=" + str(month) + " day=" + str(day))


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

    lat_range = (75, 77)  # deg
    lon_width = 10  # deg
    in_file = "risr_champ_500km_conjunctions.csv"

    fig = highly_averaged_comparison(in_file=in_file, lat_range=lat_range, lon_width=lon_width)

    plt.show()

    loc_root = str(pathlib.Path().parent.absolute())
    out_dir = loc_root + "/out"
    out_fig = out_dir + "/risr_champ_highly_averaged_comparisons-all_events - lon_width" + str(lon_width)

    print("Saving plot as " + out_fig)
    fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
