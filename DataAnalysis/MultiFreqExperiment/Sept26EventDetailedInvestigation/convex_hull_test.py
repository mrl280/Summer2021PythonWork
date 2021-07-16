import calendar
import os
import pathlib
import warnings

import numpy as np
import pandas as pd
import pydarn

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, NullFormatter
from scipy.spatial import ConvexHull, convex_hull_plot_2d

from DataAnalysis.EchoOccurrence.lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch
from lib.basic_SD_df_filter import basic_SD_df_filter


def convex_hull_test(single_day_df, beam_range, gate_range, hour_range=None, area=None):
    """

    To see how the gate-time profiles compare at the different frequencies, try wrapping each in a convex hull and
    plotting them ontop of each other

    :param single_day_df: pandas.Data_Frame:
            A single days worth of data that we want to range profile
    :param gate_range: (<int>, <int>):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param hour_range: (<int>, <int>) (optional, default is None):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not quite inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
            If None, then hour range is determined automatically from the data
    :param area: int:
            The numbered area of interest.

    :return: matplotlib.pyplot.figure:
            The figure, it can then be viewed, modified, or saved to file
    """

    hull_colours = {10: 'red',
                    12: 'blue',
                    13: 'green',
                    14: 'magenta'}

    if len(single_day_df) <= 0:
        warnings.warn("There is no data in the provided dataframe", category=Warning)

    df = single_day_df.copy()  # Convenience
    time_units = 'ut'
    frequencies = [10, 12, 13, 14]

    print("Filtering data..")
    # Get some information about the station we are working with
    station = df['station'].iat[0]
    radar_id = pydarn.read_hdw_file(station).stid

    # Make sure all of the data is of the same spatial resolution
    spatial_resolution = df['rangeSep'].iat[0]
    df = df.loc[(df['rangeSep'] == spatial_resolution)]

    # Filter out ground, low quality, and low power
    df = basic_SD_df_filter(df)

    df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id,
                                date_time_est=(df['datetime'].iat[0]).to_pydatetime())
    if hour_range is None:
        # Then infer hour range based on the available data
        hour_range = (round(df[time_units].iat[0], 2), round(df[time_units].iat[-1], 2))

    # Put frequencies in MHz and round, this makes them easier to compare
    df['transFreq_MHz'] = round(df['transFreq'] * 1e-3, 0)
    df = df.loc[np.isin(df['transFreq_MHz'], frequencies)]  # drop any rows with unrecognized frequencies
    df.reset_index(drop=True, inplace=True)

    date = df['datetime'].iat[0]
    date_string = calendar.month_name[date.month] + " " + str(date.day) + ", " + str(date.year)

    print("Preparing the plot...")
    resolution_string = str(spatial_resolution) + " km data"
    if beam_range[0] == beam_range[1]:
        beam_string = "Beam " + str(beam_range[0])
    else:
        beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    gate_string = "Gates " + str(gate_range[0]) + "-" + str(gate_range[1])

    fig, ax = plt.subplots(figsize=[8, 5], nrows=1, ncols=1, constrained_layout=True, dpi=300)

    if area is None:
        fig.suptitle(date_string + " at " + station.upper() + "; " + beam_string + "; " + gate_string
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)
    else:
        fig.suptitle(date_string + " at " + station.upper() + "; Area " + str(area)
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)

    format_plot(ax, x_lim=hour_range, y_lim=gate_range)

    for freq in frequencies:
        df_ff = df[(df['transFreq_MHz'] == freq)].copy()

        # Try to remove outliers
        df_ff = df_ff.loc[(df_ff['vel'] >= 100)]

        # Remove anything below the median and everything above the third quartile
        pwr = df_ff['pwr']  # Convenience
        pwr_median = np.median(pwr)
        pwr_upper_quartile = pwr[pwr > pwr_median]
        Q3 = np.median(pwr_upper_quartile)

        df_ff = df_ff.loc[(df_ff['pwr'] >= pwr_median) & (df_ff['pwr'] <= Q3)]

        points = np.vstack((df_ff[time_units], df_ff['gate'])).T

        # Build and plot the convex hull
        hull = ConvexHull(points)
        # _ = convex_hull_plot_2d(hull)

        for i, simplex in enumerate(hull.simplices):
            ax.plot(points[simplex, 0], points[simplex, 1], linewidth=1, color=hull_colours[freq])

            if i == 0:
                # Over-plot one to provide a legend handle
                ax.plot(points[simplex, 0], points[simplex, 1], linewidth=1, color=hull_colours[freq],
                        label=str(freq) + " MHz")

    ax.legend(loc='upper right')

    return fig


def format_plot(ax, x_lim, y_lim):
    """
    :param ax: matplotlib.axes:
            The axis to format
    :param x_lim: str:
            x label text
    :param y_lim: str:
            y label text
    :param t_diffs: dictionary of floats keyed by integer frequencies:
            The extra time delays to add in when adjusting elevation angles, in microseconds.
    """

    label_font_size = 12

    ax.set_ylim(y_lim)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(1))

    ax.set_xlim(x_lim)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax.set_ylabel("Range Gate", fontsize=label_font_size)
    ax.set_xlabel("UT Time [hour]", fontsize=label_font_size)

    ax.grid(b=True, which='major', axis='both', linewidth=1, linestyle='-')
    ax.grid(b=True, which='minor', axis='both', linewidth=0.4, linestyle='--')


if __name__ == "__main__":
    """ Testing """

    testing = True

    area = 2

    station = "rkn"
    year = "2016"
    month = "09"
    day = "26"
    start_hour = 0
    end_hour = 4

    beam_range = (7, 7)
    gate_range = (0, 74)

    # Read in SuperDARN data
    if area is None:
        loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
        in_dir = loc_root + "/DataReading/SD/data/" + station + "/" + station + year + month + day
        in_file = in_dir + "/" + station + year + month + day + ".pkl"
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        in_dir = loc_root + "/data"
        in_file = in_dir + "/" + station + year + month + day + "_area" + str(area) + ".pkl"
    print("Reading in file: " + in_file)
    df = pd.read_pickle(in_file)

    # Restrict data to within the desired hour range
    _, start_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=start_hour)
    _, end_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=end_hour)
    df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]

    fig = convex_hull_test(single_day_df=df, beam_range=beam_range, gate_range=gate_range,
                           hour_range=(start_hour, end_hour), area=area)

    if testing:
        plt.show()
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        if area is None:
            out_fig = out_dir + "/" + "watermann-" + station + year + month + day
        else:
            out_fig = out_dir + "/" + "watermann-" + station + year + month + day + "_area" + str(area)

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
