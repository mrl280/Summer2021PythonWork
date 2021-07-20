import calendar
import os
import pathlib
import statistics
import warnings
import bz2

import numpy as np
import pandas as pd
import _pickle as cPickle

from matplotlib import pyplot as plt

from DataAnalysis.DataReading.SD.elevation_v2 import elevation_v2
from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch
from lib.basic_SD_df_filter import basic_SD_df_filter


def percent_diff_elv_trial(single_day_df, beam_range, gate_range, area=None):
    """

    Find the optimal t_diff for multi-frequency events

    The idea: An artificial time difference is used to adjust elevation angles, and, for a given band of echoes,
    elevation angle should not increase with range.

    Adding in this correctional time difference causes elevation angles near the end of the bands to wrap back
    around. This wrapping results in a large percent difference change - maybe we can use this percent difference to
    figure out how much adjustment is required.

    TODO: Might need to restrict to echoes on the far range of the bands, this is where the wrapping happens

    :param single_day_df: pandas.Data_Frame:
            A SuperDARN dataframe
    :param beam_range: (<int>, <int>):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param gate_range: (<int>, <int>):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param area: int:
            The numbered area of interest.

    :return: matplotlib.pyplot.figure:
            The figure, it can then be viewed, modified, or saved to file
    """

    if len(single_day_df) <= 0:
        warnings.warn("There is no data in the provided dataframe", category=Warning)

    df = single_day_df.copy()  # Convenience
    frequencies = [10, 12, 13, 14]  # MHz
    colours = {'all': 'black',
               10: 'red',
               12: 'blue',
               13: 'green',
               14: 'magenta'}

    # Make sure all of the data is of the same spatial resolution
    spatial_resolution = df['rangeSep'].iat[0]
    df = df.loc[(df['rangeSep'] == spatial_resolution)]

    # Filter for beam and gate ranges of interest
    df = df.loc[((df['bmnum'] >= beam_range[0]) & (df['bmnum'] <= beam_range[1])) &
                ((df['gate'] >= gate_range[0]) & (df['gate'] <= gate_range[1]))]
    df.reset_index(drop=True, inplace=True)

    # Filter out ground scatter, low quality echoes, and low power echoes
    df = basic_SD_df_filter(df)

    # Put frequencies in MHz and round, this makes them easier to compare
    df['transFreq_MHz'] = round(df['transFreq'] * 1e-3, 0)
    df = df.loc[np.isin(df['transFreq_MHz'], frequencies)]  # drop any rows with unrecognized frequencies
    df.reset_index(drop=True, inplace=True)

    t_diff_start = -30  # Make sure we have lots of buffer because the ends will be useless
    t_diff_stop = 30
    num = int(t_diff_stop - t_diff_start + 1)
    tdiffs_in_ns = np.linspace(t_diff_start, t_diff_stop, num=num)  # Array of t_diffs in nanoseconds
    tdiffs_in_us = tdiffs_in_ns / 1000  # Array of t_diffs in microseconds

    date = df['datetime'].iat[0]
    date_string = calendar.month_name[date.month] + " " + str(date.day) + ", " + str(date.year)

    resolution_string = str(spatial_resolution) + " km data"
    if beam_range[0] == beam_range[1]:
        beam_string = "Beam " + str(beam_range[0])
    else:
        beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    gate_string = "Gates " + str(gate_range[0]) + "-" + str(gate_range[1])

    print("Preparing the plot...")
    fig, axes = plt.subplots(figsize=[10, 8], dpi=300, nrows=1, ncols=2, constrained_layout=True)
    if area is None:
        fig.suptitle(date_string + " at " + station.upper() + "; " + beam_string + "; " + gate_string
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)
    else:
        fig.suptitle(date_string + " at " + station.upper() + "; Area " + str(area)
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)

    axes[0].set_ylabel("Percent Change [%]")
    axes[1].set_ylabel("First Derivative of Percent Change")

    # Apply common subplot formatting
    for ax in axes:
        ax.grid(b=True, which='both', axis='both', zorder=5)

        ax.set_xlabel("Artificial Time Delay, t_diff [\u03BCs]")

    # Compute and plot data for all frequencies
    print(" Computing for all frequencies...")
    t_diff_adjustment = plot_percent_diff_data(axes=axes, df_ff=df, tdiffs_in_us=tdiffs_in_us, color=colours['all'],
                                               label="All Frequencies")

    print("Suggested t_diff adjustment for all frequency data is:" + str(t_diff_adjustment) + "microseconds")

    # Compute and plot data for each frequency
    for freq in frequencies:
        print(" Computing data for " + str(freq) + "MHz...")
        df_ff = df[(df['transFreq_MHz'] == freq)].copy()

        t_diff_adjustment = plot_percent_diff_data(axes=axes, df_ff=df_ff, tdiffs_in_us=tdiffs_in_us,
                                                   color=colours[freq], label=str(freq) + " MHz")
        print("Suggested t_diff adjustment for " + str(freq) + " MHz data is: " + str(t_diff_adjustment) +
              " microseconds")

    for ax in axes:
        ax.legend(loc='upper right')

    return fig


def plot_percent_diff_data(axes, df_ff, tdiffs_in_us, color, label):
    """


    :param df_ff: pandas.Data_Frame: The frequency restricted dataframe to use
    :param axes: The axes to draw on
    :param tdiffs_in_us: np.array: Array of t_diffs in microseconds
    :param color: str:
    :param label: str:
    """

    delta = tdiffs_in_us[1] - tdiffs_in_us[0]
    tdiffs_for_plotting = tdiffs_in_us[1:] - delta / 2

    # Compute elevation angles for each t_diff
    for t_diff in tdiffs_in_us:
        elevation_v2(df=df_ff, t_diff=t_diff)
        df_ff["elv_for_" + str(t_diff)] = df_ff['adjElv']

    # Compute the mean percent change between elevation angles resulting from adjacent t_diffs
    mean_percent_change = []
    for t in range(len(tdiffs_in_us) - 1):
        t_diff_init = tdiffs_in_us[t]
        t_diff_final = tdiffs_in_us[t + 1]
        elv_init_key = "elv_for_" + str(t_diff_init)
        elv_final_key = "elv_for_" + str(t_diff_final)

        mean_percent_change.append(statistics.mean(
            (df_ff[elv_final_key] - df_ff[elv_init_key]) / df_ff[elv_init_key]) * 100)

    # Lightly plot the raw data
    axes[0].plot(tdiffs_for_plotting, mean_percent_change, color=color, linewidth=0.2)

    # Smooth the signal
    kernel_size = 5
    kernel = np.ones(kernel_size) / kernel_size
    smoothed_mean_percent_change = np.convolve(mean_percent_change, kernel, mode='same')

    # Plot the smoothed signal
    axes[0].plot(tdiffs_for_plotting, smoothed_mean_percent_change, color=color, linewidth=1, label=label)

    # Compute the derivative of the smoothed signal
    mean_percent_change_diff = np.diff(smoothed_mean_percent_change)

    # The derivative is one element smaller
    tdiffs_for_plotting_diff = tdiffs_for_plotting[1:] - delta / 2

    # Lightly plot the raw data
    axes[1].plot(tdiffs_for_plotting_diff, mean_percent_change_diff, color=color, linewidth=0.2)

    # Smooth the derivative
    smoothed_mean_percent_change_diff = np.convolve(mean_percent_change_diff, kernel, mode='same')

    # Find the minimum of this derivative, it may be a good estimate of the required adjustment
    # Make sure we are not way off base be restricting the candidates
    likely_candidates = np.where((tdiffs_for_plotting_diff > -0.01) &
                                 (tdiffs_for_plotting_diff < 0.01))
    min_candidates_y = smoothed_mean_percent_change_diff[likely_candidates]
    min_candidates_x = tdiffs_for_plotting_diff[likely_candidates]

    min_y = np.min(min_candidates_y)
    x_index = np.where(min_candidates_y == min_y)
    min_x = min_candidates_x[x_index]

    # Plot the smoothed signal
    axes[1].plot(tdiffs_for_plotting_diff, smoothed_mean_percent_change_diff, color=color, linewidth=1,
                 label=label + " - min=" + str(min_x))

    return min_x


if __name__ == "__main__":
    """ Testing """

    testing = False

    station = "rkn"
    year = "2016"
    month = "09"
    day = "26"
    start_hour = 0
    end_hour = 4
    gate_range = (0, 74)
    area = 5

    # station = "rkn"
    # year = "2017"
    # month = "10"
    # day = "23"
    # start_hour = 1
    # end_hour = 4
    # gate_range = (0, 40)
    # area = None

    # station = "rkn"
    # year = "2017"
    # month = "02"
    # day = "04"
    # start_hour = 4
    # end_hour = 7
    # gate_range = (0, 74)
    # area = None

    beam_range = (7, 7)

    if area is not None:
        loc_root = str((pathlib.Path().parent.absolute()))
        in_dir = loc_root + "/data"
        in_file = in_dir + "/" + station + year + month + day + "_area" + str(area) + ".pbz2"
        print("Reading in file: " + in_file)
        data_stream = bz2.BZ2File(in_file, "rb")
        df = cPickle.load(data_stream)
    else:
        loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
        in_dir = loc_root + "/DataReading/SD/data/" + station + "/" + station + year + month + day
        in_file = in_dir + "/" + station + year + month + day + ".pbz2"
        print("Reading in file: " + in_file)
        data_stream = bz2.BZ2File(in_file, "rb")
        df = cPickle.load(data_stream)

        _, start_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=start_hour)
        _, end_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=end_hour)
        df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]

    fig = percent_diff_elv_trial(single_day_df=df, beam_range=beam_range, gate_range=gate_range, area=area)

    if testing:
        plt.show()
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        if area is None:
            out_fig = out_dir + "/" + "perc_diff_analysis-" + station + year + month + day
        else:
            out_fig = out_dir + "/" + "perc_diff_analysis-" + station + year + month + day + "_area" + str(area)

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
