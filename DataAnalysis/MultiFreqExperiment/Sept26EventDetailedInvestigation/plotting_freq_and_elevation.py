import calendar
import os
import pathlib
import statistics
import warnings
import pydarn
import bz2

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from scipy import stats

from DataAnalysis.DataReading.SD.elevation_v2 import elevation_v2
from DataAnalysis.EchoOccurrence.lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch
from lib.basic_SD_df_filter import basic_SD_df_filter


def plotting_freq_and_elevation(single_day_df, beam_range, gate_range):
    """

    A program to investigate how elevation angles differ based on frequency, and how the required correctional t-diff
    may also need to be frequency dependant.

    :param single_day_df: pandas.Data_Frame:
            A single days worth of data that we want to range profile
    :param gate_range: (<int>, <int>):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :return: matplotlib.pyplot.figure:
            The figure, it can then be viewed, modified, or saved to file
    """

    if len(single_day_df) <= 0:
        warnings.warn("There is no data in the provided dataframe", category=Warning)

    df = single_day_df.copy()  # Convenience
    time_units = 'ut'
    frequencies = [10, 12, 13, 14]  # MHz
    t_diffs = [-0.003, 0.0, 0.001, 0.002, 0.003, 0.006]  # microseconds
    cmap = 'jet'
    zlim = (0, 35)

    print("Filtering data..")
    # Get some information about the station we are working with
    station = df['station'].iat[0]
    radar_id = pydarn.read_hdw_file(station).stid

    # Make sure all of the data is of the same spatial resolution
    spatial_resolution = df['rangeSep'].iat[0]
    df = df.loc[(df['rangeSep'] == spatial_resolution)]

    # Filter for beam and gate ranges of interest
    df = df.loc[((df['bmnum'] >= beam_range[0]) & (df['bmnum'] <= beam_range[1])) &
                ((df['gate'] >= gate_range[0]) & (df['gate'] <= gate_range[1]))]
    df.reset_index(drop=True, inplace=True)

    # Filter out ground scatter, low quality echoes, and low power echoes
    df = basic_SD_df_filter(df)

    df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id,
                                date_time_est=(df['datetime'].iat[0]).to_pydatetime())
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

    fig = plt.figure(figsize=[24, 10], constrained_layout=True, dpi=300)
    axes = add_axes(fig=fig)
    format_subplots(axes=axes, x_lim=hour_range, y_lim=gate_range, t_diffs=t_diffs)

    fig.suptitle(date_string + " at " + station.upper() + "; " + beam_string + "; " + gate_string +
                 "; " + resolution_string + "\n" + "Produced by " + str(os.path.basename(__file__)), fontsize=18)

    # Compute hour edges
    bins_per_hour = 60
    n_bins_x = int((hour_range[1] - hour_range[0]) * bins_per_hour)
    hour_edges = np.linspace(hour_range[0], hour_range[1], num=(n_bins_x + 1))

    # Compute gate_edges
    bins_per_gate = 1
    n_bins_y = int(((gate_range[1] + 1) - gate_range[0]) * bins_per_gate)
    gate_edges = np.linspace(gate_range[0], gate_range[1] + 1, num=(n_bins_y + 1), dtype=int)

    # Look though all of the subplots are profile the data
    for ax_idx, t_diff in enumerate(t_diffs):

        # Recompute elevation
        print("Recomputing Elevation Angles for t_diff=" + str(t_diff))
        elevation_v2(df=df, t_diff=t_diff)

        for freq in frequencies:

            df_ff = df[(df['transFreq_MHz'] == freq)]
            ax = axes[str(freq)][ax_idx]
            label = "t_diff=" + str(t_diff) + " \u03BCs"

            result, _, _, _ = stats.binned_statistic_2d(df_ff[time_units], df_ff['gate'], values=df_ff['adjElv'],
                                                        bins=[hour_edges, gate_edges])

            plot = ax.pcolormesh(hour_edges, gate_edges, result.transpose(),
                                 cmap=cmap, vmin=zlim[0], vmax=zlim[1], zorder=0)

            # ax.grid(b=True, which='both', axis='both', zorder=5)

    cbar_text_format = '%d'
    cbar = fig.colorbar(plot, cax=axes['color_bar'], orientation="horizontal", format=cbar_text_format, extend='max')
    cbar.ax.tick_params(labelsize=18)

    """"""
    # Just to see how extreme the difference in elevation angles is
    # elevation_v2(df=df, t_diff=0.003)
    # df['percent_diff'] = (df['adjElv'] - df['elv']) / df['elv']
    #
    # print(df[['elv', 'adjElv', 'percent_diff']])
    # print(statistics.mean(df['percent_diff']))
    """"""

    return fig


def add_axes(fig):
    """

    Create a sets of exes for the data.

    :param matplotlib.pyplot.figure: The figure to draw the axes on
    :return: Dictionary, matplotlib.pyplot.axis:
        A dictionary of data axes, and a colour bar axis
    """

    gs = fig.add_gridspec(ncols=18, nrows=9)

    # Remember that splices don't include last index
    axes = dict()

    axes["10"] = [fig.add_subplot(gs[0:2, 0:3]), fig.add_subplot(gs[0:2, 3:6]), fig.add_subplot(gs[0:2, 6:9]),
                  fig.add_subplot(gs[0:2, 9:12]), fig.add_subplot(gs[0:2, 12:15]), fig.add_subplot(gs[0:2, 15:18])]
    axes["12"] = [fig.add_subplot(gs[2:4, 0:3]), fig.add_subplot(gs[2:4, 3:6]), fig.add_subplot(gs[2:4, 6:9]),
                  fig.add_subplot(gs[2:4, 9:12]), fig.add_subplot(gs[2:4, 12:15]), fig.add_subplot(gs[2:4, 15:18])]
    axes["13"] = [fig.add_subplot(gs[4:6, 0:3]), fig.add_subplot(gs[4:6, 3:6]), fig.add_subplot(gs[4:6, 6:9]),
                  fig.add_subplot(gs[4:6, 9:12]), fig.add_subplot(gs[4:6, 12:15]), fig.add_subplot(gs[4:6, 15:18])]
    axes["14"] = [fig.add_subplot(gs[6:8, 0:3]), fig.add_subplot(gs[6:8, 3:6]), fig.add_subplot(gs[6:8, 6:9]),
                  fig.add_subplot(gs[6:8, 9:12]), fig.add_subplot(gs[6:8, 12:15]), fig.add_subplot(gs[6:8, 15:18])]
    axes['color_bar'] = fig.add_subplot(gs[8, 6:12])

    return axes


def format_subplots(axes, x_lim, y_lim, t_diffs):
    """

    :param axes: matplotlib.axes: The axes to format
    :param x_lim: str: x label text
    :param y_lim: str: y label text
    :param t_diffs: numpy.array of floats: The time delays to print out
    """

    label_font_size = 12
    title_font_size = 12

    for freq, ax_arr in axes.items():
        if ax_arr == axes['color_bar']:
            continue
        for ax_idx, ax in enumerate(ax_arr):

            ax.set_title(freq + " MHz Adjusted Elv, t_diff=" + str(t_diffs[ax_idx]) + " \u03BCs",
                         fontsize=title_font_size)

            if freq == "14":
                ax.set_xlabel("UT Time [hour]", fontsize=label_font_size)
            else:
                ax.set_xticklabels([])

            if ax != ax_arr[0]:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel("Range Gate", fontsize=label_font_size)

            ax.set_ylim(y_lim)
            ax.set_xlim(x_lim)


if __name__ == "__main__":
    """ Testing """

    testing = True

    station = "rkn"
    year = "2016"
    month = "09"
    day = "26"
    start_hour = 0
    end_hour = 4

    # year = "2017"
    # month = "10"
    # day = "23"
    # start_hour = 1
    # end_hour = 4

    # year = "2017"
    # month = "02"
    # day = "04"
    # start_hour = 4
    # end_hour = 7

    beam_range = (7, 7)
    gate_range = (0, 74)

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
    in_dir = loc_root + "/DataReading/SD/data/" + station + "/" + station + year + month + day
    in_file = in_dir + "/" + station + year + month + day + ".pbz2"
    data_stream = bz2.BZ2File(in_file, "rb")
    df = pd.read_pickle(data_stream)

    _, start_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=start_hour)
    _, end_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=end_hour)
    df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]

    fig = plotting_freq_and_elevation(single_day_df=df, beam_range=beam_range, gate_range=gate_range)

    if testing:
        plt.show()
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        out_fig = out_dir + "/" + "freq_and_elevation-" + station + year + month + day

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
