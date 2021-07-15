import calendar
import os
import pathlib
import warnings
import pydarn

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy import stats

from DataAnalysis.DataReading.SD.elevation_v2 import elevation_v2
from DataAnalysis.EchoOccurrence.lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch
from lib.basic_SD_df_filter import basic_SD_df_filter


def simple_range_time_profiler(single_day_df, beam_range, gate_range, t_diffs, hour_range=None, area=None):
    """

    A simple single day range profiler for the multi-frequency events

    Plot velocity and elevation angle profiles for each frequency

    :param single_day_df: pandas.Data_Frame:
            A single days worth of data that we want to range profile
    :param gate_range: (<int>, <int>):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param t_diffs: dictionary of floats keyed by integer frequencies:
            The extra time delays to add in when adjusting elevation angles, in microseconds.
    :param hour_range: (<int>, <int>) (optional):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not quite inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param area: int:
            The numbered area of interest.

    :return: matplotlib.pyplot.figure:
            The figure, it can then be viewed, modified, or saved to file
    """

    if len(single_day_df) <= 0:
        warnings.warn("There is no data in the provided dataframe", category=Warning)

    df = single_day_df.copy()  # Convenience
    time_units = 'ut'
    frequencies = [10, 12, 13, 14]
    subplot_types = ["vel", "adjElv"]
    color_maps = {'vel': 'seismic_r',
                  'adjElv': 'jet'}
    zlims = {'vel': (-600, 600),
             'adjElv': (0, 30)}

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
    # t_diff_string = "tdiff: " + str(t_diff) + " \u03BCs"

    fig = plt.figure(figsize=[12, 12], constrained_layout=True, dpi=300)
    data_axes = add_axes(fig=fig)
    format_subplots(axes=data_axes, x_lim=hour_range, y_lim=gate_range, t_diffs=t_diffs)

    if area is None:
        fig.suptitle(date_string + " at " + station.upper() + "; " + beam_string + "; " + gate_string
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)
    else:
        fig.suptitle(date_string + " at " + station.upper() + "; Area " + str(area)
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)

    print("Computing and plotting binned occurrence rates...")

    # Compute hour edges
    bins_per_hour = 60
    n_bins_x = int((hour_range[1] - hour_range[0]) * bins_per_hour)
    hour_edges = np.linspace(hour_range[0], hour_range[1], num=(n_bins_x + 1))

    # Compute gate_edges
    bins_per_gate = 1
    n_bins_y = int(((gate_range[1] + 1) - gate_range[0]) * bins_per_gate)
    gate_edges = np.linspace(gate_range[0], gate_range[1] + 1, num=(n_bins_y + 1), dtype=int)

    # Look though all of the frequencies and plot the data
    for freq in frequencies:
        for param in subplot_types:
            df_ff = df[(df['transFreq_MHz'] == freq)].copy()
            ax = data_axes[freq][param]

            if len(df_ff) <= 0:
                continue
            else:
                if param == 'adjElv':
                    # Recompute elevation with the extra t_diff
                    print("Recomputing Elevation Angles for " + str(freq) + " MHz data - t_diff=" + str(t_diffs[freq]))
                    elevation_v2(df=df_ff, t_diff=t_diffs[freq])  # t_diff is in microseconds

                result, _, _, _ = stats.binned_statistic_2d(df_ff[time_units], df_ff['gate'], values=df_ff[param],
                                                            bins=[hour_edges, gate_edges])

                plot = ax.pcolormesh(hour_edges, gate_edges, result.transpose(),
                                     cmap=color_maps[param], vmin=zlims[param][0], vmax=zlims[param][1], zorder=0)

                cbar_text_format = '%d'
                if param == 'vel':
                    cbar = fig.colorbar(plot, ax=ax, orientation="vertical", format=cbar_text_format, extend='both')
                else:
                    cbar = fig.colorbar(plot, ax=ax, orientation="vertical", format=cbar_text_format, extend='max')

            ax.grid(b=True, which='major', axis='both', linewidth=1, linestyle='-')
            ax.grid(b=True, which='minor', axis='both', linewidth=0.4, linestyle='--')

    return fig


def add_axes(fig):
    """

    Create a sets of exes for the data.

    :param matplotlib.pyplot.figure: The figure to draw the axes on
    :return: Dictionary, matplotlib.pyplot.axis:
        A dictionary of data axes, and a colour bar axis
    """

    gs = fig.add_gridspec(ncols=7, nrows=9)

    # Remember that splices don't include last index
    data_axes = dict()
    data_axes[10] = {"vel": fig.add_subplot(gs[0:2, 0:3]), "adjElv": fig.add_subplot(gs[2:4, 0:3])}
    data_axes[12] = {"vel": fig.add_subplot(gs[5:7, 0:3]), "adjElv": fig.add_subplot(gs[7:9, 0:3])}
    data_axes[13] = {"vel": fig.add_subplot(gs[0:2, 4:7]), "adjElv": fig.add_subplot(gs[2:4, 4:7])}
    data_axes[14] = {"vel": fig.add_subplot(gs[5:7, 4:7]), "adjElv": fig.add_subplot(gs[7:9, 4:7])}

    return data_axes


def format_subplots(axes, x_lim, y_lim, t_diffs):
    """
    :param axes: matplotlib.axes:
            The axes to format
    :param x_lim: str:
            x label text
    :param y_lim: str:
            y label text
    :param t_diffs: dictionary of floats keyed by integer frequencies:
            The extra time delays to add in when adjusting elevation angles, in microseconds.
    """

    subplot_types = ["vel", "adjElv"]
    label_font_size = 14
    title_font_size = 14

    for freq, axis in axes.items():

        for subplot_type in subplot_types:
            ax = axis[subplot_type]
            ax.tick_params(labeltop=True, labelright=True)

            ax.set_ylim(y_lim)
            ax.yaxis.set_major_locator(MultipleLocator(10))
            ax.yaxis.set_minor_locator(MultipleLocator(1))

            ax.set_xlim(x_lim)
            ax.xaxis.set_major_locator(MultipleLocator(0.5))
            ax.xaxis.set_minor_locator(MultipleLocator(0.1))

            ax.set_ylabel("Range Gate", fontsize=label_font_size)

            if subplot_type == "vel":
                ax.set_title(str(freq) + " MHz Velocities", fontsize=title_font_size)
            elif subplot_type == "adjElv":
                ax.set_title(str(freq) + " MHz Adjusted Elevation Angles, tdiff=" + str(t_diffs[freq]) + " \u03BCs",
                             fontsize=title_font_size)
                ax.set_xlabel("UT Time [hour]", fontsize=label_font_size)


if __name__ == "__main__":
    """ Testing """

    testing = False

    area = None

    station = "rkn"
    year = "2016"
    month = "09"
    day = "26"
    start_hour = 0
    end_hour = 4

    beam_range = (7, 7)
    gate_range = (0, 74)
    t_diffs = {10: 0.003,  # microseconds
               12: 0.003,
               13: 0.003,
               14: 0.003}

    # t_diffs = {10: -0.003,  # microseconds
    #            12: 0.001,
    #            13: 0.002,
    #            14: 0.003}

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

    fig = simple_range_time_profiler(single_day_df=df, beam_range=beam_range, gate_range=gate_range,
                                     hour_range=(start_hour, end_hour), t_diffs=t_diffs, area=area)

    if testing:
        plt.show()
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        if area is None:
            out_fig = out_dir + "/" + "simple_range_time_profile-" + station + year + month + day
        else:
            out_fig = out_dir + "/" + "simple_range_time_profile-" + station + year + month + day + "_area" + str(area)

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
