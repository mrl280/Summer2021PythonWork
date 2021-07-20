import calendar
import math
import os
import pathlib
import warnings
import bz2

import _pickle as cPickle
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch
from lib.basic_SD_df_filter import basic_SD_df_filter


def boxplot(single_day_df, beam_range, gate_range, area=None):
    """

    Make box-like plots with for a couple parameters of interest (frequency is along x)

    :param single_day_df: pandas.Data_Frame:
            A single days worth of data that we want to range profile
    :param gate_range: (<int>, <int>):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param area: int:
            The numbered area of interest.

    :return: matplotlib.pyplot.figure:
            The figure, it can then be viewed, modified, or saved to file
    """

    if len(single_day_df) <= 0:
        warnings.warn("There is no data in the provided dataframe", category=Warning)

    df = single_day_df.copy()  # Convenience
    frequencies = [10, 12, 13, 14]

    print("Filtering data..")
    # Make sure all of the data is of the same spatial resolution
    spatial_resolution = df['rangeSep'].iat[0]
    df = df.loc[(df['rangeSep'] == spatial_resolution)]

    # Filter out ground, low quality, and low power
    df = basic_SD_df_filter(df)

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

    fig, axes = plt.subplots(figsize=[12, 4], nrows=1, ncols=3, constrained_layout=True, dpi=300)

    if area is None:
        fig.suptitle(date_string + " at " + station.upper() + "; " + beam_string + "; " + gate_string
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)
    else:
        fig.suptitle(date_string + " at " + station.upper() + "; Area " + str(area)
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)

    format_plot(axes[0], labels={'x': "Frequency [MHz]", 'y': "Velocity [m/s]"})
    format_plot(axes[1], labels={'x': "Frequency [MHz]", 'y': "Spectral Width [m/s]"})
    format_plot(axes[2], labels={'x': "Frequency [MHz]", 'y': "Spectral Width [Hz]"})

    # fitACF widths are already in m/s
    df['wdt_ms'] = df['wdt']

    df['wdt_hz'] = (2 * math.pi * df['wdt']) * (df['transFreq'] / 1e5) / 3

    df = df.loc[((df['vel'] > -600) & (df['vel'] < 600)) &
                ((df['wdt'] > 0) & (df['wdt'] < 400))]

    add_data_to_plot(ax=axes[0], df=df, frequencies=np.asarray(frequencies), y_param='vel')
    add_data_to_plot(ax=axes[1], df=df, frequencies=np.asarray(frequencies), y_param='wdt_ms')
    add_data_to_plot(ax=axes[2], df=df, frequencies=np.asarray(frequencies), y_param='wdt_hz')

    return fig


def add_data_to_plot(ax, df, frequencies, y_param):
    """

    Add boxplot data.

    :param ax: matplotlib.axis:
            The axes to draw on
    :param df: pandas.DataFrame:
            The dataframe to use
    :param frequencies: numpy.array:
            List of frequencies in MHz
    :param y_param:
            The parameter that is on the y-axis
    """

    # Build up boxplot data at each frequency
    df_10 = df[(df['transFreq_MHz'] == 10)].copy()
    spread_10 = df_10[y_param]
    center10 = np.asarray([np.median(df_10[y_param])] * int(len(spread_10) / 2))
    data_10 = np.concatenate((spread_10, center10))

    df_12 = df[(df['transFreq_MHz'] == 12)].copy()
    spread_12 = df_12[y_param]
    center12 = np.asarray([np.median(df_12[y_param])] * int(len(spread_12) / 2))
    data_12 = np.concatenate((spread_12, center12))

    df_13 = df[(df['transFreq_MHz'] == 13)].copy()
    spread_13 = df_13[y_param]
    center13 = np.asarray([np.median(df_13[y_param])] * int(len(spread_13) / 2))
    data_13 = np.concatenate((spread_13, center13))

    df_14 = df[(df['transFreq_MHz'] == 14)].copy()
    spread_14 = df_14[y_param]
    center14 = np.asarray([np.median(df_14[y_param])] * int(len(spread_14) / 2))
    data_14 = np.concatenate((spread_14, center14))

    datasets = [data_10, data_12, data_13, data_14]
    ax.boxplot(datasets, labels=frequencies, showfliers=False, notch=True, showmeans=True, meanline=True)

    # Compute the medians ourself so we can draw a line between them
    medians = np.zeros(shape=frequencies.shape)
    medians[:] = np.nan

    for idx, freq in enumerate(frequencies):
        df_ff = df[(df['transFreq_MHz'] == freq)].copy()

        medians[idx] = np.median(df_ff[y_param])

    # box plots start at 1 on the x-axis, so we need to use special x values
    ax.plot([1, 2, 3, 4], medians, color='grey', linestyle='--', linewidth=1)


def format_plot(ax, labels):
    """

    :param ax: matplotlib.axes:
            The axis to format
    :param labels: dict:
            dictionary containing x and y labels, keyed with 'x' and 'y'
    """

    label_font_size = 12

    # ax.set_ylim(y_lim)
    # ax.yaxis.set_major_locator(MultipleLocator(200))
    # ax.yaxis.set_minor_locator(MultipleLocator(100))

    # ax.set_xlim(x_lim)
    ax.xaxis.set_major_locator(MultipleLocator(1))

    ax.set_xlabel(labels['x'], fontsize=label_font_size)
    ax.set_ylabel(labels['y'], fontsize=label_font_size)

    # ax.grid(b=True, which='major', axis='both', linewidth=1, linestyle='-')
    # ax.grid(b=True, which='minor', axis='both', linewidth=0.4, linestyle='--')


if __name__ == "__main__":
    """ Testing """

    testing = True

    area = None

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
        in_file = in_dir + "/" + station + year + month + day + ".pbz2"
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        in_dir = loc_root + "/data"
        in_file = in_dir + "/" + station + year + month + day + "_area" + str(area) + ".pbz2"
    print("Reading in file: " + in_file)
    data_stream = bz2.BZ2File(in_file, "rb")
    df = cPickle.load(data_stream)

    # Restrict data to within the desired hour range
    _, start_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=start_hour)
    _, end_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=end_hour)
    df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]

    fig = boxplot(single_day_df=df, gate_range=gate_range, beam_range=beam_range, area=area)

    if testing:
        plt.show()
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        if area is None:
            out_fig = out_dir + "/" + "boxplot-" + station + year + month + day
        else:
            out_fig = out_dir + "/" + "boxplot-" + station + year + month + day + "_area" + str(area)

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
