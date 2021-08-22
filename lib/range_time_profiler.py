import bz2
import pathlib

import pandas as pd
import pydarn
import numpy as np
from matplotlib import pyplot as plt

from scipy import stats

from DataAnalysis.DataReading.SD.elevation_v2 import elevation_v2
from DataAnalysis.EchoOccurrence.lib.add_decimal_hour_to_df import add_decimal_hour_to_df


def range_time_profiler(df, gate_range, hour_range, param, time_units='ut'):
    """

    A simple range time profiler for a single parameter

    :param df: pandas.DataFrame:
            A SuperDARN fit dataframe.
    :param gate_range: (int, int):
            The gate range to profile (y-axis limits)
    :param hour_range: (float, float):
            The hour range to profile (x-axis limits)
    :param param: str:
            The fit parameter to range-time profile

    :param time_units: str (optional; default is 'ut')
            The time units to use.
                'ut' for universal time
                'mlt' for magnetic local time
                'lt' for local time (based on longitude)
                'lst' for local standard time (based on time zones)
                'ast' for apparent solar time (based on the apparent angular motion of the sun across the sky)

    :return fig: matplotlib.pyplot.figure:

    """

    if len(df) <= 0:
        raise Exception("range_time_profiler(): Empty dataframe.")

    # Add UT decimal hour to df
    station = df['station'].iat[0]
    df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=pydarn.read_hdw_file(station).stid,
                                date_time_est=df['datetime'].iat[0])
    df = df.loc[(df[time_units] > hour_range[0]) & (df[time_units] < hour_range[1]) &
                (df['gate'] > gate_range[0]) & (df['gate'] < gate_range[1])]

    # Compute hour edges
    bins_per_hour = 60
    n_bins_x = int((hour_range[1] - hour_range[0]) * bins_per_hour)
    hour_edges = np.linspace(hour_range[0], hour_range[1], num=(n_bins_x + 1))

    # Compute gate_edges
    bins_per_gate = 1
    n_bins_y = int(((gate_range[1] + 1) - gate_range[0]) * bins_per_gate)
    gate_edges = np.linspace(gate_range[0], gate_range[1] + 1, num=(n_bins_y + 1), dtype=int)

    result, _, _, _ = stats.binned_statistic_2d(df[time_units], df['gate'], values=df[param],
                                                bins=[hour_edges, gate_edges])

    fig, ax = plt.subplots(figsize=(8, 4), dpi=300, constrained_layout=True, nrows=1, ncols=1)
    ax.set_xlabel("Hour [" + time_units.upper() + "]")
    ax.set_xlabel("Range Gate")
    plot = ax.pcolormesh(hour_edges, gate_edges, result.transpose(), cmap='jet', vmin=0, vmax=40,
                         zorder=0)

    cbar_text_format = '%d'
    if param == 'vel':
        cbar = fig.colorbar(plot, ax=ax, orientation="vertical", format=cbar_text_format, extend='both')
    else:
        cbar = fig.colorbar(plot, ax=ax, orientation="vertical", format=cbar_text_format, extend='max')

    return fig


if __name__ == "__main__":
    """ """

    station = "hok"
    year = "2014"
    month = "02"
    day = "23"
    t_diff = 0  # in nanoseconds
    gate_range = (0, 74)
    hour_range = (6, 12)
    beam = 7

    time_units = 'ut'
    start_hour = 6
    end_hour = 12

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()))
    in_dir = loc_root + "/DataAnalysis/DataReading/SD/data/" + station + "/" + station + year + month + day
    in_file = in_dir + "/" + station + year + month + day + ".pbz2"
    data_stream = bz2.BZ2File(in_file, "rb")
    df = pd.read_pickle(data_stream)

    df = df.loc[df['bmnum'] == beam]
    df.reset_index(drop=True, inplace=True)


    print("Recomputing elevation angles...")
    elevation_v2(df=df, t_diff=t_diff / 1000)

    fig = range_time_profiler(df=df, gate_range=gate_range, hour_range=hour_range, param='adjElv', time_units='ut')
    plt.show()
