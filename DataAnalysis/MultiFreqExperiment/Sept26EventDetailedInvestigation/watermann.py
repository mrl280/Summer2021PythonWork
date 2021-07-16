import calendar
import os
import pathlib
import warnings
import pydarn

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, NullFormatter
from scipy import stats

from DataAnalysis.DataReading.SD.elevation_v2 import elevation_v2
from DataAnalysis.EchoOccurrence.lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch
from lib.basic_SD_df_filter import basic_SD_df_filter


def watermann(single_day_df, axes_vlimits, area=None):
    """

    Produce three plots:
        - normalized echo count vs doppler histogram
        - width vs doppler
        - power vs doppler

    :param single_day_df: pandas.Data_Frame:
            A single days worth of data that we want to range profile
    :param axes_vlimits: dict:
            Dictionary of vlim (zlim) ranges keyed by plot type.
    :param area: int:
            The numbered area of interest.

    :return: matplotlib.pyplot.figure:
            The figure, it can then be viewed, modified, or saved to file
    """

    if len(single_day_df) <= 0:
        warnings.warn("There is no data in the provided dataframe", category=Warning)

    df = single_day_df.copy()  # Convenience
    frequencies = [10, 12, 13, 14]
    axes_limits = {'vel': (-600, 600),
                   'count': (0, 1.0),
                   'wdt': (0, 200),
                   'pwr': (0, 50)}
    axes_labels = {'vel': "Velocities [m/s]",
                   'count': "Normalized Echo Count",
                   'wdt': "Spectral Width [m/s]",
                   'pwr': "Power [dB]"}
    cbar_labelsize = 12

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

    # Put frequencies in MHz and round, this makes them easier to access them
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

    fig = plt.figure(figsize=[12, 8], constrained_layout=True, dpi=300)
    axes = add_axes(fig=fig)
    format_subplots(axes=axes, axes_limits=axes_limits, axes_labels=axes_labels)

    if area is None:
        fig.suptitle(date_string + " at " + station.upper() + "; " + beam_string + "; " + gate_string
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)
    else:
        fig.suptitle(date_string + " at " + station.upper() + "; Area " + str(area)
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)

    for freq in frequencies:
        df_ff = df[(df['transFreq_MHz'] == freq)].copy()

        add_in_count_data(ax=axes[freq]['count'], df_ff=df_ff)

        # Add width contours
        wdt_reference_plot = add_in_wdt_data(ax=axes[freq]['wdt'], df_ff=df_ff, vlims=axes_vlimits['wdt'])
        cbar_wdt = fig.colorbar(wdt_reference_plot, cax=axes['cbar']['wdt'], shrink=0.75, format='%d')
        cbar_wdt.ax.tick_params(labelsize=cbar_labelsize)

        # Add power contours
        pwr_reference_plot = add_in_pwr_data(ax=axes[freq]['pwr'], df_ff=df_ff, vlims=axes_vlimits['pwr'])
        cbar_pwr = fig.colorbar(pwr_reference_plot, cax=axes['cbar']['pwr'], shrink=0.75, format='%d')
        cbar_pwr.ax.tick_params(labelsize=cbar_labelsize)

        # Add in gridlines
        for keys, ax in axes[freq].items():
            ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.7, zorder=4)

    return fig


def add_in_pwr_data(ax, df_ff, vlims):
    """
    Plot width vs doppler

    :param vlims: (int, int):
            Colour bar limits
    :param: fig: matplotlib.pyplot.figure:
            The figure to use
    :param ax: matplotlib.pyplot.axis:
            The axis to draw on
    :param df_ff: pandas.DataFrame:
            The frequency restricted dataframe to use

    :return plot: matplotlib.pyplot.plot:
            An example plot that can be referenced when adding the colour bar
    """

    if len(df_ff) <= 0:
        return

    vel_range = ax.get_xlim()
    pwr_range = ax.get_ylim()

    n_levels = 12
    levels = np.linspace(start=vlims[0], stop=vlims[1], num=(n_levels + 1))

    # Compute velocity edges
    vel_bin_width = 10  # m/s
    n_bins_vel = int((vel_range[1] - vel_range[0]) / vel_bin_width)
    vel_edges = np.linspace(vel_range[0], vel_range[1], num=(n_bins_vel + 1))

    # Compute power edges
    pwr_bin_width = 2
    n_bins_y = int((pwr_range[1] - pwr_range[0]) / pwr_bin_width)
    pwr_edges = np.linspace(pwr_range[0], pwr_range[1], num=(n_bins_y + 1))

    # Compute bin centers
    delta_vel = vel_edges[1] - vel_edges[0]
    delta_pwr = pwr_edges[1] - pwr_edges[0]
    binned_xcenters = vel_edges[1:] - delta_vel / 2
    binned_ycenters = pwr_edges[1:] - delta_pwr / 2

    result, _, _, _ = stats.binned_statistic_2d(df_ff['vel'], df_ff['pwr'], values=None, statistic='count',
                                                bins=[vel_edges, pwr_edges])

    # Replace the zeros with nans so no data shows up as white
    result[result == 0] = np.nan

    plot = ax.contourf(binned_xcenters, binned_ycenters, result.transpose(), cmap='jet', levels=levels, extend='max')

    # Mark our 3 dB power cutoff
    ax.axhline(y=3, color='black', linestyle='-', linewidth=1, zorder=3)

    return plot


def add_in_wdt_data(ax, df_ff, vlims):
    """
    Plot width vs doppler

    :param vlims: (int, int):
            Colour bar limits
    :param: fig: matplotlib.pyplot.figure:
            The figure to use
    :param ax: matplotlib.pyplot.axis:
            The axis to draw on
    :param df_ff: pandas.DataFrame:
            The frequency restricted dataframe to use

    :return plot: matplotlib.pyplot.plot:
            An example plot that can be referenced when adding the colour bar
    """

    if len(df_ff) <= 0:
        return

    vel_range = ax.get_xlim()
    wdt_range = ax.get_ylim()

    n_levels = 12
    levels = np.linspace(start=vlims[0], stop=vlims[1], num=(n_levels + 1))

    # Compute velocity edges
    vel_bin_width = 10  # m/s
    n_bins_vel = int((vel_range[1] - vel_range[0]) / vel_bin_width)
    vel_edges = np.linspace(vel_range[0], vel_range[1], num=(n_bins_vel + 1))

    # Compute width edges
    wdt_bin_width = 5
    n_bins_y = int((wdt_range[1] - wdt_range[0]) / wdt_bin_width)
    wdt_edges = np.linspace(wdt_range[0], wdt_range[1], num=(n_bins_y + 1))

    # Compute bin centers
    delta_vel = vel_edges[1] - vel_edges[0]
    delta_wdt = wdt_edges[1] - wdt_edges[0]
    binned_xcenters = vel_edges[1:] - delta_vel / 2
    binned_ycenters = wdt_edges[1:] - delta_wdt / 2

    result, _, _, _ = stats.binned_statistic_2d(df_ff['vel'], df_ff['wdt'], values=None, statistic='count',
                                                bins=[vel_edges, wdt_edges])

    # Replace the zeros with nans so no data shows up as white
    result[result == 0] = np.nan

    plot = ax.contourf(binned_xcenters, binned_ycenters, result.transpose(), cmap='jet', levels=levels, extend='max')

    return plot


def add_in_count_data(ax, df_ff):
    """

    Plot normalized echo count vs doppler

    :param ax: matplotlib.axes:
            The axes to draw or
    :param df_ff: pandas.DataFrame:
            The frequency restricted dataframe to use
    """

    pass


def add_axes(fig):
    """

    Create a sets of exes for the data.

    :param matplotlib.pyplot.figure:
        The figure to draw the axes on

    :return axes: Dictionary of matplotlib.pyplot.axis:
        A dictionary of data axes
    """

    gs = fig.add_gridspec(nrows=3, ncols=25)

    axes = dict()
    axes[10] = {"count": fig.add_subplot(gs[0, 0:6]),
                "wdt": fig.add_subplot(gs[1, 0:6]),
                "pwr": fig.add_subplot(gs[2, 0:6])}
    axes[12] = {"count": fig.add_subplot(gs[0, 6:12]),
                "wdt": fig.add_subplot(gs[1, 6:12]),
                "pwr": fig.add_subplot(gs[2, 6:12])}
    axes[13] = {"count": fig.add_subplot(gs[0, 12:18]),
                "wdt": fig.add_subplot(gs[1, 12:18]),
                "pwr": fig.add_subplot(gs[2, 12:18])}
    axes[14] = {"count": fig.add_subplot(gs[0, 18:24]),
                "wdt": fig.add_subplot(gs[1, 18:24]),
                "pwr": fig.add_subplot(gs[2, 18:24])}

    # Leave some axes for the colour bars
    axes['cbar'] = {"wdt": fig.add_subplot(gs[1, 24]),
                    "pwr": fig.add_subplot(gs[2, 24])}

    return axes


def format_subplots(axes, axes_limits, axes_labels):
    """
    :param axes: matplotlib.axes:
            The axes to format
    :param axes_limits: Dictionary of ranges keyed by parameter:
            The axes limits
    """

    subplot_types = ["count", "wdt", "pwr"]
    label_font_size = 14
    title_font_size = 16

    x_lim = axes_limits['vel']
    x_label = axes_labels['vel']

    for freq, axis in axes.items():
        if axis == axes['cbar']:
            continue

        for subplot_type in subplot_types:
            ax = axis[subplot_type]

            y_lim = axes_limits[subplot_type]
            ax.set_ylim(y_lim)
            # ax.yaxis.set_major_locator(MultipleLocator(10))
            # ax.yaxis.set_minor_locator(MultipleLocator(1))

            ax.set_xlim(x_lim)
            ax.xaxis.set_major_locator(MultipleLocator(200))
            ax.xaxis.set_minor_locator(MultipleLocator(100))

            # Turn off labels on all but the left border plots
            if freq == 10:
                y_label = axes_labels[subplot_type]
                ax.set_ylabel(y_label, fontsize=label_font_size)
            else:
                ax.yaxis.set_major_formatter(NullFormatter())

            # Turn off labels on all but the bottom border plots
            if subplot_type == 'pwr':
                ax.set_xlabel(x_label, fontsize=label_font_size)
            else:
                ax.xaxis.set_major_formatter(NullFormatter())

            if subplot_type == 'count':
                ax.set_title(str(freq) + " MHz", fontsize=title_font_size)


if __name__ == "__main__":
    """ Testing """

    testing = True

    area = None
    axes_vlimits = {'wdt': (0, 25),
                    'pwr': (0, 50)}

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

    fig = watermann(single_day_df=df, axes_vlimits=axes_vlimits, area=area)

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
