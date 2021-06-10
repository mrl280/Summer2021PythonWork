import pathlib
import pydarn

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import radar_fov, SuperDARNRadars
from scipy import stats

from lib.add_mlt_to_df import add_mlt_to_df
from lib.cm.modified_viridis import modified_viridis_2
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.z_min_max_defaults import z_min_max_defaults
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def occ_year_vs_ut(station, year_range, month_range=None, time_units='mlt', hour_range=None, gate_range=None, beam_range=None,
                   parameter=None, local_testing=False):
    """

    Produce a contour plot with year on the y-axis and time along the x-axis.
    Plots a simple echo count

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only considers 45 km data
        - Does not distinguish frequency
        - This program uses fitACF 3.0 data.  To change this, modify the source code.
        - year_range is assumed UT, hour_range is either MLT of UT depending on time_units

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year_range: (<int>, <int>):
            Inclusive. The year range to consider.
    :param month_range: (<int>, <int>) (optional):
            Inclusive. The months of the year to consider.  If omitted (or None), then all days will be considered.
    :param time_units: str: 'ut' for universal time or 'mlt' for magnetic local time:
            The time units to plot along x.  Default is 'mlt'
    :param hour_range: (<int>, <int>) (optional):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not quite inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param gate_range: (<int>, <int>) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param parameter: str (optional):
            Parameter to be averaged (e.g. 'v' or 'p_l')
            If omitted, then a simple echo count will be plotted.
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :return: matplotlib.pyplot.figure
            The figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    time_units = check_time_units(time_units)
    hour_range = check_hour_range(hour_range)

    print("Retrieving data...")
    # TODO: Decide what to do about month and day ranges
    df = get_data_handler(station, year_range=year_range, month_range=month_range, day_range=(1, 31),
                          hour_range=hour_range, gate_range=gate_range, beam_range=beam_range,
                          local_testing=local_testing)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid

    print("Filtering data...")
    df = df.loc[(df['p_l'] >= 3)]  # Restrict to points with at least 3 dB
    if parameter is not None:
        zmin, zmax = z_min_max_defaults(parameter)
        df = df.loc[(df[parameter] >= zmin) & (df[parameter] <= zmax)]

        if parameter == 'v':
            cmap = 'seismic_r'
            levels = np.linspace(zmin, zmax, 13, endpoint=True)
        else:
            cmap = modified_viridis_2()
            levels = 6
    else:
        cmap = modified_viridis_2()
        levels = 12

    df.reset_index(drop=True, inplace=True)
    df = only_keep_45km_res_data(df)

    print("Preparing the plot...")
    n_rows = (year_range[1] + 1) - year_range[0]  # We will have one subplot for every year of observations
    fig, ax = plt.subplots(figsize=(8, 9), sharex='col', dpi=300, nrows=n_rows, ncols=1)
    plt.subplots_adjust(hspace=0.05)

    # Apply common subplot formatting
    ax[n_rows - 1].set_xlabel("Time, " + time_units.upper(), fontsize=18)
    for row in reversed(range(n_rows)):
        ax[row].set_ylim([0, 12])
        ax[row].set_xlim(hour_range)
        ax[row].tick_params(axis='y', which='major', labelleft=False, direction='in')
        ax[row].tick_params(axis='x', which='major', direction='in')
        plt.xticks(fontsize=16)
        ax[row].xaxis.set_major_locator(MultipleLocator(2))  # Every 2 hours
        ax[row].yaxis.set_major_locator(MultipleLocator(6.5))  # Half way through the year
        ax[row].grid(b=True, which='major', axis='x', linestyle='--', linewidth=0.5, zorder=4)
        ax[row].grid(b=True, which='major', axis='y', linestyle='--', linewidth=0.5, zorder=4)
        ax[row].set_ylabel(year_range[1] - row, fontsize=18)

    # Loop through the years, adding in the data
    n_bins_y = 48  # half month bins
    n_bins_x = (hour_range[1] - hour_range[0]) * 2  # half hour bins
    contour_range = [ax[0].get_xlim(), ax[0].get_ylim()]  # All plots are the same size
    for row in reversed(range(n_rows)):
        year = year_range[1] - row

        # Build a restricted dataframe with the year's data
        start_datetime, start_epoch = build_datetime_epoch(year=year, month=1, day=1, hour=0)
        end_datetime, end_epoch = build_datetime_epoch(year=year, month=12, day=31, hour=24)
        df_yy = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)].copy()
        df_yy.reset_index(drop=True, inplace=True)

        df_yy['ut_time'] = df_yy['hour'] + df_yy['minute'] / 60 + df_yy['second'] / 3600

        if time_units == "mlt":
            print("Computing MLTs for " + str(year) + " data...")

            # To compute mlt we need longitudes..
            # we will use the middle of the year and assume magnetic longitudes don't change much over the year
            mid_datetime = start_datetime + (end_datetime - start_datetime) / 2
            cell_corners_aacgm_lats, cell_corners_aacgm_lons = \
                radar_fov(stid=radar_id, coords='aacgm', date=mid_datetime)

            df_yy = add_mlt_to_df(cell_corners_aacgm_lons=cell_corners_aacgm_lons,
                                  cell_corners_aacgm_lats=cell_corners_aacgm_lats, df=df_yy)

            df_yy['xdata'] = df_yy['mlt']
        else:
            df_yy['xdata'] = df_yy['ut_time']

        df_yy = df_yy.loc[(df_yy['xdata'] >= hour_range[0]) & (df_yy['xdata'] <= hour_range[1])]
        df_yy.reset_index(drop=True, inplace=True)

        # Compute decimal datetime to plot along y
        df_yy['ydata'] = (df_yy['month'] - 1) + (df_yy['day'] - 1) / 31 + df_yy['ut_time'] / 730

        if parameter is None:
            # We just want a simple echo count
            binned_counts, bin_xedges, bin_yedges, bin_numbers = stats.binned_statistic_2d(
                df_yy['xdata'], df_yy['ydata'], values=None,
                statistic='count', bins=[n_bins_x, n_bins_y], range=contour_range)
        else:
            # We want to the median of the chosen parameter
            binned_counts, bin_xedges, bin_yedges, bin_numbers = stats.binned_statistic_2d(
                df_yy['xdata'], df_yy['ydata'], values=df_yy[parameter],
                statistic='median', bins=[n_bins_x, n_bins_y], range=contour_range)
            binned_counts = np.nan_to_num(binned_counts)

        # Compute bin centers
        bin_xwidth = (bin_xedges[1] - bin_xedges[0])
        bin_ywidth = (bin_yedges[1] - bin_yedges[0])
        bin_xcenters = bin_xedges[1:] - bin_xwidth / 2
        bin_ycenters = bin_yedges[1:] - bin_ywidth / 2

        # Plot the data
        cont = ax[row].contourf(bin_xcenters, bin_ycenters, binned_counts.transpose(),
                                cmap=cmap, levels=levels)
        cbar = fig.colorbar(cont, ax=ax[row], shrink=0.75)
        cbar.ax.tick_params(labelsize=16)

    print("Returning the figure...")
    return fig


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    station = "rkn"
    fig = occ_year_vs_ut(station=station, time_units='mlt', year_range=(2011, 2012), month_range=(1, 12),
                         gate_range=(30, 74), beam_range=(7, 7),
                         parameter='v', local_testing=local_testing)

    if local_testing:
        a = 1
        plt.show()
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        out_file = out_dir + "/occ_year_vs_time_" + station
        print("Saving plot as " + out_file)
        fig.savefig(out_file + ".jpg", format='jpg', dpi=300)
