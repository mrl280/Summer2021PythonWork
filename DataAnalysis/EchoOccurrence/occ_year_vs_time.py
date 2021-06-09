import pathlib

import aacgmv2
import numpy as np
import pandas as pd
import pydarn

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import radar_fov
from scipy import stats

from lib.z_min_max_defaults import z_min_max_defaults
from lib.centroid import centroid
from lib.build_date_epoch import build_date_epoch
from lib.data_getters.get_data import get_data
from lib.data_getters.get_local_dummy_data import get_local_dummy_data
from lib.range_checkers import *


def occ_year_vs_ut(station, year_range, time_units='mlt', hour_range=None, gate_range=None, beam_range=None,
                   local_testing=False, parameter=None):
    """

    Produce a contour plot with year on the y-axis and time along the x-axis.

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only considers 45 km data
        - Does not distinguish frequency
        - This program uses fitACF 3.0 data.  To change this, modify the source code.
        - All times and dates are assumed UT

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year_range: (<int>, <int>):
            Inclusive. The year range to consider.
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
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :param parameter: str (optional):
            Parameter to be averaged (e.g. 'v' or 'p_l')
            If omitted, then a simple echo count will be plotted.
    :return: matplotlib.pyplot.figure
            The figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    if time_units.upper() == "MLT":
        time_units = "mlt"
    elif time_units.upper() == "UT":
        time_units = "ut"
    else:
        raise Exception("invalid time units passed to occ_year_vs_ut().  Options are 'mlt' or 'ut'.")

    if parameter is not None:
        zmin, zmax = z_min_max_defaults(parameter)

    year_range = check_year_range(year_range)
    hour_range = check_hour_range(hour_range)
    month_range = (1, 12)
    day_range = (1, 31)  # TODO: fix

    if isinstance(station, str):
        hdw_info = pydarn.read_hdw_file(station)  # Get the hardware file, there is lots of good stuff in there
    else:
        raise Exception("Error: Please enter the station as a 3 character string.  e.g. 'rkn'")

    beam_range = check_beam_range(beam_range, hdw_info)
    gate_range = check_gate_range(gate_range, hdw_info)

    print("Retrieving data...")
    if local_testing:
        # Just read in some test data
        warnings.warn("Running in local testing mode, just going to use local dummy data", category=Warning)
        # df = get_local_dummy_data(station=station, year=2011, month=9, day=29, start_hour_UT=0, end_hour_UT=23)
        df = get_local_dummy_data(station=station, year=2011, month=11, day=12, start_hour_UT=0, end_hour_UT=2)
        df_2 = get_local_dummy_data(station=station, year=2011, month=9, day=29, start_hour_UT=0, end_hour_UT=24)
        df_3 = get_local_dummy_data(station=station, year=2012, month=10, day=15, start_hour_UT=0, end_hour_UT=24)
        df = pd.concat([df, df_2, df_3])
        df = df.loc[(df['bmnum'] >= beam_range[0]) & (df['bmnum'] <= beam_range[1]) &
                    (df['slist'] >= gate_range[0]) & (df['slist'] <= gate_range[1])]
        # print(df.keys())
    else:
        df = get_data(station=station, year_range=year_range, month_range=month_range, day_range=day_range,
                      hour_range=hour_range, gate_range=gate_range, beam_range=beam_range)

    df = df.loc[(df['p_l'] >= 3)]  # Restrict to points with at least 3 dB
    if parameter is not None:
        df = df.loc[(df[parameter] >= zmin) & (df[parameter] <= zmax)]
    df.reset_index(drop=True, inplace=True)

    # Restrict to 45 km mode data, you need to modify the fan if you want to use data of a different spatial resolution
    number_of_records_before_spatial_resolution_check = df.shape[0]
    df = df.loc[(df['frang'] == 180) & (df['rsep'] == 45)]
    df.reset_index(drop=True, inplace=True)
    number_of_records_after_spatial_resolution_check = df.shape[0]
    if number_of_records_before_spatial_resolution_check > number_of_records_after_spatial_resolution_check:
        warnings.warn("Not all data within the specified range is 45 km resolution data.  "
                      "All non-45 km data was discarded.", category=Warning)

    # We will have one subplot for every year of observations
    n_rows = (year_range[1] + 1) - year_range[0]

    # Set up the plot
    fig, ax = plt.subplots(figsize=(8, 9), sharex='col', dpi=300, nrows=n_rows, ncols=1)
    plt.subplots_adjust(hspace=0.05)

    if time_units == "mlt":
        ax[n_rows - 1].set_xlabel('Time, MLT')
    elif time_units == "ut":
        ax[n_rows - 1].set_xlabel('Time, UT')

    # Apply common subplot formatting
    for row in reversed(range(n_rows)):
        ax[row].set_ylim([0, 13])
        ax[row].set_xlim(hour_range)
        ax[row].tick_params(axis='y', which='major', labelleft=False, direction='in')
        ax[row].tick_params(axis='x', which='major', direction='in')
        ax[row].xaxis.set_major_locator(MultipleLocator(2))  # Every 2 hours
        ax[row].yaxis.set_major_locator(MultipleLocator(6.5))  # Half way through the year
        ax[row].grid(b=True, which='major', axis='x', linestyle='--', linewidth=0.5, zorder=4)
        ax[row].grid(b=True, which='major', axis='y', linestyle='--', linewidth=0.5, zorder=4)
        ax[row].set_ylabel(year_range[1] - row)

    # Loop through the years, adding in the data
    n_bins_y = 24  # half month bins
    n_bins_x = (hour_range[1] - hour_range[0]) * 2  # half hour bins
    contour_range = [ax[0].get_xlim(), ax[0].get_ylim()]  # All plots are the same size
    for row in reversed(range(n_rows)):
        year = year_range[1] - row

        # Build a restricted dataframe with the year's data
        start_datetime, start_epoch = build_date_epoch(year=year, month=1, day=1, hour=0)
        end_datetime, end_epoch = build_date_epoch(year=year, month=12, day=31, hour=24)
        df_yy = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)].copy()
        df_yy.reset_index(drop=True, inplace=True)

        if time_units == "mlt":
            print("Computing MLTs for " + str(year) + " data...")

            # To compute mlt we need longitudes..
            # we will use the middle of the year and assume magnetic longitudes don't change much over the year
            mid_datetime = start_datetime + (end_datetime - start_datetime) / 2
            beam_corners_aacgm_lats, beam_corners_aacgm_lons = \
                radar_fov(stid=hdw_info.stid, coords='aacgm', date=mid_datetime)

            aacgm_lons = []
            dates = []
            for i in range(len(df_yy)):
                date = datetime.datetime(year, df_yy['month'][i], df_yy['day'][i],
                                         df_yy['hour'][i], df_yy['minute'][i], int(df_yy['second'][i]))
                dates.append(date)

                # TODO: Figure out if you need to get beam_corners_aacgm_lons for every time
                # beam_corners_aacgm_lats, beam_corners_aacgm_lons = \
                #     radar_fov(stid=hdw_info.stid, coords='aacgm', date=date)

                gate_corner = df_yy['slist'][i]
                beam_corner = df_yy['bmnum'][i]

                # estimate the cell with the centroid
                cent_lon, cent_lat = centroid([(beam_corners_aacgm_lons[gate_corner, beam_corner],
                                                beam_corners_aacgm_lats[gate_corner, beam_corner]),
                                               (beam_corners_aacgm_lons[gate_corner + 1, beam_corner],
                                                beam_corners_aacgm_lats[gate_corner + 1, beam_corner]),
                                               (beam_corners_aacgm_lons[gate_corner, beam_corner + 1],
                                                beam_corners_aacgm_lats[gate_corner, beam_corner + 1]),
                                               (beam_corners_aacgm_lons[gate_corner + 1, beam_corner + 1],
                                                beam_corners_aacgm_lats[gate_corner + 1, beam_corner + 1])])
                aacgm_lons.append(cent_lon)

            df_yy['mlt'] = aacgmv2.convert_mlt(arr=aacgm_lons, dtime=dates, m2a=False)
            df_yy = df_yy.loc[(df_yy['mlt'] >= hour_range[0]) & (df_yy['mlt'] <= hour_range[1])]
            df_yy['xdata'] = df_yy['mlt']
        else:
            # Use UT
            df_yy = df_yy.loc[(df_yy['hour'] >= hour_range[0]) & (df_yy['hour'] <= hour_range[1])]
            df_yy['xdata'] = df_yy['hour'] + df_yy['minute'] / 60 + df_yy['second'] / 3600

        df_yy.reset_index(drop=True, inplace=True)

        # Compute decimal datetime to plot along y
        df_yy['ydata'] = (df_yy['month'] - 1) + (df_yy['day'] - 1) / 31 + df_yy['xdata'] / 730

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
        cont = ax[row].contourf(bin_xcenters, bin_ycenters, binned_counts.transpose(), 5)
        fig.colorbar(cont, ax=ax[row])

    return fig


if __name__ == '__main__':
    """ Testing """

    local_testing = False

    station = "rkn"
    fig = occ_year_vs_ut(station=station, time_units='mlt', year_range=(2012, 2013),
                         gate_range=(20, 30), beam_range=(7, 7),
                         parameter='v', local_testing=local_testing)

    if local_testing:
        plt.show()
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        out_file = out_dir + "/occ_year_vs_time_" + station
        print("Saving plot as " + out_file)
        fig.savefig(out_file + ".jpg", format='jpg', dpi=300)
