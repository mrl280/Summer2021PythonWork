import math
import pathlib
import pydarn
import calendar

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import radar_fov, SuperDARNRadars
import matplotlib.ticker as mticker

from lib.add_mlt_to_df import add_mlt_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def occ_diurnal_variation(station, year, day_range=None, hour_range=None,
                          gate_range=None, beam_range=None, freq_range=None,
                          time_units='mlt', local_testing=False):
    """

    Produce a series of occurrence plots that showcase diurnal variation.

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only considers 45 km data
        - To check which fitACF program is being used, refer to the data readers in lib.data_getters

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year: int:
            The year to consider
    :param day_range: (<int>, <int>) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param hour_range: (<int>, <int>) (optional):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not quite inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param gate_range: (<int>, <int>) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param freq_range: (<float>, <float>) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.
    :param time_units: str: 'ut' for universal time or 'mlt' for magnetic local time:
            The time units to plot along x.
    :param local_testing: bool (optional): default is False.
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.

    :return: pandas.DataFrame, matplotlib.pyplot.figure
            The dataframe used and the figure created.
            The figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    month_offset = 3

    time_units = check_time_units(time_units)
    year = check_year(year)
    hour_range = check_hour_range(hour_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid
    hemisphere = this_radars_info.hemisphere

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    gate_string = "Gates " + str(gate_range[0]) + "-" + str(gate_range[1])
    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"

    print("Preparing the figure...")
    fig = plt.figure(figsize=[10, 12], constrained_layout=True, dpi=300)
    gs = fig.add_gridspec(8, 6)

    # Build all of the spring time axes
    mar_axis = fig.add_subplot(gs[0, 0])
    apr_axis = fig.add_subplot(gs[0, 1])
    may_axis = fig.add_subplot(gs[0, 2])
    spring_axis = fig.add_subplot(gs[1:4, 0:3])  # Splices don't include last index
    spring_month_range = (3, 5)

    # Build all of the summer time axes
    jun_axis = fig.add_subplot(gs[4, 0])
    jul_axis = fig.add_subplot(gs[4, 1])
    aug_axis = fig.add_subplot(gs[4, 2])
    summer_axis = fig.add_subplot(gs[5:8, 0:3])  # Splices don't include last index
    summer_month_range = (6, 8)

    # Build all of the fall time axes
    sep_axis = fig.add_subplot(gs[0, 3])
    oct_axis = fig.add_subplot(gs[0, 4])
    nov_axis = fig.add_subplot(gs[0, 5])
    fall_axis = fig.add_subplot(gs[1:4, 3:6])  # Splices don't include last index
    fall_month_range = (9, 11)

    # Build all of the winter time axes
    dec_axis = fig.add_subplot(gs[4, 3])
    jan_axis = fig.add_subplot(gs[4, 4])
    feb_axis = fig.add_subplot(gs[4, 5])
    winter_axis = fig.add_subplot(gs[5:8, 3:6])  # Splices don't include last index

    # Group axes for easy iteration
    month_axes = [mar_axis, apr_axis, may_axis,
                  jun_axis, jul_axis, aug_axis,
                  sep_axis, oct_axis, nov_axis,
                  dec_axis, jan_axis, feb_axis]
    season_axes = [spring_axis, summer_axis, fall_axis, winter_axis]

    y_lim = [0, 0.8]
    y_axis_major_labels = [0.0, 0.2, 0.4, 0.6, 0.8]

    # Apply common subplot formatting
    for ax in month_axes + season_axes:
        ax.set_xlim(hour_range)
        ax.xaxis.set_major_locator(MultipleLocator(6))
        ax.xaxis.set_minor_locator(MultipleLocator(2))

        ax.set_ylim(y_lim)
        ax.yaxis.set_major_locator(mticker.FixedLocator(y_axis_major_labels))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))

        ax.tick_params(axis='both', which='both', direction='in', color='black', left=True, right=True)
        ax.yaxis.set_ticklabels([])

        if ax in season_axes:
            ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5, color='black')
            ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2, color='black')
        else:
            ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.2, color='black')

    # Print out the months of the axis
    for i, ax in enumerate(month_axes):
        month_num = i + month_offset
        if month_num > 12:
            month_num = month_num - 12
        ax.text(ax.get_xlim()[0] + 0.1 * (ax.get_xlim()[1] - ax.get_xlim()[0]), 0.9 * ax.get_ylim()[1],
                calendar.month_name[month_num], ha='left', va='top')

    if hemisphere.value == 1:
        spring_string = "Spring"
        summer_string = "Summer"
        fall_string = "Autumn"
        winter_string = "Winter"
    elif hemisphere.value == -1:
        # The seasons are reversed in the southern hemisphere
        spring_string = "Autumn"
        summer_string = "Winter"
        fall_string = "Spring"
        winter_string = "Summer"
    else:
        raise Exception("Error: Season not recognized.")

    season_txt_fontsize = 16
    spring_axis.text(spring_axis.get_xlim()[0] + 0.05 * (spring_axis.get_xlim()[1] - spring_axis.get_xlim()[0]),
                     0.9 * spring_axis.get_ylim()[1], spring_string, ha='left', va='center', fontsize=season_txt_fontsize)
    summer_axis.text(summer_axis.get_xlim()[0] + 0.05 * (summer_axis.get_xlim()[1] - summer_axis.get_xlim()[0]),
                     0.9 * summer_axis.get_ylim()[1], summer_string, ha='left', va='center', fontsize=season_txt_fontsize)
    fall_axis.text(fall_axis.get_xlim()[0] + 0.05 * (fall_axis.get_xlim()[1] - fall_axis.get_xlim()[0]),
                   0.9 * fall_axis.get_ylim()[1], fall_string, ha='left', va='center', fontsize=season_txt_fontsize)
    winter_axis.text(winter_axis.get_xlim()[0] + 0.05 * (winter_axis.get_xlim()[1] - winter_axis.get_xlim()[0]),
                     0.9 * winter_axis.get_ylim()[1], winter_string, ha='left', va='center', fontsize=season_txt_fontsize)

    fig.suptitle(str(year) + " at " + station.upper() + "; " + gate_string +
                 "\n" + beam_string + "; " + freq_string, fontsize=18)

    # Put ticks and labels around the edge
    mar_axis.yaxis.set_ticklabels(y_axis_major_labels, zorder=3)
    jun_axis.yaxis.set_ticklabels(y_axis_major_labels, zorder=3)
    spring_axis.yaxis.set_ticklabels(y_axis_major_labels, zorder=3)
    summer_axis.yaxis.set_ticklabels(y_axis_major_labels, zorder=3)

    summer_axis.set_xlabel("Time, " + time_units.upper(), fontsize=14)
    summer_axis.set_ylabel("Echo Occurrence Rate", fontsize=14)

    print("Retrieving data...")
    df = get_data_handler(station, year_range=(year, year), month_range=None, day_range=day_range,
                          gate_range=gate_range, beam_range=beam_range, freq_range=freq_range, occ_data=True,
                          local_testing=local_testing)
    df = only_keep_45km_res_data(df)

    # Get our raw x-data
    if time_units == "mlt":
        print("Computing MLTs...")

        # To compute mlt we need longitudes.. use the middle of the year as magnetic field estimate
        date_time_est, _ = build_datetime_epoch(year, 6, 15, 0)
        cell_corners_aacgm_lats, cell_corners_aacgm_lons = radar_fov(stid=radar_id, coords='aacgm', date=date_time_est)

        df = add_mlt_to_df(cell_corners_aacgm_lons=cell_corners_aacgm_lons,
                           cell_corners_aacgm_lats=cell_corners_aacgm_lats, df=df)

        df['xdata'] = df['mlt']

    else:
        print("Computing UTs...")

        ut_time = []
        for i in range(len(df)):
            ut_time_here = df['datetime'].iat[i].hour + df['datetime'].iat[i].minute / 60 + \
                           df['datetime'].iat[i].second / 3600

            ut_time.append(ut_time_here)

        df['xdata'] = np.asarray(ut_time)

    df = df.loc[(df['xdata'] >= hour_range[0]) & (df['xdata'] <= hour_range[1])]
    df.reset_index(drop=True, inplace=True)

    # Compute month, we will need it to filter the data
    month = []
    for i in range(len(df)):
        month.append(df['datetime'].iat[i].month)
    df['month'] = month

    # Compute hour_edges
    bins_per_hour = 4
    n_bins_x = (hour_range[1] - hour_range[0]) * bins_per_hour  # quarter hour bins
    hour_edges = np.linspace(hour_range[0], hour_range[1], num=(n_bins_x + 1))
    delta_hour = hour_edges[1] - hour_edges[0]

    # Compute bin centers
    bin_xwidth = (hour_edges[1] - hour_edges[0])
    bin_xcenters = hour_edges[1:] - bin_xwidth / 2

    add_month_data_to_plot(df=df, month_axes=month_axes, month_offset=month_offset,
                           hour_edges=hour_edges, bin_xcenters=bin_xcenters)

    # TODO: Add the seasonal plot additions to a separate function.
    """ Loop through all the seasonal plots """
    for i, ax in enumerate(season_axes):
        if ax == spring_axis:
            month_rest_range = spring_month_range
        elif ax == summer_axis:
            month_rest_range = summer_month_range
        elif ax == fall_axis:
            month_rest_range = fall_month_range
        elif ax == winter_axis:
            # We want Dec, Jan, and Feb
            # Addition and subtraction of 1 are because month filter is inclusive
            month_rest_range = (fall_month_range[1] + 1, spring_month_range[0] - 1)  # usually (12, 2)
        else:
            raise Exception("Season not recognized")

        if ax == winter_axis:
            df_season = df[(df['month'] >= month_rest_range[0]) | (df['month'] <= month_rest_range[1])]
        else:
            df_season = df[(df['month'] >= month_rest_range[0]) & (df['month'] <= month_rest_range[1])]

        print("Computing occurrence data for season number " + str(i))

        # Loop through the time bins and compute occurrence rate
        occurrence_data_is = np.empty(shape=n_bins_x)
        occurrence_data_gs = np.empty(shape=n_bins_x)
        for hour_idx, hour_start in enumerate(hour_edges):
            if hour_start == hour_edges[-1]:
                continue  # The last edge is not a starting hour

            hour_end = hour_start + delta_hour
            df_season_hh = df_season[(df_season['xdata'] >= hour_start) & (df_season['xdata'] <= hour_end)]

            try:
                occurrence_data_is[hour_idx] = sum(df_season_hh['good_iono_echo']) / len(df_season_hh)
                occurrence_data_gs[hour_idx] = sum(df_season_hh['good_grndscat_echo']) / len(df_season_hh)
            except ZeroDivisionError:
                # There are no points in this interval
                occurrence_data_is[hour_idx] = math.nan
                occurrence_data_gs[hour_idx] = math.nan
            except BaseException as e:
                print("Something bad happened while making the month plots..")
                print("Season: " + str(i))
                print("Hour index: " + str(hour_idx))
                raise e

        # Plot the data
        ax.plot(bin_xcenters, occurrence_data_is, color="blue", linestyle='-', label='IS')
        ax.plot(bin_xcenters, occurrence_data_gs, color="red", linestyle='-', label='GS')

    print("Returning the dataframe and figure...")
    return df, fig


def add_month_data_to_plot(df, month_axes, month_offset, hour_edges, bin_xcenters):
    """
    Fill in all of the monthly subplots.

    :param df: pandas.DataFrame:
    :param month_axes: numpy.array of matplotlib.axes: All of the individual month axes
    :param month_offset: int: January is not necessarily the first set of axis, so we need this offset
    :param hour_edges: numpy.array: The edges of all the hour (x-data) bins.
    :param bin_xcenters: numpy.array: The x-axis bin centers, this is the x data for plotting.
                        There is one value for every bin.
    """

    n_bins_x = len(bin_xcenters)
    delta_hour = hour_edges[1] - hour_edges[0]

    for i, ax in enumerate(month_axes):
        month_num = i + month_offset
        if month_num > 12:
            month_num = month_num - 12

        print("Computing occurrence data for " + calendar.month_name[month_num])

        df_mm = df[df['month'] == month_num]

        # Loop through the time bins and compute occurrence rate
        occurrence_data_is = np.empty(shape=n_bins_x)
        occurrence_data_gs = np.empty(shape=n_bins_x)
        for hour_idx, hour_start in enumerate(hour_edges):
            if hour_start == hour_edges[-1]:
                continue  # The last edge is not a starting hour

            hour_end = hour_start + delta_hour
            df_mm_hh = df_mm[(df_mm['xdata'] >= hour_start) & (df_mm['xdata'] <= hour_end)]

            try:
                occurrence_data_is[hour_idx] = sum(df_mm_hh['good_iono_echo']) / len(df_mm_hh)
                occurrence_data_gs[hour_idx] = sum(df_mm_hh['good_grndscat_echo']) / len(df_mm_hh)
            except ZeroDivisionError:
                # There are no points in this interval
                occurrence_data_is[hour_idx] = math.nan
                occurrence_data_gs[hour_idx] = math.nan
            except BaseException as e:
                print("Something bad happened while making the month plots..")
                print("Month: " + calendar.month_name[month_num])
                print("Hour index: " + str(hour_idx))
                raise e

        # Plot the data
        ax.plot(bin_xcenters, occurrence_data_is, color="blue", linestyle='-', linewidth=0.5, label='IS')
        ax.plot(bin_xcenters, occurrence_data_gs, color="red", linestyle='-', linewidth=0.5, label='GS')


if __name__ == '__main__':
    """ Testing """

    local_testing = False

    if local_testing:
        station = "rkn"

        # Note: year, month, and day don't matter for local testing
        df, fig = occ_diurnal_variation(station=station, year=2011, day_range=(12, 12),
                                        gate_range=(10, 30), beam_range=(6, 8), freq_range=(11, 13),
                                        time_units='ut', local_testing=local_testing)

        plt.show()


    else:
        station = "dcn"
        freq_range = (8, 10)

        datetime_now = datetime.datetime.now()
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        for year in range(2019, 2022, 1):

            df, fig = occ_diurnal_variation(station=station, year=year, day_range=None,
                                           gate_range=(10, 30), beam_range=(6, 8), freq_range=freq_range,
                                           time_units='ut', local_testing=local_testing)

            out_fig = out_dir + "/occ_diurnalVariation_" + station + "-" + str(year) + "_" + \
                      str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz"

            print("Saving plot as " + out_fig)
            fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)

            # print("Saving df as " + out_fig)
            # df.to_pickle(out_fig + ".pkl")
