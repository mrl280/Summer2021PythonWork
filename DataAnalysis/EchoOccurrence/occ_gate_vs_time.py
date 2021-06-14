import math
import pathlib
import pydarn
import calendar

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import radar_fov, SuperDARNRadars

from lib.add_mlt_to_df import add_mlt_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def occ_gate_vs_time(station, year, month, day_range=None, hour_range=None, time_units='mlt',
                     gate_range=None, beam_range=None, local_testing=False):
    """

    Produce a contour plot with gate on the y-axis and time along the x-axis.
    Plots echo occurrence rate

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only considers 45 km data
        - Does not distinguish frequency
        - This program uses fitACF 3.0 data.  To change this, modify the source code.

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year: int:
            The year to consider
    :param month: int:
            The month to consider
    :param day_range: (<int>, <int>) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
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
    :return: matplotlib.pyplot.figure
            The figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    time_units = check_time_units(time_units)
    year = check_year(year)
    month = check_month(month)
    hour_range = check_hour_range(hour_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    print("Retrieving data...")
    df = get_data_handler(station, year_range=(year, year), month_range=(month, month), hour_range=hour_range,
                          day_range=day_range, gate_range=gate_range, beam_range=beam_range, occ_data=True,
                          local_testing=local_testing)
    df = only_keep_45km_res_data(df)

    # Get our raw x-data
    if time_units == "mlt":
        print("Computing MLTs for " + str(year) + " data...")

        # To compute mlt we need longitudes.. use the middle of the month as magnetic field estimate
        date_time_est, _ = build_datetime_epoch(year, month, 15, 0)
        cell_corners_aacgm_lats, cell_corners_aacgm_lons = \
            radar_fov(stid=radar_id, coords='aacgm', date=date_time_est)

        df = add_mlt_to_df(cell_corners_aacgm_lons=cell_corners_aacgm_lons,
                           cell_corners_aacgm_lats=cell_corners_aacgm_lats, df=df)

        df['xdata'] = df['mlt']

    else:
        print("Computing UTs for " + str(year) + " data...")

        # I tried to speed this up by only recomputing when necessary
        ut_time = []
        datetime_save = df['datetime'].iat[0]
        ut_save = df['datetime'].iat[0].hour + df['datetime'].iat[0].minute / 60 + df['datetime'].iat[0].second / 3600
        for i in range(len(df)):
            if df['datetime'].iat[i] == datetime_save:
                ut_time.append(ut_save)
            else:
                datetime_save = df['datetime'].iat[i]
                ut_save = df['datetime'].iat[i].hour + df['datetime'].iat[i].minute / 60 + df['datetime'].iat[
                    i].second / 3600
                ut_time.append(ut_save)

        df['xdata'] = ut_time

    df = df.loc[(df['xdata'] >= hour_range[0]) & (df['xdata'] <= hour_range[1])]
    df.reset_index(drop=True, inplace=True)

    print("Preparing the plot...")
    fig, ax = plt.subplots(figsize=[8, 6], dpi=300)
    ax.set_ylim(gate_range)
    ax.set_xlim(hour_range)
    ax.xaxis.set_major_locator(MultipleLocator(4))
    ax.tick_params(axis='both', which='major', direction='in', color='white')
    ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5, zorder=4, color='white')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_ylabel("Range Gate", fontsize=16)
    ax.set_xlabel("Time, " + time_units.upper(), fontsize=16)
    plt.title(calendar.month_name[month] + " " + str(year) + " at " + station.upper() +
              "; Beams " + str(beam_range[0]) + "-" + str(beam_range[1]), fontsize=18)

    # Compute hour_edges
    bins_per_hour = 4
    n_bins_x = (hour_range[1] - hour_range[0]) * bins_per_hour  # quarter hour bins
    delta_hour = (hour_range[1] - hour_range[0]) / n_bins_x
    hour_edges = np.linspace(hour_range[0], hour_range[1], num=(n_bins_x + 1))

    # Compute gate_edges
    n_bins_y = ((gate_range[1] + 1) - gate_range[0])  # Single gate bins
    gate_edges = np.linspace(gate_range[0], gate_range[1] + 1, num=(n_bins_y + 1), dtype=int)

    print("Computing binned occ rates...")
    contour_data = np.empty(shape=(n_bins_x, n_bins_y))
    contour_data[:] = math.nan

    for hour_idx, hour_start in enumerate(hour_edges):
        if hour_start == hour_edges[-1]:
            continue  # The last edge is not a starting hour
        hour_end = hour_start + delta_hour
        df_hh = df[(df['xdata'] >= hour_start) & (df['xdata'] >= hour_end)]

        for gate in gate_edges:
            if gate == gate_edges[-1]:
                continue  # The last edge is not a starting gate
            df_hh_gg = df_hh[df_hh['slist'] == gate]

            try:
                contour_data[hour_idx][gate] = sum(df_hh_gg['good_echo']) / len(df_hh_gg)
            except ZeroDivisionError:
                # There are no point in this interval
                contour_data[hour_idx, gate] = math.nan
            except BaseException as e:
                print("Hour index: " + str(hour_idx))
                print("Gate: " + str(gate))
                raise e

    # Compute bin centers
    bin_xwidth = (hour_edges[1] - hour_edges[0])
    bin_ywidth = (gate_edges[1] - gate_edges[0])
    bin_xcenters = hour_edges[1:] - bin_xwidth / 2
    bin_ycenters = gate_edges[1:] - bin_ywidth / 2

    # Plot the data
    levels = 12
    levels = np.linspace(start=0, stop=1, num=(levels + 1))
    cont = ax.contourf(bin_xcenters, bin_ycenters, contour_data.transpose(),
                       cmap='jet', levels=levels)
    cbar = fig.colorbar(cont, ax=ax, shrink=0.75, format='%.2f')
    cbar.ax.tick_params(labelsize=14)

    print("Returning the figure...")
    return fig


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    if local_testing:
        station = "rkn"

        fig = occ_gate_vs_time(station=station, year=2011, month=9, day_range=None, time_units='ut',
                               gate_range=(0, 74), beam_range=None, local_testing=local_testing)

        plt.show()


    else:
        station = "inv"
        year = 2014
        month = 2

        fig = occ_gate_vs_time(station=station, year=year, month=month, day_range=(1, 3), time_units='ut',
                               gate_range=(0, 74), beam_range=(13, 15), local_testing=local_testing)

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        out_file = out_dir + "/occ_gate_vs_time_" + station + str(year) + str(month)
        print("Saving plot as " + out_file)
        fig.savefig(out_file + ".jpg", format='jpg', dpi=300)
