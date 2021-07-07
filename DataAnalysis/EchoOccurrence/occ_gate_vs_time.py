import math
import os
import pathlib

import pandas as pd
import pydarn
import calendar

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import SuperDARNRadars

from lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def occ_gate_vs_time(station, year, month_range=None, day_range=None, hour_range=None,
                     gate_range=None, beam_range=None, freq_range=None,
                     season=None, time_units='mlt', plot_type='contour',
                     local_testing=False):
    """

    Produce plots with gate on the y-axis and time along the x-axis.
    There are 2 subplots, one for ionospheric scatter and one for ground scatter
    Plots echo occurrence rates

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
    :param month_range: (<int>, <int>) (optional):
            Inclusive. The month range to consider.  If omitted (or None), then all months will be considered.
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
    :param season: str (optional):
            "Summer", "Spring", "Autumn", or "Winter"
            If a season is not None, then month_range is ignored
    :param plot_type: str (optional): default is 'contour'.
            The type of plot, either 'contour' or 'pixel'.
    :param time_units: str:
            The time units to plot along x.  Default is 'mlt'.
                'ut' for universal time
                'mlt' for magnetic local time
                'lt' for local time (based on longitude)
                'lst' for local standard time (based on time zones).
    :param local_testing: bool (optional): default is False.
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :return: pandas.DataFrame, matplotlib.pyplot.figure
            The dataframe used and the figure created.
            The figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    time_units = check_time_units(time_units)
    year = check_year(year)
    hour_range = check_hour_range(hour_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid
    hemisphere = this_radars_info.hemisphere

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    winter = False
    if season is not None:
        season = season.lower()
        # Then, as advertised, the season will override the provided month range
        if (hemisphere.value == 1 and season == "spring") or (hemisphere.value == -1 and season == "autumn"):
            month_range = (3, 5)
            mid_month = 4
        elif (hemisphere.value == 1 and season == "summer") or (hemisphere.value == -1 and season == "winter"):
            month_range = (6, 8)
            mid_month = 7
        elif (hemisphere.value == 1 and season == "autumn") or (hemisphere.value == -1 and season == "spring"):
            month_range = (9, 11)
            mid_month = 10
        elif (hemisphere.value == 1 and season == "winter") or (hemisphere.value == -1 and season == "summer"):
            winter = True
            mid_month = 1
        else:
            raise Exception("Season " + str(season) + " not recognized.")
    else:
        # Consider the provided month range
        month_range = check_month_range(month_range)
        mid_month = int(month_range[0] + (month_range[1] - month_range[0]) / 2)

    print("Retrieving data...")
    if winter:
        # We need Dec, Jan, Feb
        df1 = get_data_handler(station, year_range=(year, year), month_range=(1, 2), day_range=day_range,
                               gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                               occ_data=True, local_testing=local_testing)
        # Get Dec data separately in df2
        df2 = get_data_handler(station, year_range=(year, year), month_range=(12, 22), day_range=day_range,
                               gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                               occ_data=True, local_testing=local_testing)
        df = pd.concat([df1, df2])

    else:
        # We can go ahead and get data normally
        print("Using month range: " + str(month_range))
        df = get_data_handler(station, year_range=(year, year), month_range=month_range, day_range=day_range,
                              gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                              occ_data=True, local_testing=local_testing)
    df = only_keep_45km_res_data(df)

    # Add decimal hour to df in whatever units were requested
    date_time_est, _ = build_datetime_epoch(year=year, month=mid_month, day=15, hour=0)
    df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id, date_time_est=date_time_est)

    df = df.loc[(df[time_units] >= hour_range[0]) & (df[time_units] <= hour_range[1])]
    df.reset_index(drop=True, inplace=True)

    print("Preparing the plot...")
    if season is not None:
        month_string = season.capitalize()
    else:
        # We want to state the month range
        if month_range[1] == month_range[0]:
            month_string = calendar.month_name[month_range[0]]
        else:
            month_string = calendar.month_name[month_range[0]] + " to " + calendar.month_name[month_range[1]]
    beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"

    # Setup the plot
    fig, ax = plt.subplots(figsize=[10, 8], dpi=300, nrows=2, ncols=1, constrained_layout=True)
    fig.suptitle(month_string + " " + str(year) + " at " + station.upper() + "; " + beam_string + "; " + freq_string +
                 "\nProduced by " + str(os.path.basename(__file__)), fontsize=18)

    ax[0].set_title("Ionospheric Scatter", fontsize=16)
    ax[1].set_title("Ground Scatter", fontsize=16)

    # Apply common subplot formatting
    for i in range(ax.size):
        ax[i].set_ylim(gate_range)
        ax[i].set_xlim(hour_range)
        ax[i].xaxis.set_major_locator(MultipleLocator(4))
        ax[i].yaxis.set_major_locator(MultipleLocator(10))
        ax[i].tick_params(axis='both', which='major', direction='in', color='white', labelsize=14)
        ax[i].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5, zorder=4, color='white')
        ax[i].set_ylabel("Range Gate", fontsize=16)
        ax[i].set_xlabel("Time, " + time_units.upper(), fontsize=16)

    # Compute hour_edges
    bins_per_hour = 4  # quarter hour bins
    n_bins_x = int((hour_range[1] - hour_range[0]) * bins_per_hour)
    delta_hour = (hour_range[1] - hour_range[0]) / n_bins_x
    hour_edges = np.linspace(hour_range[0], hour_range[1], num=(n_bins_x + 1))

    # Compute gate_edges
    bins_per_gate = 1  # Single gate bins
    n_bins_y = int(((gate_range[1] + 1) - gate_range[0]) * bins_per_gate)
    gate_edges = np.linspace(gate_range[0], gate_range[1] + 1, num=(n_bins_y + 1), dtype=int)

    print("Computing binned occ rates...")
    contour_data_is = np.empty(shape=(n_bins_x, n_bins_y))
    contour_data_is[:] = math.nan

    contour_data_gs = np.empty(shape=(n_bins_x, n_bins_y))
    contour_data_gs[:] = math.nan

    for hour_idx, hour_start in enumerate(hour_edges):
        if hour_start == hour_edges[-1]:
            continue  # The last edge is not a starting hour
        hour_end = hour_start + delta_hour
        df_hh = df[(df[time_units] >= hour_start) & (df[time_units] <= hour_end)]

        for gate in gate_edges:
            if gate == gate_edges[-1]:
                continue  # The last edge is not a starting gate
            df_hh_gg = df_hh[df_hh['slist'] == gate]

            try:
                contour_data_is[hour_idx][gate] = sum(df_hh_gg['good_iono_echo']) / len(df_hh_gg)
                contour_data_gs[hour_idx][gate] = sum(df_hh_gg['good_grndscat_echo']) / len(df_hh_gg)
            except ZeroDivisionError:
                # There are no points in this interval
                contour_data_is[hour_idx, gate] = math.nan
                contour_data_gs[hour_idx, gate] = math.nan
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

    if plot_type == "contour":
        plot0 = ax[0].contourf(bin_xcenters, bin_ycenters, contour_data_is.transpose(), cmap='jet', levels=levels)
        plot1 = ax[1].contourf(bin_xcenters, bin_ycenters, contour_data_gs.transpose(), cmap='jet', levels=levels)

    elif plot_type == "pixel":
        plot0 = ax[0].imshow(np.flip(contour_data_is.transpose(), axis=0), aspect='auto', cmap="jet",
                             extent=(hour_range[0], hour_range[1], gate_range[0], gate_range[1] + 1), vmin=0, vmax=1)
        plot1 = ax[1].imshow(np.flip(contour_data_gs.transpose(), axis=0), aspect='auto', cmap="jet",
                             extent=(hour_range[0], hour_range[1], gate_range[0], gate_range[1] + 1), vmin=0, vmax=1)
    else:
        raise Exception("plot_type not recognized")

    cbar0 = fig.colorbar(plot0, ax=ax[0], shrink=0.75, format='%.2f')
    cbar0.ax.tick_params(labelsize=14)

    cbar1 = fig.colorbar(plot1, ax=ax[1], shrink=0.75, format='%.2f')
    cbar1.ax.tick_params(labelsize=14)

    print("Returning the is and gs figures...")
    return df, fig


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    if local_testing:
        station = "rkn"

        # Note: year, month, and day don't matter for local testing
        df, fig = occ_gate_vs_time(station=station, year=2011, month_range=None, day_range=(12, 12), hour_range=None,
                                   gate_range=(0, 74), beam_range=None, freq_range=(11, 13),
                                   season="Spring", time_units='ut', plot_type='contour',
                                   local_testing=local_testing)

        plt.show()


    else:
        station = "dcn"
        freq_range = (8, 10)
        seasons = ["Spring", "Summer", "Autumn", "Winter"]

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        for year in range(2019, 2022, 1):
            for season in seasons:

                _, fig = occ_gate_vs_time(station=station, year=year, day_range=None,
                                          gate_range=(0, 74), beam_range=(6, 8), freq_range=freq_range,
                                          season=season, time_units='ut', plot_type='pixel',
                                          local_testing=local_testing)

                out_fig = out_dir + "/occ_gateVtime_" + station + "-" + str(year) + "-" + season + "_" + \
                          str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz"

                print("Saving plot as " + out_fig)
                fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
