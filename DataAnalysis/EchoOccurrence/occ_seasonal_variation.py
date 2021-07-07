import math
import os
import pathlib

import pandas as pd
import pydarn

import numpy as np
import matplotlib.ticker as mticker

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import SuperDARNRadars

from lib.boxcar_smooth import boxcar_smooth
from lib.plot_solar_flux_data import plot_solar_flux_data
from lib.plot_sunspot_data import plot_sunspot_data
from lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def occ_seasonal_variation(station, year_range=None, day_range=None, hour_range=None,
                           gate_range=None, beam_range=None, freq_range=None,
                           time_sector=None, time_units='lt', local_testing=False):
    """

    Produces an occurrence rate versus year plot that is meant to showcase the seasonal variation in echo occurrence.

    Also creates subplots for:
     F10.7 Flux (https://www.spaceweather.gc.ca/forecast-prevision/solar-solaire/solarflux/sx-3-en.php)
     Sunspot data (http://sidc.oma.be/silso/home)

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only considers 45 km data
        - To check which fitACF program is being used, refer to the data readers in lib.data_getters

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year_range: (int, int):
            Inclusive. The year range to consider. If omitted (or None), then all years will be considered.
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
    :param time_sector: str (optional):
            Time sectors are used to automatically consider predefined hour ranges that are of common interest
            "Day" (12 +/- 3 h), "Dusk" (18 +/- 3 h), "Night" (24 +/- 3 h), or "Dawn" (6 +/- 3 h)
            Each time sector contains 6 hours of observations, time sectors do not overlap
            If a time_sector is not None, then hour_range is ignored
    :param time_units: str:
            The provided hour_range is assumed to be in these time units.  Default is 'lt'.
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
    year_range = check_year_range(year_range)
    freq_range = check_freq_range(freq_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    night = False
    if time_sector is not None:
        # Then, as advertised, the time_sector will override the provided hour range
        time_sector = time_sector.lower()
        if time_sector == "day":  # From 9 AM to 3 PM
            hour_range = (9, 15)
        elif time_sector == "dusk":  # From 3 PM to 9 PM
            hour_range = (16, 21)
        elif time_sector == "night":  # From 9 PM to 3 AM (21 - 2)
            night = True
        elif time_sector == "dawn":  # From 3 AM to 9 AM
            hour_range = (3, 9)
        else:
            raise Exception("Time sector " + str(time_sector) + " not recognized.")
    else:
        # Consider the provided hour range
        hour_range = check_hour_range(hour_range)

    if year_range[1] == year_range[0]:
        year_string = str(year_range[0])
    else:
        year_string = str(year_range[0]) + " to " + str(year_range[1])
    beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    gate_string = "Gates " + str(gate_range[0]) + "-" + str(gate_range[1])
    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"
    if time_sector is not None:
        hour_string = "Time Sector: " + time_sector.capitalize() + " (" + time_units.upper() + ")"
    else:
        hour_string = "Hours " + str(hour_range[0]) + "-" + str(hour_range[1]) + " " + time_units.upper()

    print("Retrieving data...")
    df = get_data_handler(station, year_range=year_range, month_range=None, day_range=day_range,
                          gate_range=gate_range, beam_range=beam_range,
                          freq_range=freq_range, occ_data=True, local_testing=local_testing)
    df = only_keep_45km_res_data(df)

    # Add decimal hour to df in whatever units were requested
    # Use the middle of the mid year as magnetic field estimate
    mid_year = int(year_range[0] + (year_range[1] - year_range[0]) / 2)
    date_time_est, _ = build_datetime_epoch(mid_year, 6, 15, 0)
    df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id, date_time_est=date_time_est)

    # Restrict for valid hour ranges
    if night:
        # From 9 PM to 3 AM (21 - 3)
        df = df.loc[(df[time_units] >= 21) | (df[time_units] <= 3)]
    else:
        # We can restrict based on hour range normally
        df = df.loc[(df[time_units] >= hour_range[0]) & (df[time_units] <= hour_range[1])]

    df.reset_index(drop=True, inplace=True)

    print("Preparing the figure...")
    y_lim = [0, 0.8]
    y_axis_major_labels = [0.0, 0.2, 0.4, 0.6, 0.8]
    x_lim = [year_range[0], year_range[1] + 1]  # Make sure we go the end of the final year

    fig, axes = plt.subplots(figsize=[8, 12], dpi=300, constrained_layout=True, nrows=3, ncols=1)
    fig.suptitle(year_string + " at " + station.upper() + "; " + gate_string + "; " + beam_string +
                 "\n" + freq_string + "; " + hour_string +
                 "\nProduced by " + str(os.path.basename(__file__)), fontsize=18)

    occ_ax = axes[0]

    # Format the first plot for echo occurrence data
    occ_ax.set_xlim(x_lim)
    occ_ax.xaxis.set_major_locator(MultipleLocator(1))
    occ_ax.xaxis.set_minor_locator(MultipleLocator(0.25))

    occ_ax.set_ylim(y_lim)
    occ_ax.yaxis.set_major_locator(mticker.FixedLocator(y_axis_major_labels))
    occ_ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    occ_ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5, color='black')
    occ_ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2, color='black')

    occ_ax.set_xlabel("Year", fontsize=14)
    occ_ax.set_ylabel("Echo Occurrence Rate", fontsize=14)

    print("Computing UT decimal years/hours...")
    minutes_in_an_hour = 60
    seconds_in_an_hour = 3600
    months_in_a_year = 12
    days_in_a_year = 365
    hours_in_a_year = 8760

    decimal_hours, decimal_years = [], []
    for i in range(len(df)):
        datetime_obj_here = df['datetime'].iat[i]

        decimal_hour_here = datetime_obj_here.hour + datetime_obj_here.minute / minutes_in_an_hour + \
                            datetime_obj_here.second / seconds_in_an_hour

        decimal_year_here = datetime_obj_here.year + (datetime_obj_here.month - 1) / months_in_a_year + \
                            (datetime_obj_here.day - 1) / days_in_a_year + decimal_hour_here / hours_in_a_year

        decimal_hours.append(decimal_hour_here)
        decimal_years.append(decimal_year_here)

    df['decimal_hour'] = np.asarray(decimal_hours)
    df['decimal_year'] = np.asarray(decimal_years)

    # Compute year_edges
    bins_per_year = 180
    n_bins_x = ((year_range[1] + 1) - year_range[0]) * bins_per_year
    year_slice_edges = np.linspace(year_range[0], (year_range[1] + 1), num=(n_bins_x + 1))
    delta_year_slice = year_slice_edges[1] - year_slice_edges[0]

    # Compute bin centers
    bin_xwidth = delta_year_slice
    bin_xcenters = year_slice_edges[1:] - bin_xwidth / 2

    # Loop through the time bins and compute occurrence rate
    occurrence_data_is = np.empty(shape=n_bins_x)
    occurrence_data_gs = np.empty(shape=n_bins_x)
    for year_idx, year_slice_start in enumerate(year_slice_edges):
        if year_slice_start == year_slice_edges[-1]:
            continue  # The last edge is not a starting edge

        year_slice_end = year_slice_start + delta_year_slice

        df_yy = df[(df['decimal_year'] >= year_slice_start) & (df['decimal_year'] <= year_slice_end)]

        try:
            occurrence_data_is[year_idx] = sum(df_yy['good_iono_echo']) / len(df_yy)
            occurrence_data_gs[year_idx] = sum(df_yy['good_grndscat_echo']) / len(df_yy)
        except ZeroDivisionError:
            # There are no points in this interval
            occurrence_data_is[year_idx] = math.nan
            occurrence_data_gs[year_idx] = math.nan
        except BaseException as e:
            print("Something bad happened while computing binned occurrence rates..")
            print("Year index: " + str(year_idx))
            print("Year slice start: " + str(year_slice_start))
            raise e

    if local_testing:
        # Plot as a point so we can see it even while testing with little data
        occ_ax.plot(bin_xcenters, occurrence_data_is, 'bo', label='IS')
        occ_ax.plot(bin_xcenters, occurrence_data_gs, 'ro', label='GS')

    else:
        # Plot as a faint line
        occ_ax.plot(bin_xcenters, occurrence_data_is, linewidth=0.25, color="blue", linestyle='-', label='IS (Raw)')
        occ_ax.plot(bin_xcenters, occurrence_data_gs, linewidth=0.25, color="red", linestyle='-', label='GS (Raw)')

        # Compute and plot smoothed data
        # Before we smooth, we have to remove all nan values - otherwise we get discontinuities in the smoothed data
        is_mini_df = pd.DataFrame({'bin_xcenters': bin_xcenters, 'occurrence_data': occurrence_data_is})
        gs_mini_df = pd.DataFrame({'bin_xcenters': bin_xcenters, 'occurrence_data': occurrence_data_gs})

        is_mini_df = is_mini_df.loc[is_mini_df['occurrence_data'].notna()]
        gs_mini_df = gs_mini_df.loc[gs_mini_df['occurrence_data'].notna()]

        # 30 bins is about a 60 day boxcar filter
        is_mini_df['occurrence_data'] = boxcar_smooth(is_mini_df['occurrence_data'], window_size=30)
        gs_mini_df['occurrence_data'] = boxcar_smooth(gs_mini_df['occurrence_data'], window_size=30)

        occ_ax.plot(is_mini_df['bin_xcenters'], is_mini_df['occurrence_data'], linewidth=1, color="blue", linestyle='-',
                    label='IS (Smoothed)')
        occ_ax.plot(gs_mini_df['bin_xcenters'], gs_mini_df['occurrence_data'], linewidth=1, color="red", linestyle='-',
                    label='GS (Smoothed)')

    occ_ax.legend(loc='upper right')

    plot_solar_flux_data(ax=axes[1], year_range=(year_range[0], year_range[1] + 1), smoothing_window_size_in_days=30)
    plot_sunspot_data(ax=axes[2], year_range=(year_range[0], year_range[1] + 1), smoothing_window_size_in_days=30)

    return df, fig


if __name__ == '__main__':
    """ Testing """

    local_testing = False

    if local_testing:
        station = "rkn"

        # # Note: year, month, and day don't matter for local testing
        # df, fig = occ_seasonal_variation(station=station, year_range=(2011, 2012), day_range=(12, 12),
        #                                  gate_range=(10, 30), beam_range=(6, 8), freq_range=(11, 13),
        #                                  time_sector="night", time_units='lt', local_testing=local_testing)

        # Note: year, month, and day don't matter for local testing
        df, fig = occ_seasonal_variation(station=station, year_range=(2011, 2012), day_range=(12, 12),
                                         gate_range=(10, 30), beam_range=(6, 8), freq_range=(11, 13),
                                         local_testing=local_testing)

        plt.show()


    else:
        station = "dcn"
        year_range = (2019, 2021)
        freq_range = (8, 10)

        time_units = "lt"

        datetime_now = datetime.datetime.now()
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        _, fig = occ_seasonal_variation(station=station, year_range=year_range, day_range=None,
                                        gate_range=(10, 30), beam_range=(6, 8), freq_range=freq_range,
                                        time_units=time_units,
                                        local_testing=local_testing)

        out_fig = out_dir + "/occ_seasonalVariation_" + station + \
                  "_" + str(year_range[0]) + "-" + str(year_range[1]) + \
                  "_" + str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz_" + "_" + time_units

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)


        # time_sectors = ["day", "dusk", "night", "dawn"]
        # for time_sector in time_sectors:
        #     _, fig = occ_seasonal_variation(station=station, year_range=year_range, day_range=None,
        #                                     gate_range=(10, 30), beam_range=(6, 8), freq_range=freq_range,
        #                                     time_sector=time_sector, time_units=time_units,
        #                                     local_testing=local_testing)
        #
        #     out_fig = out_dir + "/occ_seasonalVariation_" + station + \
        #               "_" + str(year_range[0]) + "-" + str(year_range[1]) + \
        #               "_" + str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz_" + time_sector + "_" + time_units
        #
        #     print("Saving plot as " + out_fig)
        #     fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
