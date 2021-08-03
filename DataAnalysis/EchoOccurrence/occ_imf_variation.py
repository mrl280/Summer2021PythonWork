import calendar
import math
import os
import pathlib
import pydarn

import pandas as pd
import numpy as np
import matplotlib.ticker as mticker

from bisect import bisect
from matplotlib.ticker import MultipleLocator
from pydarn import SuperDARNRadars
from matplotlib import pyplot as plt

from lib.data_getters.get_IMF_data import get_IMF_data
from lib.boxcar_smooth import boxcar_smooth
from lib.plot_solar_flux_data import plot_solar_flux_data
from lib.plot_sunspot_data import plot_sunspot_data
from lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def occ_imf_variation(station, year=None, day_range=None, hour_range=None,
                      gate_range=None, beam_range=None, freq_range=None,
                      time_units='mlt', local_testing=False, even_odd_days=None):
    """

    TODO: Plot df['By_nT_GSM'] along x, df['Bz_nT_GSM'] along y, and then echo occurrence rate along z

    Plot the y-component of the interplanetary magnetic field (By) along x, and the z-component along y.
      Occurrence is then plotted as contours or pixels.  One plot for every month.

    Idea is to show that IMF drives echo occurrence rates

    Uses OMNI IMF data.  IMF data should be pickled with pickle_imf() before running this program.
     An error is raised if no pickled data is found.

    More on OMNI data here: https://omniweb.gsfc.nasa.gov/html/omni_min_data.html
    OMNI data files were created and downloaded from here: https://omniweb.gsfc.nasa.gov/form/omni_min.html

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - To check which fitACF program is being used, refer to the data readers in lib.data_getters

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year: int:
            The year to consider.
    :param day_range: (int, int) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param hour_range: (int, int) (optional):
            This parameter only effects the all hours plot.
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not quite inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param gate_range: (int, int) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (int, int) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param freq_range: (float, float) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.

    :param time_units: str:
            The time units to plot along x.  Default is 'mlt'.
                'ut' for universal time
                'mlt' for magnetic local time
                'lt' for local time (based on longitude)
                'lst' for local standard time (based on time zones)
                'ast' for apparent solar time (based on the apparent angular motion of the sun across the sky)
    :param even_odd_days: (optional; default is None)
            'even': only even days are read in
            'odd': only odd days are read in
            None: all days are read in
    :param local_testing: bool (optional; default is False)
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.

    :return fig: matplotlib.pyplot.figure
            The figure created, it can then be modified, added to, printed out, or saved in whichever file format is
            desired.
    """

    vmaxs = {"is": 0.6,  # The top of the colour bar
             "gs": 0.2}

    year = check_year(year)
    freq_range = check_freq_range(freq_range)
    day_range = check_day_range(day_range)
    hour_range = check_hour_range(hour_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    print("     Retrieving SuperDARN data...")
    # We will force build a new dataframe because the pickled one probably won't be for all beams/gates
    df = get_data_handler(station, year_range=(year, year), month_range=None, day_range=day_range,
                          gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                          occ_data=True, local_testing=local_testing, even_odd_days=even_odd_days)
    df = only_keep_45km_res_data(df)
    df.sort_values(by='datetime', ascending=True, inplace=True)  # Just to be safe, needed for bisect matching

    # Perform hour restriction if required
    if hour_range != (0, 24):
        # Add decimal hour to df in whatever units were requested
        # Use the middle of the year as magnetic field estimate
        date_time_est, _ = build_datetime_epoch(year=year, month=6, day=15, hour=0)
        df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id, date_time_est=date_time_est)
        df = df.loc[(df[time_units] >= hour_range[0]) & (df[time_units] <= hour_range[1])]
        df.reset_index(drop=True, inplace=True)

    print("     Retrieving IMF data...")
    imf_df = get_IMF_data(year_range=(year, year))

    print("     Preparing the figure...")
    fig = plt.figure(figsize=[10, 14], constrained_layout=True, dpi=300)
    month_axes, cbar_axes = add_axes(fig=fig)
    apply_subplot_formatting(axes=month_axes)

    beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    gate_string = "Gates " + str(gate_range[0]) + "-" + str(gate_range[1])
    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"
    if day_range == (1, 31):
        day_string = "All Days"
    else:
        day_string = "Days " + str(day_range[0]) + "-" + str(day_range[1])
    if hour_range == (0, 24):
        hour_string = "All Hours"
    else:
        hour_string = "Hours " + str(hour_range[0]) + "-" + str(hour_range[1]) + " " + time_units.upper()

    fig.suptitle(str(year) + " at " + station.upper() + "; " + day_string + "; " + hour_string +
                 "\n" + gate_string + "; " + beam_string + "; " + freq_string +
                 "\nProduced by " + str(os.path.basename(__file__)), fontsize=18)

    print("     Adding month to df...")
    month = []
    for i in range(len(df)):
        month.append(df['datetime'].iat[i].month)
    df['month'] = month

    for month_str, month_dict_item in month_axes.items():

        month = list(calendar.month_abbr).index(month_str[:3])  # Pull the month number from the month string

        if month != 11:  # TODO: Remove this condition
            continue

        # We need to assign an IMF value to every value in df
        # We will use the bisect method proposed here because it is apparently O(log n) on a sorted list:
        #  https://stackoverflow.com/questions/29700214/get-the-closest-datetime-from-a-list
        # I am not really sure why this works but it seems to
        print("     Assigning IMF values for " + month_str)

        # Month restrict both dataframes, this should speed up the date matching
        df_mm = df[df['month'] == month].copy()
        df_mm.reset_index(drop=True, inplace=True)
        imf_df_mm = imf_df[imf_df['month'] == month].copy()
        imf_df_mm.reset_index(drop=True, inplace=True)

        By_nT_GSM, Bz_nT_GSM = [], []
        print("The dataframe here is of length: " + str(len(df_mm)))
        for i in range(len(df)):
            index = bisect(imf_df_mm['datetime'], df_mm['datetime'].iat[i], hi=len(imf_df_mm['datetime']) - 1)

            # print("")
            # print("Looking for a matching date to " + str(df_mm['datetime'].iat[i]))
            # print("What we found was " + str(imf_df_mm['datetime'].iat[index]))

            By_nT_GSM.append(imf_df_mm['By_nT_GSM'].iat[index])
            Bz_nT_GSM.append(imf_df_mm['Bz_nT_GSM'].iat[index])

        df_mm['By_nT_GSM'] = By_nT_GSM
        df_mm['Bz_nT_GSM'] = Bz_nT_GSM

        df_mm = df_mm.loc[df_mm['By_nT_GSM'].notna() & df_mm['Bz_nT_GSM'].notna()]
        print("after removing nan IMF values, the dataframe here is of length: " + str(len(df_mm)))
        print(df_mm.head())

        # TODO: Plot the data as both contours and pixels
        for subplot_type, ax in month_dict_item.items():

            pass

    return fig


def add_axes(fig):
    """
    :param fig: The figure to draw the axes on

    :return month_axes, cbar_axes: dictionary of axes, dictionary of axes
            Month axes will be keyed by month
    """

    # Leave a blank column in the middle to separate the is and gs subplots
    gs = fig.add_gridspec(ncols=9, nrows=7)

    month_axes, cbar_axes = dict(), dict()

    month_axes["January"] = {"is": fig.add_subplot(gs[0, 0:2]),
                             "gs": fig.add_subplot(gs[0, 5:7])}
    month_axes["February"] = {"is": fig.add_subplot(gs[0, 2:4]),
                              "gs": fig.add_subplot(gs[0, 7:9])}

    month_axes["March"] = {"is": fig.add_subplot(gs[1, 0:2]),
                           "gs": fig.add_subplot(gs[1, 5:7])}
    month_axes["April"] = {"is": fig.add_subplot(gs[1, 2:4]),
                           "gs": fig.add_subplot(gs[1, 7:9])}

    month_axes["May"] = {"is": fig.add_subplot(gs[2, 0:2]),
                         "gs": fig.add_subplot(gs[2, 5:7])}
    month_axes["June"] = {"is": fig.add_subplot(gs[2, 2:4]),
                          "gs": fig.add_subplot(gs[2, 7:9])}

    month_axes["July"] = {"is": fig.add_subplot(gs[3, 0:2]),
                          "gs": fig.add_subplot(gs[3, 5:7])}
    month_axes["August"] = {"is": fig.add_subplot(gs[3, 2:4]),
                            "gs": fig.add_subplot(gs[3, 7:9])}

    month_axes["September"] = {"is": fig.add_subplot(gs[4, 0:2]),
                               "gs": fig.add_subplot(gs[4, 5:7])}
    month_axes["October"] = {"is": fig.add_subplot(gs[4, 2:4]),
                             "gs": fig.add_subplot(gs[4, 7:9])}

    month_axes["November"] = {"is": fig.add_subplot(gs[5, 0:2]),
                              "gs": fig.add_subplot(gs[5, 5:7])}
    month_axes["December"] = {"is": fig.add_subplot(gs[5, 2:4]),
                              "gs": fig.add_subplot(gs[5, 7:9])}

    # # Add two colour bar axes, one for is and one for gs
    cbar_axes["is"] = fig.add_subplot(gs[6, 1:3])
    cbar_axes["gs"] = fig.add_subplot(gs[6, 6:8])

    return month_axes, cbar_axes


def apply_subplot_formatting(axes):
    """
    :param axes: matplotlib.axes: All of the smaller monthly axis requiring formatting
    """

    x_lim = (-10, 10)
    y_lim = (-10, 10)

    label_font_size = 10
    title_font_size = 10

    for month, month_dict_item in axes.items():
        for subplot_type, ax in month_dict_item.items():

            ax.set_title(month[:3] + "; " + subplot_type.upper(), fontsize=title_font_size)
            ax.set_xlabel("By [nT] (GSM)", fontsize=label_font_size)
            ax.set_ylabel("Bz [nT] (GSM)", fontsize=label_font_size)

            ax.set_xlim(x_lim)
            ax.xaxis.set_major_locator(MultipleLocator(5))
            ax.xaxis.set_minor_locator(MultipleLocator(1))

            ax.set_ylim(y_lim)
            ax.yaxis.set_major_locator(MultipleLocator(5))
            ax.yaxis.set_minor_locator(MultipleLocator(1))

            ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5, color='black')
            ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2, color='black')


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    if local_testing:
        station = "rkn"

        # Note: year, month, and day don't matter for local testing
        fig = occ_imf_variation(station=station, year=2011, day_range=(12, 12), hour_range=None,
                                gate_range=(10, 30), beam_range=(6, 8), freq_range=(11, 13),
                                local_testing=local_testing)

        plt.show()


    else:
        station = "dcn"
        even_odd_days = None
        year = 2019
        freq_range = (8, 14)

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        fig = occ_imf_variation(station=station, year=year, day_range=None, hour_range=None,
                                gate_range=(10, 30), beam_range=(6, 8), freq_range=freq_range,
                                local_testing=local_testing, even_odd_days=even_odd_days)

        out_fig = out_dir + "/occ_imf_variation_" + station + \
                  "_" + str(year) + "_" + str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz"

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
