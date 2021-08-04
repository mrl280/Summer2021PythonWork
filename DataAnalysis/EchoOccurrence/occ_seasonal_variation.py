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
                           time_units='lt', local_testing=False, even_odd_days=None):
    """

    Produces an occurrence rate versus year plot that is meant to showcase the seasonal variation in echo occurrence.

    There is one subplot for all of hour_range, and then subplots for each 6 h time sector
     Time sectors are used to automatically consider predefined hour ranges that are of common interest
     "Day" (12 +/- 3 h), "Dusk" (18 +/- 3 h), "Night" (24 +/- 3 h), are "Dawn" (6 +/- 3 h). Time sectors do not overlap.

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
            This parameter only effects the all hours plot.
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
    :param time_units: str:
            The provided hour_range is assumed to be in these time units.  Default is 'lt'.
                'ut' for universal time
                'mlt' for magnetic local time
                'lt' for local time (based on longitude)
                'lst' for local standard time (based on time zones)
                'ast' for apparent solar time (based on the apparent angular motion of the sun across the sky)
    :param even_odd_days: (optional; default is None)
            'even': only even days are read in
            'odd': only odd days are read in
            None: all days are read in
    :param local_testing: bool (optional): default is False.
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.

    :return: pandas.DataFrame, matplotlib.pyplot.figure
            The dataframe used and the figure created.
            The figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    y_lim = [0, 0.8]
    title_font_size = 12

    time_units = check_time_units(time_units)
    year_range = check_year_range(year_range)
    freq_range = check_freq_range(freq_range)
    hour_range = check_hour_range(hour_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    print("Retrieving data...")
    df = get_data_handler(station, year_range=year_range, month_range=None, day_range=day_range,
                          gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                          occ_data=True, local_testing=local_testing, even_odd_days=even_odd_days)
    df = only_keep_45km_res_data(df)

    # Add decimal hour to df in whatever units were requested
    # Use the middle of the mid year as magnetic field estimate
    mid_year = int(year_range[0] + (year_range[1] - year_range[0]) / 2)
    date_time_est, _ = build_datetime_epoch(mid_year, 6, 15, 0)
    df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id, date_time_est=date_time_est)

    print("Preparing the figure...")
    if year_range[1] == year_range[0]:
        year_string = str(year_range[0])
    else:
        year_string = str(year_range[0]) + " to " + str(year_range[1])
    beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    gate_string = "Gates " + str(gate_range[0]) + "-" + str(gate_range[1])
    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"

    fig = plt.figure(figsize=[12, 8], dpi=300, constrained_layout=True)
    axes = add_axes(fig=fig)
    format_subplots(axes=axes, x_lim=(year_range[0], year_range[1] + 1), y_lim=y_lim)
    fig.suptitle(year_string + " at " + station.upper() + "; " + gate_string + "; " + beam_string + "; " + freq_string
                 + "\nProduced by " + str(os.path.basename(__file__)), fontsize=18)

    print("Computing UT decimal years...")
    months_in_a_year = 12
    days_in_a_year = 365
    hours_in_a_year = 8760

    decimal_years = []
    for i in range(len(df)):
        datetime_obj_here = df['datetime'].iat[i]
        decimal_hour_here = df[time_units].iat[i]

        decimal_year_here = datetime_obj_here.year + (datetime_obj_here.month - 1) / months_in_a_year + \
                            (datetime_obj_here.day - 1) / days_in_a_year + decimal_hour_here / hours_in_a_year

        decimal_years.append(decimal_year_here)

    df['decimal_year'] = np.asarray(decimal_years)

    # Compute year_edges
    bins_per_year = 180
    n_bins_x = ((year_range[1] + 1) - year_range[0]) * bins_per_year
    year_slice_edges = np.linspace(year_range[0], (year_range[1] + 1), num=(n_bins_x + 1))

    day_hour_range = (9, 15)  # From 9 AM to 3 PM
    dusk_hour_range = (15, 21)  # From 3 PM to 9 PM
    night_hour_range = (21, 3)  # Note that this one is different!  From 9 PM to 3 AM (21 - 2)
    dawn_hour_range = (3, 9)  # From 3 AM to 9 AM

    # Add in custom hour range data
    df_custom = df[(df[time_units] >= hour_range[0]) & (df[time_units] <= hour_range[1])]
    add_data_to_plot(df=df_custom, ax=axes['custom'], year_edges=year_slice_edges)
    title_string = "Hours " + str(hour_range[0]) + "-" + str(hour_range[1]) + " " + time_units.upper()
    axes['custom'].set_title(title_string, fontsize=title_font_size)
    axes['custom'].legend(loc='upper right', ncol=2, prop={'size': 6})

    # Add in dawn data
    time_sector = 'dawn'
    df_dawn = df[(df[time_units] >= dawn_hour_range[0]) & (df[time_units] < dawn_hour_range[1])]
    add_data_to_plot(df=df_dawn, ax=axes[time_sector], year_edges=year_slice_edges)
    title_string = time_sector.capitalize() + " (3 AM to 9 AM; in " + time_units.upper() + ")"
    axes[time_sector].set_title(title_string, fontsize=title_font_size)

    # Add in day data
    time_sector = 'day'
    df_day = df[(df[time_units] >= day_hour_range[0]) & (df[time_units] < day_hour_range[1])]
    add_data_to_plot(df=df_day, ax=axes[time_sector], year_edges=year_slice_edges)
    title_string = time_sector.capitalize() + " (9 AM to 3 PM; in " + time_units.upper() + ")"
    axes[time_sector].set_title(title_string, fontsize=title_font_size)

    # Add in dusk data
    time_sector = 'dusk'
    df_dusk = df[(df[time_units] >= dusk_hour_range[0]) & (df[time_units] < dusk_hour_range[1])]
    add_data_to_plot(df=df_dusk, ax=axes[time_sector], year_edges=year_slice_edges)
    title_string = time_sector.capitalize() + " (3 PM to 9 PM; in " + time_units.upper() + ")"
    axes[time_sector].set_title(title_string, fontsize=title_font_size)

    # Add in night data
    # Note that this one has an OR condition because we are bridging data from either side of the midnight line
    time_sector = 'night'
    df_night = df.loc[(df[time_units] >= night_hour_range[0]) | (df[time_units] < night_hour_range[1])]
    add_data_to_plot(df=df_night, ax=axes[time_sector], year_edges=year_slice_edges)
    title_string = time_sector.capitalize() + " (9 PM to 3 AM; in " + time_units.upper() + ")"
    axes[time_sector].set_title(title_string, fontsize=title_font_size)

    # Add in solar fux and sunspot data.  We want to know if echo occurrence rates correlates with solar activity
    plot_solar_flux_data(ax=axes['f10.7'], year_range=(year_range[0], year_range[1] + 1),
                         smoothing_window_size_in_days=30)
    axes['f10.7'].set_title("Solar Flux Data from Penticton", fontsize=title_font_size)

    plot_sunspot_data(ax=axes['sunspot'], year_range=(year_range[0], year_range[1] + 1),
                      smoothing_window_size_in_days=30)
    axes['sunspot'].set_title("Sunspot Data from Belgium", fontsize=title_font_size)

    return df, fig


def add_data_to_plot(df, ax, year_edges):
    """
    :param year_edges: numpy.array:
            Array/list of the year edges including the closing edge
    :param axes: matplotlib.axes:
            The axes to draw on
    """

    n_bins_x = len(year_edges) - 1
    delta_year_slice = year_edges[1] - year_edges[0]

    # Compute bin centers
    bin_xwidth = delta_year_slice
    bin_xcenters = year_edges[1:] - bin_xwidth / 2

    # Loop through the time bins and compute occurrence rate
    occurrence_data_is = np.empty(shape=n_bins_x)
    occurrence_data_gs = np.empty(shape=n_bins_x)

    for year_idx, year_slice_start in enumerate(year_edges):
        if year_slice_start == year_edges[-1]:
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
        ax.plot(bin_xcenters, occurrence_data_is, 'bo', label='IS')
        ax.plot(bin_xcenters, occurrence_data_gs, 'ro', label='GS')

    else:
        # Plot as a faint line
        ax.plot(bin_xcenters, occurrence_data_is, linewidth=0.25, color="blue", linestyle='-', label='IS (Raw)')
        ax.plot(bin_xcenters, occurrence_data_gs, linewidth=0.25, color="red", linestyle='-', label='GS (Raw)')

        # Compute and plot smoothed data
        # Before we smooth, we have to remove all nan values - otherwise we get discontinuities in the smoothed data
        is_mini_df = pd.DataFrame({'bin_xcenters': bin_xcenters, 'occurrence_data': occurrence_data_is})
        gs_mini_df = pd.DataFrame({'bin_xcenters': bin_xcenters, 'occurrence_data': occurrence_data_gs})

        is_mini_df = is_mini_df.loc[is_mini_df['occurrence_data'].notna()]
        gs_mini_df = gs_mini_df.loc[gs_mini_df['occurrence_data'].notna()]

        # 30 bins is about a 60 day boxcar filter
        is_mini_df['occurrence_data'] = boxcar_smooth(is_mini_df['occurrence_data'], window_size=30)
        gs_mini_df['occurrence_data'] = boxcar_smooth(gs_mini_df['occurrence_data'], window_size=30)

        ax.plot(is_mini_df['bin_xcenters'], is_mini_df['occurrence_data'], linewidth=1, color="blue", linestyle='-',
                label='IS (Smoothed)')
        ax.plot(gs_mini_df['bin_xcenters'], gs_mini_df['occurrence_data'], linewidth=1, color="red", linestyle='-',
                label='GS (Smoothed)')


def format_subplots(axes, x_lim, y_lim):
    """
    :param axes: matplotlib.axes:
            The axes to format
    :param x_lim: (int, int): x-axis limits
    :param y_lim: (int, int): y-axis limits
    """

    label_font_size = 12

    for freq, ax in axes.items():
        if ax == axes['f10.7'] or ax == axes['sunspot']:
            continue

        ax.set_xlim(x_lim)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))

        ax.set_ylim(y_lim)
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))

        ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5, color='black')
        ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2, color='black')

        ax.set_xlabel("Year", fontsize=label_font_size)
        ax.set_ylabel("Occurrence", fontsize=label_font_size)


def add_axes(fig):
    """

    Create a sets of exes for the data.

    :param matplotlib.pyplot.figure:
        The figure to draw the axes on
    :return: Dictionary of matplotlib.pyplot.axis:
        A dictionary of data axes
    """

    gs = fig.add_gridspec(ncols=2, nrows=12)

    axes = dict()
    axes['custom'] = fig.add_subplot(gs[0:4, 0])
    axes['f10.7'] = fig.add_subplot(gs[4:8, 0])
    axes['sunspot'] = fig.add_subplot(gs[8:12, 0])

    axes['dawn'] = fig.add_subplot(gs[0:3, 1])
    axes['day'] = fig.add_subplot(gs[3:6, 1])
    axes['dusk'] = fig.add_subplot(gs[6:9, 1])
    axes['night'] = fig.add_subplot(gs[9:12, 1])

    return axes


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
        df, fig = occ_seasonal_variation(station=station, year_range=(2011, 2021), day_range=(12, 12),
                                         gate_range=(10, 30), beam_range=(6, 8), freq_range=(11, 13),
                                         time_units='lt', local_testing=local_testing)

        plt.show()


    else:
        station = "mcm"
        even_odd_days = "odd"
        year_range = (2013, 2021)
        freq_range = (8, 14)

        time_units = "lt"

        datetime_now = datetime.datetime.now()
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        _, fig = occ_seasonal_variation(station=station, year_range=year_range, day_range=None, hour_range=None,
                                        gate_range=(10, 30), beam_range=(6, 8), freq_range=freq_range,
                                        time_units=time_units, local_testing=local_testing, even_odd_days=even_odd_days)

        out_fig = out_dir + "/occ_seasonalVariation_" + station + \
                  "_" + str(year_range[0]) + "-" + str(year_range[1]) + \
                  "_" + str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz_" + time_units

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
