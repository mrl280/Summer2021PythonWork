import math
import os
import pathlib
import pydarn

from matplotlib.ticker import MultipleLocator
from cartopy.util import add_cyclic_point
from matplotlib import pyplot as plt
from pydarn import radar_fov, SuperDARNRadars

import datetime as datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import numpy as np

from DataAnalysis.EchoOccurrence.lib.only_keep_overlap import only_keep_overlap
from occ_clock_diagram import add_lat_labels_to_clock_diagram

from lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from lib.cm.modified_jet import modified_jet
from lib.add_mlt_to_df import add_mlt_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def two_station_vel_comparison(station1, station2, start_epoch, end_epoch,
                               gate_range=None, beam_range=None, freq_range=None,
                               local_testing=False):
    """

    Produce a series of plots that is meant to showcase the elevation/velocity comparison between two SuperDARN radars.

    The plan was to use find_high_vel_overlap_events() to find event of interest and then to feed those events into
     this program.

    Contour comparison plots will have station1 along the x-axis.

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            Maxwell holds all SuperDARN data.
        - Only considers 45 km data.
            (a warning will be printed if other spatial resolution data is stripped from the dataset)
        - To check which fitACF program is being used, refer to the data readers in lib.data_getters

    :param station1: str:
            The first radar station to consider, as 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param station2: str:
            The second radar station to consider, again as a 3 character string.
    :param start_epoch: int:
            The starting epoch of the event to consider
    :param end_epoch:
            The ending epoch of the event to consider.  Usually this is 4 hours after the starting epoch.
    :param gate_range: (int, int) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that early gates probably are not seeing F region echoes, so you probably don't want to consider them
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (int, int) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param freq_range: (float, float) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    """

    seconds_in_an_hour = 3600
    t_diffs = {station1: 0.000,
               station2: 0.000}
    param_ranges = {'vel': (-600, 600),
                    'height': (100, 300)}

    if end_epoch < start_epoch:
        raise Exception("two_station_vel_comparison(): start epoch must be before end epoch.")

    # When writing this I was planning on 4 hour events, warn if this is not what we were given
    four_hours_worth_of_seconds = 4 * seconds_in_an_hour
    if ((end_epoch - start_epoch) - four_hours_worth_of_seconds) > 5:
        warnings.warn("two_station_vel_comparison() was designed to run on four-hour events, "
                      "events of different lengths might not work")

    # Grab the radars info from the hardware files
    all_radars_info = SuperDARNRadars()
    station1_stid = pydarn.read_hdw_file(station1).stid
    station2_stid = pydarn.read_hdw_file(station2).stid
    first_radars_info = all_radars_info.radars[station1_stid]
    second_radars_info = all_radars_info.radars[station2_stid]

    reference_hemisphere = first_radars_info.hemisphere
    if second_radars_info.hemisphere != reference_hemisphere:
        raise Exception("two_station_vel_comparison(): " + station1 + " and " + station2
                        + " are not in the same hemisphere.  We can't compare them because they have no overlap.")

    starting_datetime = datetime.datetime.fromtimestamp(start_epoch)
    ending_datetime = datetime.datetime.fromtimestamp(end_epoch)

    year_range = (starting_datetime.year, ending_datetime.year)
    month_range = (starting_datetime.month, ending_datetime.month)
    day_range = (starting_datetime.day, ending_datetime.day)

    # print("Getting SuperDARN data for " + station1.upper())
    # df1 = get_data_handler(station1, year_range=year_range, month_range=month_range, day_range=day_range,
    #                        gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
    #                        local_testing=local_testing, even_odd_days=None)
    # df1 = only_keep_45km_res_data(df1)
    #
    # print("Restricting " + station1.upper() + "'s data - only keep data for those cells that overlap with "
    #       + station2.upper())
    # df1 = only_keep_overlap(station=station1, df=df1, other_station=station2, gate_min=gate_range[0])
    #
    # print("Getting SuperDARN data for " + station2.upper())
    # df2 = get_data_handler(station1, year_range=year_range, month_range=month_range, day_range=day_range,
    #                        gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
    #                        local_testing=local_testing, even_odd_days=None)
    # df2 = only_keep_45km_res_data(df2)
    #
    # print("Restricting " + station2.upper() + "'s data - only keep data for those cells that overlap with "
    #       + station1.upper())
    # df2 = only_keep_overlap(station=station2, df=df2, other_station=station1, gate_min=gate_range[0])

    # station1_stid = pydarn.read_hdw_file(station1).stid
    # df1 = add_decimal_hour_to_df(df=df1, time_units='ut', stid=station1_stid,
    #                              date_time_est=(df1['datetime'].iat[0]).to_pydatetime())
    #
    # station2_stid = pydarn.read_hdw_file(station2).stid
    # df2 = add_decimal_hour_to_df(df=df2, time_units='ut', stid=station2_stid,
    #                              date_time_est=(df2['datetime'].iat[0]).to_pydatetime())

    print("     Preparing the figure...")
    fig = plt.figure(figsize=[20, 8], constrained_layout=True, dpi=300)
    axes = add_axes(fig=fig, reference_hemisphere=reference_hemisphere)

    apply_subplot_formatting(axes=axes, station1=station1, station2=station2,
                             reference_hemisphere=reference_hemisphere, t_diffs=t_diffs,
                             gate_range=gate_range, vel_range=param_ranges['vel'], height_range=param_ranges['height'])

    # title_figure(fig=fig, station=station, year=year, time_units=time_units, beam_range=beam_range,
    #              gate_range=gate_range, freq_range=freq_range, hour_range=hour_range, day_range=day_range,
    #              title_fontsize=title_fontsize)

    return fig


def apply_subplot_formatting(axes, station1, station2, reference_hemisphere, t_diffs=None,
                             gate_range=(0, 74), vel_range=(-600, 600), height_range=(0, 300)):
    """


    :param axes: dictionary of matplotlib.axes:
            A dictionary containing all the axes items.
    :param station1: str:
            The first radar station to consider, as 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param station2: str:
            The second radar station to consider, again as a 3 character string.
    :param reference_hemisphere: pydarn hemisphere object:
            The radars' hemisphere

    :param t_diffs: Dictionary of floats: (Optional, default is {station1: 0.000, station2: 0.000})
            The extra time delay added in when computing adjusted elevation angle, in microseconds.
    :param gate_range: (int, int) (Optional; default is (0, 74))
            See two_station_vel_comparison() docstring
    :param height_range: (float, float) (Optional; default is (0, 300))
            The height range in km
    :param vel_range: (float, float) (Optional; default is (-600, 600))
            The velocity range in m/s

    """

    if t_diffs is None:
        t_diffs = {station1: 0.000, station2: 0.000}

    label_font_size = 10
    title_font_size = 10
    bisector_colour = "purple"

    for axis_key, axis_item in axes.items():

        if axis_key == 'station1_rti' or axis_key == 'station2_rti':
            # We are formatting the large RTI plots
            for subplot_type, ax in axis_item.items():

                y_lim = gate_range

                ax.set_ylim(y_lim)
                ax.yaxis.set_major_locator(MultipleLocator(10))  # Every 10 gates
                ax.yaxis.set_minor_locator(MultipleLocator(1))

                # ax.set_xlim(x_lim)
                ax.xaxis.set_major_locator(MultipleLocator(0.5))  # Every half hour
                ax.xaxis.set_minor_locator(MultipleLocator(0.1))

                ax.grid(b=True, which='major', axis='both', linewidth=1, linestyle='-')
                ax.grid(b=True, which='minor', axis='both', linewidth=0.4, linestyle='--')

                ax.set_ylabel("Range Gate", fontsize=label_font_size)

                if axis_key == 'station1_rti':
                    station_mnemonic = station1.upper()
                    t_diff = t_diffs[station1]
                else:
                    station_mnemonic = station2.upper()
                    t_diff = t_diffs[station2]

                if subplot_type == "vel":
                    ax.set_title(station_mnemonic.upper() + " Velocities", fontsize=title_font_size)
                elif subplot_type == "adjElv":
                    ax.set_title(station_mnemonic.upper() + " Adjusted Elevation Angles, tdiff=" + str(t_diff)
                                 + " \u03BCs", fontsize=title_font_size)
                    ax.set_xlabel("UT Time [hour]", fontsize=label_font_size)

        elif axis_key == 'map':
            # We are formatting the large map
            lat_extreme = reference_hemisphere.value * 60  # deg
            ax = axis_item
            ax.set_extent([-180, 180, reference_hemisphere.value * 90, lat_extreme], crs=ccrs.PlateCarree())
            ax.gridlines()
            ax.add_feature(cfeature.OCEAN)
            ax.add_feature(cfeature.LAND)

        else:
            # We are formatting the comparison plots

            for subplot_type, ax in axis_item.items():

                if subplot_type == "vel":
                    ax.set_ylim(vel_range)
                    ax.set_xlim(vel_range)

                    ax.set_xlabel(station1.upper() + " Velocities [m/s]")
                    ax.set_ylabel(station2.upper() + " Velocities [m/s]")

                    ax.yaxis.set_major_locator(MultipleLocator(600))
                    ax.xaxis.set_major_locator(MultipleLocator(600))
                    ax.yaxis.set_minor_locator(MultipleLocator(100))
                    ax.xaxis.set_minor_locator(MultipleLocator(100))
                    if axis_key == "full_event_comparison":
                        plt.sca(ax)
                        plt.xticks([-600, -400, -200, 0, 200, 400, 600])
                        plt.yticks([-600, -400, -200, 0, 200, 400, 600])

                elif subplot_type == "height":
                    ax.set_ylim(height_range)
                    ax.set_xlim(height_range)

                    ax.set_xlabel(station1.upper() + " Virtual Heights [km]")
                    ax.set_ylabel(station2.upper() + " Virtual Heights [km]")

                    ax.yaxis.set_major_locator(MultipleLocator(100))
                    ax.xaxis.set_major_locator(MultipleLocator(100))
                    ax.yaxis.set_minor_locator(MultipleLocator(50))
                    ax.xaxis.set_minor_locator(MultipleLocator(50))
                    if axis_key == "full_event_comparison":
                        plt.sca(ax)
                        plt.xticks([100, 150, 200, 250, 300])
                        plt.yticks([100, 150, 200, 250, 300])

                ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
                ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
                ax.plot(ax.get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
                ax.plot([0, 0], ax.get_xlim(), linestyle='-', linewidth=0.5, color='black')
                ax.plot([ax.get_ylim()[0], ax.get_ylim()[1]], [ax.get_xlim()[0], ax.get_xlim()[1]],
                        linestyle='--', linewidth=1, color=bisector_colour)


def add_axes(fig, reference_hemisphere):
    """
    :param fig: matplotlib.pyplot.figure:
            The figure to draw the axes on
    :param reference_hemisphere: pydarn hemisphere object:
            The radars' hemisphere

    :return axes: dictionary of dictionary of axes:
            A dictionary containing all the axes items.
    """

    gs = fig.add_gridspec(ncols=11, nrows=4)

    axes = dict()  # Remember that splices don't include last index

    # We will have two large range-time plots - one for each station
    axes['station1_rti'] = {"vel": fig.add_subplot(gs[0, 0:3]), "adjElv": fig.add_subplot(gs[1, 0:3])}
    axes['station2_rti'] = {"vel": fig.add_subplot(gs[0, 6:9]), "adjElv": fig.add_subplot(gs[1, 6:9])}

    # We will have a map that shows each radar's fan and highlights the overlap
    if reference_hemisphere.value == 1:
        proj = ccrs.NorthPolarStereo()
    elif reference_hemisphere.value == -1:
        proj = ccrs.SouthPolarStereo()
    else:
        raise Exception("Hemisphere not recognized.")
    axes['map'] = fig.add_subplot(gs[0:2, 3:6], projection=proj)

    # We will have a couple larger plots that show select comparisons for the whole event duration
    axes['full_event_comparison'] = {"vel": fig.add_subplot(gs[2:4, 0:3]), "height": fig.add_subplot(gs[2:4, 3:6])}

    # We will have some small plots that show select comparisons for parts of the event
    axes['interval1'] = {"vel": fig.add_subplot(gs[0, 9]), "height": fig.add_subplot(gs[0, 10])}
    axes['interval2'] = {"vel": fig.add_subplot(gs[1, 9]), "height": fig.add_subplot(gs[1, 10])}
    axes['interval3'] = {"vel": fig.add_subplot(gs[2, 9]), "height": fig.add_subplot(gs[2, 10])}
    axes['interval4'] = {"vel": fig.add_subplot(gs[3, 9]), "height": fig.add_subplot(gs[3, 10])}

    return axes


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    if local_testing:
        station1 = "dcn"
        station2 = "mcm"

        start_epoch = 1321070400  # Saturday, November 12, 2011 4:00:00 AM = 1321070400
        end_epoch = 1321084800  # Saturday, November 12, 2011 8:00:00 AM = 1321084800

        gate_range = (20, 74)
        beam_range = (0, 15)

        # Note: year, month, and day don't matter for local testing
        fig = two_station_vel_comparison(station1=station1, station2=station2,
                                         start_epoch=start_epoch, end_epoch=end_epoch,
                                         gate_range=gate_range, beam_range=beam_range, freq_range=None,
                                         local_testing=local_testing)

        plt.show()


    else:
        pass
        # station = "dcn"
        # even_odd_days = None
        # years = [2019, 2020, 2021]
        # plot_type = 'contour'
        # freq_range = (8, 11)
        #
        # loc_root = str((pathlib.Path().parent.absolute()))
        # out_dir = loc_root + "/out"
        #
        # for year in years:
        #     fig = occ_imf_variation(station=station, year=year, day_range=(1, 15), hour_range=None,
        #                             gate_range=(10, 30), beam_range=(6, 8), freq_range=freq_range,
        #                             local_testing=local_testing, plot_type=plot_type, even_odd_days=even_odd_days)
        #
        #     out_fig = out_dir + "/occ_imf_variation_" + station + \
        #               "_" + str(year) + "_" + str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz_" + plot_type
        #
        #     print("Saving plot as " + out_fig)
        #     fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
