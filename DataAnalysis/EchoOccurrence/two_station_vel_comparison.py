import calendar
import os
import pathlib
import pydarn

import cartopy.crs as ccrs
import numpy as np

from matplotlib import pyplot as plt
from pydarn import SuperDARNRadars

from lib.adjust_velocities_los import adjust_velocities_los
from lib.get_data_handler import get_data_handler
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from two_station_vel_comparison_file_events import format_data_axis, format_map_axis
from two_station_vel_comparison_single_event import complete_map_subplot, complete_comparison_plot
from lib.build_two_radar_matched_data import build_two_radar_matched_data
from lib.data_getters.input_checkers import *


def two_station_vel_comparison(station1, station2, year, month_range=None, day_range=None, even_odd_days=None,
                               use_only_coinciding_beams=True, station1_ref_beam=None, station2_ref_beam=None,
                               gate_range1=None, beam_range1=None, gate_range2=None, beam_range2=None,
                               freq_range=None, plot_type='contour', local_testing=False):
    """
    Produce a velocity comparison plot comparing two SuperDARN radars.
    This program runs on the provided period of time (year, month, and day range).

    Other version of similar programs that you may want to check out include:
        - two_station_vel_comparison_file_events() runs on a bunch of events summarized in a file.  This file of events
            is often produced by find_high_vel_overlap_events().
        - two_station_vel_comparison_single_event() runs on a single event from start_epoch to end_epoch.

    Contour/pixel comparison plots will have station1 along the x-axis.

    Originally this program was designed to use the full area of mutual overlap, and take the cosine of all velocity
     measurements to orient them along the same reference beam (same line-of-sight).  However, from Koustov:
     "E region velocity does not follow the cosine rule. How well this rule works in the F region is also an issue".
    So, there are two options:
        - If use_only_coinciding_beams=True, then gate and beam range restrictions will be applied, but velocities will
           not be modified.
        - If use_only_coinciding_beams=False, then gate and beam range restrictions will still be applied, but
           velocities will be multiplied by cosine of their angular separation from their reference beam and reference
           beams will be added to the map plot.

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
    :param year: int:
            The year to consider.  At most this program can produce a plot for a year's worth of comparisons.
            Anything more and it would take too long to run.
    :param month_range: (<int>, <int>) (optional):
            Inclusive. The months of the year to consider.  If omitted (or None), then all days will be considered.
    :param day_range: (<int>, <int>) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param even_odd_days: (optional; default is None)
            'even': only even days are read in
            'odd': only odd days are read in
            None: all days are read in

    :param use_only_coinciding_beams: bool (optional; default is True).
            - If True; apply gate and beam range restrictions, but don't adjust velocities
            - If False; still apply gate and beam range restrictions, but also modify velocities (requires reference
            beams be passed in)
    :param station1_ref_beam: int:
            Only referenced if use_only_coinciding_beams == False.
            The beam used to adjust station1 velocity measurements so they are as if they are looking along this beam.
            We need to modify LoS velocities so that line-of-sight velocity measurements from different directions can
            be compared directly - to do this select two beam (one from each station) that point the same direction and
            adjust all velocities so it is as if all velocity measurements are looking along this same line-of-sight.
    :param station2_ref_beam: int:
            Only referenced if use_only_coinciding_beams == False.
            The beam used to adjust station2 velocity measurements so they are as if they are looking along this beam.

    :param gate_range1: (int, int) (optional):
            Inclusive. The gate range for the first station.  If omitted (or None), then all the gates will be considered.
            Note that early gates probably are not seeing F region echoes, so you probably don't want to consider them
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range1: (int, int) (optional):
            Inclusive. The beam range for the first station.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param gate_range2: (int, int) (optional):
            Inclusive. The gate range for the second station.  If omitted (or None), then all the gates will be considered.
            Note that early gates probably are not seeing F region echoes, so you probably don't want to consider them
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range2: (int, int) (optional):
            Inclusive. The beam range for the second station.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param freq_range: (float, float) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.

    :param plot_type: str (optional; default is 'contour'):
            The type of plot, either 'contour' or 'pixel', default is 'contour'
    :param local_testing: bool (optional; default is False):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.

    :return fig:  matplotlib.pyplot.figure:
            The figure. It can then be modified, added to, printed out, or saved in whatever format is desired.
    """

    vel_range = (-1000, 1000)
    n_vel_bins = 50

    count_min = 1  # Minimum number of points needed in a spatial/temporal 'clump' in order plot the matched point
    title_fontsize = 18
    time_interval_s = 180  # Controls spatial resolution of matched data

    if use_only_coinciding_beams == False:
        # We would like to use the entire overlap, make sure we have reference beams
        if station1_ref_beam is None or station2_ref_beam is None:
            raise Exception("two_station_vel_comparison(): In order to use the entire overlap, you must pass in "
                            "reference beams.")

    year = check_year(year=year)
    month_range = check_month_range(month_range=month_range)
    day_range = check_day_range(day_range=day_range)
    if month_range[0] == month_range[1]:
        month_string = calendar.month_name[month_range[0]]
    else:
        month_string = calendar.month_name[month_range[0]] + " to " + calendar.month_name[month_range[1]]

    # Compute parameter edges used in the contour/pixel comparisons
    vel_edges = np.linspace(vel_range[0], vel_range[1], num=(n_vel_bins + 1))

    # Grab the radars info from the hardware files
    all_radars_info = SuperDARNRadars()
    station1_stid = pydarn.read_hdw_file(station1).stid
    station2_stid = pydarn.read_hdw_file(station2).stid
    first_radars_info = all_radars_info.radars[station1_stid]
    second_radars_info = all_radars_info.radars[station2_stid]

    # Check input ranges, make sure they work for both radars
    freq_range = check_freq_range(freq_range=freq_range)
    gate_range1 = check_gate_range(gate_range=gate_range1, hdw_info=first_radars_info.hardware_info)
    gate_range2 = check_gate_range(gate_range=gate_range2, hdw_info=second_radars_info.hardware_info)
    beam_range1 = check_beam_range(beam_range=beam_range1, hdw_info=first_radars_info.hardware_info)
    beam_range2 = check_beam_range(beam_range=beam_range2, hdw_info=second_radars_info.hardware_info)

    # Make sure that both of the radars provided are in the same hemisphere
    reference_hemisphere = first_radars_info.hemisphere
    if second_radars_info.hemisphere != reference_hemisphere:
        raise Exception("two_station_vel_comparison(): " + station1 + " and " + station2
                        + " are not in the same hemisphere.  We can't compare them because they have no overlap.")

    if reference_hemisphere.value == 1:
        proj = ccrs.NorthPolarStereo()
    elif reference_hemisphere.value == -1:
        proj = ccrs.SouthPolarStereo()
    else:
        raise Exception("Hemisphere not recognized.")

    print("Preparing the figure...")
    fig = plt.figure(figsize=[12, 6], constrained_layout=True, dpi=300)
    gs = fig.add_gridspec(ncols=2, nrows=1)

    print("Working on the map...")
    map_axis = fig.add_subplot(gs[0, 0], projection=proj)
    format_map_axis(ax=map_axis, reference_hemisphere=reference_hemisphere)
    complete_map_subplot(ax=map_axis, station1=station1, station2=station2,
                         station1_ref_beam=station1_ref_beam, station2_ref_beam=station2_ref_beam,
                         gate_range1=gate_range1, gate_range2=gate_range2,
                         beam_range1=beam_range1, beam_range2=beam_range2,
                         use_only_coinciding_beams=use_only_coinciding_beams)

    data_axis = fig.add_subplot(gs[0, 1])
    format_data_axis(ax=data_axis, vel_range=vel_range, station1=station1, station2=station2)

    title_figure(fig=fig, station1=station1, station2=station2, year=year, month_string=month_string,
                 gate_range1=gate_range1, gate_range2=gate_range2, beam_range1=beam_range1, beam_range2=beam_range2,
                 freq_range=freq_range, title_fontsize=title_fontsize)

    print("     Retrieving SuperDARN data for station " + station1.upper() + "...")
    df1 = get_data_handler(station1, year_range=(year, year), month_range=month_range, day_range=day_range,
                           gate_range=gate_range1, beam_range=beam_range1, freq_range=freq_range,
                           even_odd_days=even_odd_days, local_testing=local_testing)
    df1 = only_keep_45km_res_data(df=df1)

    print("     Retrieving SuperDARN data for station " + station2.upper() + "...")
    df2 = get_data_handler(station2, year_range=(year, year), month_range=month_range, day_range=day_range,
                           gate_range=gate_range2, beam_range=beam_range2, freq_range=freq_range,
                           even_odd_days=even_odd_days, local_testing=local_testing)
    df2 = only_keep_45km_res_data(df=df2)

    if not use_only_coinciding_beams:
        # We want to use the whole overlap range - cosine adjust velocities
        df1 = adjust_velocities_los(station=station1, df=df1, ref_beam=station1_ref_beam)
        df2 = adjust_velocities_los(station=station2, df=df2, ref_beam=station2_ref_beam)

    df1 = df1.loc[(np.abs(df1['v']) > 50) & (np.abs(df1['v']) < 1000)]  # Remove extreme values
    df2 = df2.loc[(np.abs(df1['v']) > 50) & (np.abs(df2['v']) < 1000)]  # Low velocities are usually not of interest

    df1.reset_index(drop=True, inplace=True)
    df2.reset_index(drop=True, inplace=True)

    print("Building a matched dataframe...")
    matched_df = build_two_radar_matched_data(station1=station1, df1=df1, gate_range=gate_range1,
                                              beam_range=beam_range1, station2=station2, df2=df2,
                                              other_gate_range=gate_range2, other_beam_range=beam_range2,
                                              time_interval_s=time_interval_s)
    matched_df = matched_df.loc[(matched_df['count1'] >= count_min) & (matched_df['count2'] >= count_min)]

    print("Completing the large comparison plot...")
    complete_comparison_plot(fig=fig, ax=data_axis, matched_df=matched_df,
                             param='v', param_edges=vel_edges, cbar=True, plot_type=plot_type)

    return fig


def title_figure(fig, station1, station2, year, month_string, gate_range1, beam_range1, gate_range2, beam_range2,
                 freq_range, title_fontsize):
    """
    :param fig: The figure to title
    :param station1: str:
            The first radar station to consider, as 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param station2: str:
            The second radar station to consider, again as a 3 character string.
    :param year: int:
            The year we are considering
    :param gate_range1: See two_station_vel_comparison() docstring
    :param beam_range1: See two_station_vel_comparison() docstring
    :param gate_range2: See two_station_vel_comparison() docstring
    :param beam_range2: See two_station_vel_comparison() docstring
    :param freq_range: See two_station_vel_comparison() docstring
    :param title_fontsize:
    """

    first_station_string = "Station 1: " + station1.upper() + \
                           " (Allowed Beams " + str(beam_range1[0]) + "-" + str(beam_range1[1]) + \
                           "; Allowed Gates: " + str(gate_range1[0]) + "-" + str(gate_range1[1]) + ")"
    second_station_string = "Station 2: " + station2.upper() + \
                            " (Allowed Beams " + str(beam_range2[0]) + "-" + str(beam_range2[1]) + \
                            "; Allowed Gates: " + str(gate_range2[0]) + "-" + str(gate_range2[1]) + ")"

    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"

    fig.suptitle(month_string + " " + str(year) + " Velocity Comparison; " + freq_string +
                 "\n" + first_station_string +
                 "\n" + second_station_string +
                 "\n" + "Produced by " + str(os.path.basename(__file__)), fontsize=title_fontsize)


if __name__ == '__main__':
    """ Testing """

    local_testing = False
    use_only_coinciding_beams = True  # See note from Koustov regarding validity of cosine approach

    if local_testing:
        station1 = "dcn"
        gate_range1 = (20, 74)
        beam_range1 = (14, 15)
        station1_ref_beam = 15
        station2 = "mcm"
        gate_range2 = (20, 74)
        beam_range2 = (8, 8)
        station2_ref_beam = 8

        # station1 = "inv"
        # gate_range1 = (20, 74)
        # beam_range1 = (12, 14)
        # station1_ref_beam = 13
        # station2 = "kod"
        # gate_range2 = (20, 74)
        # beam_range2 = (7, 8)
        # station2_ref_beam = 7

        # station1 = "pgr"
        # gate_range1 = (20, 74)
        # beam_range1 = (5, 6)
        # station1_ref_beam = 6
        # station2 = "cvw"
        # gate_range2 = (20, 74)
        # beam_range2 = (15, 15)
        # station2_ref_beam = 15

        # station1 = "inv"
        # gate_range1 = (15, 74)
        # beam_range1 = (15, 15)
        # station1_ref_beam = 15
        # station2 = "cly"
        # gate_range2 = (15, 74)
        # beam_range2 = (4, 5)
        # station2_ref_beam = 4

        # Note: year, month, and day don't matter for local testing
        fig = two_station_vel_comparison(station1=station1, station2=station2, year=2011,
                                         use_only_coinciding_beams=use_only_coinciding_beams,
                                         station1_ref_beam=station1_ref_beam,
                                         station2_ref_beam=station2_ref_beam,
                                         gate_range1=gate_range1, beam_range1=beam_range1,
                                         gate_range2=gate_range2, beam_range2=beam_range2,
                                         freq_range=None, plot_type='contour',
                                         local_testing=local_testing)

        plt.show()


    else:

        station1 = "dcn"
        gate_range1 = (20, 74)
        beam_range1 = (14, 15)
        station1_ref_beam = 15
        station2 = "mcm"
        gate_range2 = (20, 74)
        beam_range2 = (8, 8)
        station2_ref_beam = 8

        # station1 = "inv"
        # gate_range1 = (20, 74)
        # beam_range1 = (12, 14)
        # station1_ref_beam = 13
        # station2 = "kod"
        # gate_range2 = (20, 74)
        # beam_range2 = (7, 8)
        # station2_ref_beam = 7

        # station1 = "pgr"
        # gate_range1 = (20, 74)
        # beam_range1 = (5, 6)
        # station1_ref_beam = 6
        # station2 = "cvw"
        # gate_range2 = (20, 74)
        # beam_range2 = (15, 15)
        # station2_ref_beam = 15

        # station1 = "inv"
        # gate_range1 = (15, 74)
        # beam_range1 = (15, 15)
        # station1_ref_beam = 15
        # station2 = "cly"
        # gate_range2 = (15, 74)
        # beam_range2 = (4, 5)
        # station2_ref_beam = 4

        year = 2019
        freq_range = (8, 11)  # Koustov suggests < 11
        plot_type = 'contour'
        month_range = None
        day_range = None
        even_odd_days = 'odd'

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        fig = two_station_vel_comparison(station1=station1, station2=station2, year=year,
                                         month_range=month_range, day_range=day_range,
                                         use_only_coinciding_beams=use_only_coinciding_beams,
                                         station1_ref_beam=station1_ref_beam,
                                         station2_ref_beam=station2_ref_beam,
                                         gate_range1=gate_range1, beam_range1=beam_range1,
                                         gate_range2=gate_range2, beam_range2=beam_range2,
                                         freq_range=freq_range, plot_type=plot_type, even_odd_days=even_odd_days,
                                         local_testing=local_testing)

        out_fig = out_dir + "/two_station_comparison-" + station1 + "_" + station2 + "-" + str(year)

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
