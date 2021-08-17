import os
import pathlib
import pydarn
import deprecation

from matplotlib.ticker import MultipleLocator
from matplotlib import pyplot as plt
from pydarn import radar_fov, SuperDARNRadars
from scipy import stats

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

from lib.two_station_overlap_events.get_select_dcn_mcm_overlap_events import get_select_dcn_mcm_overlap_events
from lib.elevation_v2 import elevation_v2
from lib.only_keep_overlap import only_keep_overlap
from lib.build_two_radar_matched_data import build_two_radar_matched_data
from lib.cm.single_colour_cmap import single_colour_cmap
from lib.compute_sd_radar_overlap import compute_df_radar_overlap
from lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.data_getters.input_checkers import *


@deprecation.deprecated(details="Elevation angle data is not suitable for this purpose, use "
                                "two_station_vel_comparison() instead.")
def two_station_vel_adjElv_comparison(station1, station2, start_epoch, end_epoch,
                                      gate_range=None, beam_range=None, freq_range=None,
                                      local_testing=False, plot_type='contour'):
    """

    # TODO: Out of all beams that overlap spatially, need to select only those beams that are in the same
       direction/azimuth.

    Produce a series of plots that is meant to showcase the elevation/velocity comparison between two SuperDARN radars.

    The plan was to use find_high_vel_overlap_events() to find event of interest and then to feed those events into
     this program.

    Contour/pixel comparison plots will have station1 along the x-axis.

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
            The starting epoch of the event to consider.
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
    :param plot_type: str (optional):
            The type of plot, either 'contour' or 'pixel', default is 'contour'.

    :return fig:  matplotlib.pyplot.figure:
            The figure. It can then be modified, added to, printed out, or saved in whatever format is desired.
    """

    Re = 6370  # Radius of the Earth, [km]
    seconds_in_an_hour = 3600
    minutes_in_an_hour = 60
    t_diffs = {station1: 0.000,
               station2: 0.000}
    param_ranges = {'v': (-1000, 1000),
                    'height': (0, 1000),
                    'adjElv': (0, 30)}
    param_n_bins = {'v': 50,
                    'height': 50
                    }
    count_min = 1  # Minimum number of points needed in a spatial/temporal 'clump' in order plot the matched point
    color_maps = {'v': 'seismic_r',  # For the range-time profiles
                  'adjElv': 'jet'}
    title_fontsize = 18
    time_interval_s = 60  # Controls spatial resolution of matched data

    if end_epoch < start_epoch:
        raise Exception("two_station_vel_comparison(): start epoch must be before end epoch.")

    # Compute parameter edges used in the contour/pixel comparisons
    param_edges = {'v': np.linspace(param_ranges['v'][0], param_ranges['v'][1], num=(param_n_bins['v'] + 1)),
                   'height': np.linspace(param_ranges['height'][0], param_ranges['height'][1],
                                         num=(param_n_bins['height'] + 1))}

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

    # Check input ranges, make sure they work for both radars
    freq_range = check_freq_range(freq_range=freq_range)
    gate_range = check_gate_range(gate_range=gate_range, hdw_info=first_radars_info.hardware_info)
    gate_range = check_gate_range(gate_range=gate_range, hdw_info=second_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range=beam_range, hdw_info=first_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range=beam_range, hdw_info=second_radars_info.hardware_info)

    # Make sure that both of the radars provided are in the same hemisphere
    reference_hemisphere = first_radars_info.hemisphere
    if second_radars_info.hemisphere != reference_hemisphere:
        raise Exception("two_station_vel_comparison(): " + station1 + " and " + station2
                        + " are not in the same hemisphere.  We can't compare them because they have no overlap.")

    starting_datetime = datetime.datetime.utcfromtimestamp(start_epoch)  # Note we are grabbing UTC datetime objs
    ending_datetime = datetime.datetime.utcfromtimestamp(end_epoch)

    start_hour = starting_datetime.hour + starting_datetime.minute / minutes_in_an_hour
    end_hour = ending_datetime.hour + ending_datetime.minute / minutes_in_an_hour
    hour_range = (start_hour, end_hour)  # TODO: Figure out if it is okay if this wraps through midnight
    # print("Starting hour: " + str(start_hour))
    # print("Ending hour: " + str(end_hour))

    year_range = (starting_datetime.year, ending_datetime.year)
    year_range = check_year_range(year_range=year_range)

    month_range = (starting_datetime.month, ending_datetime.month)
    month_range = check_month_range(month_range=month_range)

    day_range = (starting_datetime.day, ending_datetime.day)
    day_range = check_day_range(day_range=day_range)

    print("Getting SuperDARN data for " + station1.upper())
    df1 = get_data_handler(station=station1, year_range=year_range, month_range=month_range, day_range=day_range,
                           gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                           local_testing=local_testing, even_odd_days=None)
    df1 = only_keep_45km_res_data(df1)

    df1 = df1.loc[(df1['epoch'] >= start_epoch) & (df1['epoch'] <= end_epoch)]
    df1 = df1.loc[(df1['v'] > -1000) & (df1['v'] < 1000)]  # Remove extreme values
    df1.reset_index(drop=True, inplace=True)

    print("Restricting " + station1.upper() + "'s data - only keeping data for those cells that overlap with "
          + station2.upper())
    df1 = only_keep_overlap(station=station1, df=df1, other_station=station2, gate_min=gate_range[0])

    print("Getting SuperDARN data for " + station2.upper())
    df2 = get_data_handler(station=station2, year_range=year_range, month_range=month_range, day_range=day_range,
                           gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                           local_testing=local_testing, even_odd_days=None)
    df2 = only_keep_45km_res_data(df2)

    df2 = df2.loc[(df2['epoch'] >= start_epoch) & (df2['epoch'] <= end_epoch)]
    df2 = df2.loc[(df2['v'] > -1000) & (df2['v'] < 1000)]  # Remove extreme values
    df2.reset_index(drop=True, inplace=True)

    print("Restricting " + station2.upper() + "'s data - only keep data for those cells that overlap with "
          + station1.upper())
    df2 = only_keep_overlap(station=station2, df=df2, other_station=station1, gate_min=gate_range[0])

    print("Recomputing Elevation Angles for " + station1.upper() + " data - t_diff=" + str(t_diffs[station1]))
    elevation_v2(station=station1, df=df1, t_diff=t_diffs[station1])  # t_diff is in microseconds

    print("Recomputing Elevation Angles for " + station2.upper() + " data - t_diff=" + str(t_diffs[station2]))
    elevation_v2(station=station2, df=df2, t_diff=t_diffs[station2])  # t_diff is in microseconds

    # print(df1[['slist', 'adjElv']])
    # print(df2[['slist', 'adjElv']])

    print("Computing virtual heights for " + station1.upper())
    slant_range1 = 120 + 45 / 2 + 45 * df1['slist']  # Slant range [km]
    df1['height'] = np.sqrt(Re * Re + slant_range1 * slant_range1 + 2 * Re * slant_range1
                            * np.sin(np.radians(np.asarray(df1['adjElv'])))) - Re

    print("Computing virtual heights for " + station2.upper())
    slant_range2 = 120 + 45 / 2 + 45 * df2['slist']  # Slant range [km]
    df2['height'] = np.sqrt(Re * Re + slant_range2 * slant_range2 + 2 * Re * slant_range2
                            * np.sin(np.radians(np.asarray(df2['adjElv'])))) - Re

    print("Adding decimal time to " + station1.upper() + "'s data.")
    station1_stid = pydarn.read_hdw_file(station1).stid
    df1 = add_decimal_hour_to_df(df=df1, time_units='ut', stid=station1_stid,
                                 date_time_est=(df1['datetime'].iat[0]).to_pydatetime())

    print("Adding decimal time to " + station2.upper() + "'s data.")
    station2_stid = pydarn.read_hdw_file(station2).stid
    df2 = add_decimal_hour_to_df(df=df2, time_units='ut', stid=station2_stid,
                                 date_time_est=(df2['datetime'].iat[0]).to_pydatetime())

    print("Preparing the figure...")
    fig = plt.figure(figsize=[23, 10], constrained_layout=True, dpi=300)
    axes = add_axes(fig=fig, reference_hemisphere=reference_hemisphere)

    apply_subplot_formatting(axes=axes, station1=station1, station2=station2,
                             reference_hemisphere=reference_hemisphere, t_diffs=t_diffs,
                             gate_range=gate_range, vel_range=param_ranges['v'], height_range=param_ranges['height'],
                             hour_range=hour_range)

    title_figure(fig=fig, station1=station1, station2=station2,
                 starting_datetime=starting_datetime, ending_datetime=ending_datetime,
                 beam_range=beam_range, gate_range=gate_range, freq_range=freq_range,
                 title_fontsize=title_fontsize)

    print("Working on the range-time profile for " + station1.upper())
    complete_range_time_profile(fig=fig, axes=axes['station1_rti'], df=df1, gate_range=gate_range,
                                hour_range=hour_range, param_ranges=param_ranges, color_maps=color_maps)

    print("Working on the range-time profile for " + station2.upper())
    complete_range_time_profile(fig=fig, axes=axes['station2_rti'], df=df2, gate_range=gate_range,
                                hour_range=hour_range, param_ranges=param_ranges, color_maps=color_maps)

    print("Working on the map")
    complete_map_subplot(ax=axes['map'], station1=station1, station2=station2,
                         gate_range=gate_range, beam_range=beam_range)

    print("Building a matched dataframe...")
    matched_df = build_two_radar_matched_data(station1=station1, df1=df1, station2=station2, df2=df2,
                                              time_interval_s=time_interval_s, gate_min=gate_range[0])
    matched_df = matched_df.loc[(matched_df['count1'] >= count_min) & (matched_df['count2'] >= count_min)]

    print("Completing the large comparison plots")
    complete_comparison_plot(fig=fig, axes=axes['full_event_comparison'], matched_df=matched_df,
                             param_edges=param_edges, cbar=True, plot_type=plot_type)

    # Divide the event into subsections, we will produce small comparison plots for each sub-section
    number_of_sub_sections = 4
    sub_event_edges = np.linspace(start=start_epoch, stop=end_epoch, num=(number_of_sub_sections + 1), dtype=int)
    sub_event_duration = sub_event_edges[1] - sub_event_edges[0]

    for i, starting_sub_event_edge in enumerate(sub_event_edges):
        if starting_sub_event_edge == sub_event_edges[-1]:
            continue  # The last edge is not a slice start

        ending_sub_event_edge = starting_sub_event_edge + sub_event_duration

        matched_df_sub_event = matched_df[(matched_df['matched_epoch'] >= starting_sub_event_edge) &
                                          (matched_df['matched_epoch'] <= ending_sub_event_edge)]

        print("Completing the comparison plots for interval " + str(i))
        complete_comparison_plot(fig=fig, axes=axes['interval' + str(i)], matched_df=matched_df_sub_event,
                                 param_edges=param_edges, cbar=False, plot_type=plot_type)

        # We also need to title each of the small plots with the sub-interval times
        for param, ax in axes['interval' + str(i)].items():
            sub_event_starting_datetime = datetime.datetime.utcfromtimestamp(starting_sub_event_edge)
            sub_event_ending_datetime = datetime.datetime.utcfromtimestamp(ending_sub_event_edge)

            ax.set_title(str(sub_event_starting_datetime.time())[:-3] + " to "
                         + str(sub_event_ending_datetime.time())[:-3] + " (UTC)")

    return fig


def complete_comparison_plot(fig, axes, matched_df, param_edges, cbar=False, plot_type='contour'):
    """

    Complete the comparison plot subplots.

    :param fig: matplotlib.pyplot.figure:
            The figure to draw on.
    :param axes: dictionary of matplotlib.pyplot.axis:
            The axes to draw on.  Dictionary keys are the parameters to plot.
    :param matched_df: pandas.DataFrame:
            A dataframe with one row for each temporal/spatial chunk
            As returned from build_two_radar_matched_data()
    :param param_edges: dictionary of numpy.array objects:
            Dictionary containing an array of edges for each parameter
    :param cbar: bool (optional; default is False)
            Whether or not to add a colour bar to plot
    :param plot_type: (optional; default is 'contour')
            The type of comparison plot wanted: either pixel or contour
    """

    if len(matched_df) <= 0:
        return

    cmap = 'jet'

    for param, ax in axes.items():
        result, _, _, _ = stats.binned_statistic_2d(matched_df[param + '1'], matched_df[param + '2'], values=None,
                                                    statistic='count', bins=[param_edges[param], param_edges[param]])
        if plot_type == 'pixel':
            plot = ax.pcolormesh(param_edges[param], param_edges[param], result.transpose(), cmap=cmap, zorder=0)

        elif plot_type == 'contour':

            bin_width = (param_edges[param][1] - param_edges[param][0])
            bin_centers = param_edges[param][1:] - bin_width / 2

            plot = ax.contourf(bin_centers, bin_centers, result.transpose(), cmap=cmap, zorder=0)

        else:
            raise Exception("complete_comparison_plot(): plot_type " + str(plot_type) + " not recognized.")

        if cbar:
            cbar = fig.colorbar(plot, ax=ax, format='%.1f', orientation="vertical")


def complete_range_time_profile(fig, axes, df, gate_range, hour_range, param_ranges, color_maps):
    """

    Complete the range time profile subplots.

    :param fig: matplotlib.pyplot.figure:
            The figure to draw on.
    :param axes: dictionary of matplotlib.pyplot.axis:
            The axes to draw on.  Dictionary keys are the parameters to plot.
    :param df: pandas.DataFrame:
            A SuperDARN fit dataframe.  Must contain parameters for each of the axes' keys.
    :param gate_range: (int, int):
            The gate range to profile.
    :param hour_range: (float, float):
            The decimal-hour range to profile in UT.
    :param color_maps: dictionary of matplotlib.colors.Colormap:
            The colour maps to use.  Must be an entry for each of the axes' keys.
    """

    cbar_text_format = '%d'

    # Compute hour edges
    bins_per_hour = 60
    n_bins_x = int((hour_range[1] - hour_range[0]) * bins_per_hour)
    hour_edges = np.linspace(hour_range[0], hour_range[1], num=(n_bins_x + 1))

    # Compute gate_edges
    bins_per_gate = 1
    n_bins_y = int(((gate_range[1] + 1) - gate_range[0]) * bins_per_gate)
    gate_edges = np.linspace(gate_range[0], gate_range[1] + 1, num=(n_bins_y + 1), dtype=int)

    for param, ax in axes.items():

        if len(df) <= 0:
            continue

        else:
            result, _, _, _ = stats.binned_statistic_2d(df['ut'], df['slist'], values=df[param],
                                                        bins=[hour_edges, gate_edges])

            plot = ax.pcolormesh(hour_edges, gate_edges, result.transpose(),
                                 cmap=color_maps[param], vmin=param_ranges[param][0], vmax=param_ranges[param][1],
                                 zorder=0)

            if param == 'v':
                cbar = fig.colorbar(plot, ax=ax, orientation="vertical", format=cbar_text_format, extend='both')
            else:
                cbar = fig.colorbar(plot, ax=ax, orientation="vertical", format=cbar_text_format, extend='max')

        # ax.grid(b=True, which='major', axis='both', linewidth=1, linestyle='-')
        # ax.grid(b=True, which='minor', axis='both', linewidth=0.4, linestyle='--')


def title_figure(fig, station1, station2, starting_datetime, ending_datetime,
                 beam_range, gate_range, freq_range,
                 title_fontsize):
    """
    :param fig: The figure to title
    :param station1: str:
            The first radar station to consider, as 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param station2: str:
            The second radar station to consider, again as a 3 character string.
    :param starting_datetime: datetime.datetime:
            The start of the event
    :param ending_datetime:
            The end of the event
    :param beam_range: See two_station_vel_comparison() docstring
    :param gate_range: See two_station_vel_comparison() docstring
    :param freq_range: See two_station_vel_comparison() docstring
    :param title_fontsize:
    """

    beam_string = "Allowed Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    gate_string = "Allowed Gates " + str(gate_range[0]) + "-" + str(gate_range[1])
    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"

    fig.suptitle(station1.upper() + " and " + station2.upper() + " Velocity Comparison.  From "
                 + str(starting_datetime) + " to " + str(ending_datetime) + " (UTC)."
                                                                            "\n" + gate_string + "; " + beam_string + "; " + freq_string +
                 "\nProduced by " + str(os.path.basename(__file__)), fontsize=title_fontsize)


def complete_map_subplot(ax, station1, station2, gate_range, beam_range, colours=None):
    """

    Complete the map subplot.
    It has both fans, with the area of overlap highlighted.

    :param ax: matplotlib.axes:
            The axes to draw on
    :param station1: str:
            The first radar station to consider, as 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param station2: str:
            The second radar station to consider, again as a 3 character string.
    :param gate_range: See two_station_vel_comparison() docstring
    :param beam_range: See two_station_vel_comparison() docstring

    :param colours: dictionary of named colours: (optional; default is {station1: "red", station2: "blue"}):
            The colour to use for each station's dot/shading
    """

    if colours == None:
        colours = {station1: "red",
                   station2: "blue"}

    for station in [station1, station2]:

        """ Plot both radar's fans"""
        all_radars_info = SuperDARNRadars()
        station_stid = pydarn.read_hdw_file(station).stid
        radar_info = all_radars_info.radars[station_stid]
        radar_lon = radar_info.hardware_info.geographic.lon
        radar_lat = radar_info.hardware_info.geographic.lat

        print("Plotting the fan for " + station)
        print("     Radar lat: " + str(radar_lat))
        print("     Radar lon: " + str(radar_lon))

        # Plot the radar as a dot
        ax.plot([radar_lon, radar_lon], [radar_lat, radar_lat], marker="o", markersize=3, transform=ccrs.Geodetic(),
                label=station, color=colours[station])

        cell_corners_lats, cell_corners_lons = radar_fov(stid=station_stid, coords='geo')

        # plot all the beam boundary lines
        for beam_line in range(beam_range[0], beam_range[1] + 2):
            ax.plot(cell_corners_lons[gate_range[0]:gate_range[1] + 2, beam_line],
                    cell_corners_lats[gate_range[0]:gate_range[1] + 2, beam_line],
                    color='black', linewidth=0.1, transform=ccrs.Geodetic(), zorder=4)

        # plot the arcs boundary lines
        for range_ in range(gate_range[0], gate_range[1] + 2):
            ax.plot(cell_corners_lons[range_, beam_range[0]:beam_range[1] + 2],
                    cell_corners_lats[range_, beam_range[0]:beam_range[1] + 2],
                    color='black', linewidth=0.1, transform=ccrs.Geodetic(), zorder=4)

        """ Shade the overlapping areas"""
        if station == station1:
            other_station = station2
        else:
            other_station = station1
        # Build reduced arrays containing only the cells in the specified gate/beam range
        reduced_cell_corners_lons = cell_corners_lons[gate_range[0]: gate_range[1] + 2,
                                    beam_range[0]: beam_range[1] + 2]
        reduced_cell_corners_lats = cell_corners_lats[gate_range[0]: gate_range[1] + 2,
                                    beam_range[0]: beam_range[1] + 2]

        # Loop through all valid gates and beams - mark valid overlaps
        num_gates = (gate_range[1] + 1) - gate_range[0]
        num_beams = (beam_range[1] + 1) - beam_range[0]
        scan = np.zeros(shape=(num_gates, num_beams))
        scan[:] = np.nan

        for gate_idx in range(num_gates):
            for beam_idx in range(num_beams):
                gate = gate_range[0] + gate_idx
                beam = beam_range[0] + beam_idx
                # print("Gate: " + str(gate) + ", beam : " + str(beam))
                station2_beam, station2_gate = compute_df_radar_overlap(station1=station, station1_beam=beam,
                                                                        station1_gate=gate, station2=other_station,
                                                                        gate_min=gate_range[0])

                if station2_beam is not None and station2_gate is not None:
                    # Then we have a valid overlap
                    scan[gate_idx, beam_idx] = 1

        cmap = single_colour_cmap(named_color=colours[station])
        ax.pcolormesh(reduced_cell_corners_lons, reduced_cell_corners_lats, scan,
                      transform=ccrs.PlateCarree(), cmap=cmap, alpha=0.2, zorder=3)

    ax.legend(loc='upper right', ncol=2)


def apply_subplot_formatting(axes, station1, station2, reference_hemisphere, hour_range, t_diffs=None,
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
    :param hour_range: (float, float):
            The hour range to consider.

    :param t_diffs: Dictionary of floats: (Optional, default is {station1: 0.000, station2: 0.000})
            The extra time delay added in when computing adjusted elevation angle, in microseconds.
    :param gate_range: (int, int) (Optional; default is (0, 74))
            See two_station_vel_comparison() docstring
    :param height_range: (float, float) (Optional; default is (0, 300))
            The height range in km
    :param vel_range: (float, float) (Optional; default is (-600, 600))
            The velocity range in m/s


    """

    print("Using a vel range of " + str(vel_range))

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
                x_lim = hour_range

                ax.set_ylim(y_lim)
                ax.yaxis.set_major_locator(MultipleLocator(10))  # Every 10 gates
                ax.yaxis.set_minor_locator(MultipleLocator(1))

                ax.set_xlim(x_lim)
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

                if subplot_type == "v":
                    ax.set_title(station_mnemonic.upper() + " Velocities", fontsize=title_font_size)
                elif subplot_type == "adjElv":
                    ax.set_title(station_mnemonic.upper() + " Adjusted Elevation Angles, tdiff=" + str(t_diff)
                                 + " \u03BCs", fontsize=title_font_size)
                    ax.set_xlabel("UT Time [hour]", fontsize=label_font_size)

        elif axis_key == 'map':
            # We are formatting the large map
            lat_extreme = reference_hemisphere.value * 40  # deg
            ax = axis_item
            ax.set_extent([-180, 180, reference_hemisphere.value * 90, lat_extreme], crs=ccrs.PlateCarree())
            ax.gridlines()
            # ax.add_feature(cfeature.OCEAN)
            ax.add_feature(cfeature.LAND)

        else:
            # We are formatting the comparison plots

            for subplot_type, ax in axis_item.items():

                if subplot_type == "v":
                    ax.set_ylim(vel_range)
                    ax.set_xlim(vel_range)

                    ax.set_xlabel(station1.upper() + " Velocities [m/s]")
                    ax.set_ylabel(station2.upper() + " Velocities [m/s]")

                    ax.yaxis.set_major_locator(MultipleLocator(500))
                    ax.xaxis.set_major_locator(MultipleLocator(500))
                    ax.yaxis.set_minor_locator(MultipleLocator(100))
                    ax.xaxis.set_minor_locator(MultipleLocator(100))
                    if axis_key == "full_event_comparison":
                        ax.set_title("Full Event Velocity Comparison", fontsize=title_font_size)
                        # plt.sca(ax)
                        # plt.xticks([-600, -400, -200, 0, 200, 400, 600])
                        # plt.yticks([-600, -400, -200, 0, 200, 400, 600])

                elif subplot_type == "height":
                    ax.set_ylim(height_range)
                    ax.set_xlim(height_range)

                    ax.set_xlabel(station1.upper() + " Virtual Heights [km]")
                    ax.set_ylabel(station2.upper() + " Virtual Heights [km]")

                    ax.yaxis.set_major_locator(MultipleLocator(500))
                    ax.xaxis.set_major_locator(MultipleLocator(500))
                    ax.yaxis.set_minor_locator(MultipleLocator(100))
                    ax.xaxis.set_minor_locator(MultipleLocator(100))
                    if axis_key == "full_event_comparison":
                        ax.set_title("Full Event Virtual Height Comparison", fontsize=title_font_size)
                        # plt.sca(ax)
                        # plt.xticks([100, 150, 200, 250, 300])
                        # plt.yticks([100, 150, 200, 250, 300])

                ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
                ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
                ax.plot(ax.get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
                ax.plot([0, 0], ax.get_xlim(), linestyle='-', linewidth=0.5, color='black')
                ax.plot([ax.get_ylim()[0], ax.get_ylim()[1]], [ax.get_xlim()[0], ax.get_xlim()[1]],
                        linestyle='--', linewidth=2, color=bisector_colour)


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
    axes['station1_rti'] = {"v": fig.add_subplot(gs[0, 0:3]), "adjElv": fig.add_subplot(gs[1, 0:3])}
    axes['station2_rti'] = {"v": fig.add_subplot(gs[0, 6:9]), "adjElv": fig.add_subplot(gs[1, 6:9])}

    # We will have a map that shows each radar's fan and highlights the overlap
    if reference_hemisphere.value == 1:
        proj = ccrs.NorthPolarStereo()
    elif reference_hemisphere.value == -1:
        proj = ccrs.SouthPolarStereo()
    else:
        raise Exception("Hemisphere not recognized.")
    axes['map'] = fig.add_subplot(gs[0:2, 3:6], projection=proj)

    # We will have a couple larger plots that show select comparisons for the whole event duration
    axes['full_event_comparison'] = {"v": fig.add_subplot(gs[2:4, 0:3]), "height": fig.add_subplot(gs[2:4, 3:6])}

    # We will have some small plots that show select comparisons for parts of the event
    axes['interval0'] = {"v": fig.add_subplot(gs[0, 9]), "height": fig.add_subplot(gs[0, 10])}
    axes['interval1'] = {"v": fig.add_subplot(gs[1, 9]), "height": fig.add_subplot(gs[1, 10])}
    axes['interval2'] = {"v": fig.add_subplot(gs[2, 9]), "height": fig.add_subplot(gs[2, 10])}
    axes['interval3'] = {"v": fig.add_subplot(gs[3, 9]), "height": fig.add_subplot(gs[3, 10])}

    return axes


if __name__ == '__main__':
    """ Testing """

    local_testing = False

    if local_testing:
        station1 = "dcn"
        station2 = "mcm"

        start_epoch = 1321070400  # Saturday, November 12, 2011 4:00:00 AM UTC = 1321070400
        end_epoch = 1321084800  # Saturday, November 12, 2011 8:00:00 AM UTC = 1321084800

        gate_range = (20, 74)
        beam_range = (0, 15)

        # Note: year, month, and day don't matter for local testing
        fig = two_station_vel_adjElv_comparison(station1=station1, station2=station2,
                                                start_epoch=start_epoch, end_epoch=end_epoch,
                                                gate_range=gate_range, beam_range=beam_range, freq_range=None,
                                                local_testing=local_testing, plot_type='contour')

        plt.show()


    else:

        event_df = get_select_dcn_mcm_overlap_events()

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        for i in range(len(event_df)):
            station1 = event_df['station1'].iat[i]
            station2 = event_df['station2'].iat[i]
            start_epoch = event_df['start_epoch'].iat[i]
            end_epoch = event_df['end_epoch'].iat[i]
            gate_range = event_df['gate_range'].iat[i]
            beam_range = event_df['beam_range'].iat[i]
            freq_range = event_df['freq_range'].iat[i]
            plot_type = event_df['plot_type'].iat[i]

            print("Running from " + str(start_epoch) + " to " + str(end_epoch) + " for "
                  + station1.upper() + " and " + station2.upper())

            fig = two_station_vel_adjElv_comparison(station1=station1, station2=station2,
                                                    start_epoch=start_epoch, end_epoch=end_epoch,
                                                    gate_range=gate_range, beam_range=beam_range, freq_range=None,
                                                    local_testing=local_testing, plot_type=plot_type)

            out_fig = out_dir + "/two_station_comparison_vel_adjElv-" + station1 + "_" + station2 + \
                      "-from_" + str(start_epoch) + "_to_" + str(end_epoch)

            print("Saving plot as " + out_fig)
            fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
