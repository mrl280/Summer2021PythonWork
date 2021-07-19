import calendar
import math
import os
import pathlib
import pydarn

from aacgmv2 import get_aacgm_coord
from cartopy.util import add_cyclic_point
from matplotlib import pyplot as plt
from pydarn import radar_fov, SuperDARNRadars

import cartopy.crs as ccrs
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import numpy as np

from lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from lib.cm.modified_jet import modified_jet
from lib.add_mlt_to_df import add_mlt_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def occ_clock_diagram(station, year, month_range=None, day_range=None, gate_range=None, beam_range=None,
                      freq_range=None, time_units='mlt', plot_type='contour', local_testing=False):
    """

    Produce a full circle stereographic occurrence plot (12 at the top).

    # TODO: Figure out how to fill in the center of the circle Cyclic point is not working (adds in a 95 deg point),
       and setting the last point to 90 deg is not working

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only considers 45 km data.
            (a warning will be printed if other spatial resolution data is stripped from the dataset)
        - To check which fitACF program is being used, refer to the data readers in lib.data_getters

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year: int:
            The year to consider.
    :param month_range: (<int>, <int>) (optional):
            Inclusive. The months of the year to consider.  If omitted (or None), then all days will be considered.
    :param day_range: (<int>, <int>) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param gate_range: (<int>, <int>) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param freq_range: (<float>, <float>) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.
    :param time_units: str:  # TODO: Not sure if this type of plot makes sense for time units other than mlt?
            The time units to plot on the circle, 12 is always at the top.  Default is 'mlt'.
                'ut' for universal time
                'mlt' for magnetic local time
                'lt' for local time (based on longitude)
                'lst' for local standard time (based on time zones)
                'ast' for apparent solar time (based on the apparent angular motion of the sun across the sky)
    :param plot_type: str (optional):
            The type of plot, either 'contour' or 'pixel', default is 'contour'
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :return: matplotlib.pyplot.figure: The figure.
            It can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    is_vmax = 0.8
    gs_vmax = 0.4

    time_units = check_time_units(time_units)
    year = check_year(year)
    month_range = check_month_range(month_range)
    freq_range = check_freq_range(freq_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid
    radar_lon = this_radars_info.hardware_info.geographic.lon
    radar_lat = this_radars_info.hardware_info.geographic.lat

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    gate_string = "Gates " + str(gate_range[0]) + "-" + str(gate_range[1])
    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"

    print("     Retrieving data...")
    # We will force build a new dataframe because the pickled one probably won't be for all beams/gates
    df = get_data_handler(station, year_range=(year, year), month_range=month_range, day_range=day_range,
                          gate_range=gate_range, beam_range=beam_range, freq_range=freq_range, occ_data=True,
                          force_build=True, local_testing=local_testing)
    df = only_keep_45km_res_data(df)

    if month_range[0] == month_range[1]:
        month_string = calendar.month_name[month_range[0]]
        date_time_est, _ = build_datetime_epoch(year=year, month=month_range[0], day=15, hour=0)
    else:
        month_string = calendar.month_name[month_range[0]] + " to " + calendar.month_name[month_range[1]]
        # Use the middle of the year as an estimate
        date_time_est, _ = build_datetime_epoch(year=year, month=6, day=15, hour=0)

    print("     Computing MLTs for " + str(year) + " data...")
    cell_corners_aacgm_lats, cell_corners_aacgm_lons = radar_fov(stid=radar_id, coords='aacgm', date=date_time_est)

    df = add_mlt_to_df(cell_corners_aacgm_lons=cell_corners_aacgm_lons, cell_corners_aacgm_lats=cell_corners_aacgm_lats,
                       df=df)

    if time_units != 'mlt':
        # We already have mlt, otherwise add the required time units
        df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id, date_time_est=date_time_est)

    # Right now xdata is in the range 0-24, we need to put it in the range 0-360 for circular plotting
    df['xdata'] = 15 * df[time_units]

    # Compute a circle in axis coordinates which can be used as a boundary
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5  # 0.5 means middle of the circle
    vertices = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(vertices * radius + center)

    # Convert radar coordinates to aacgm
    radar_lat_aacgm, radar_lon_aacgm, radar_mlt = get_aacgm_coord(radar_lat, radar_lon, 0, date_time_est)
    radar_mlts = np.arange(0, 360, 1)
    radar_lats_aacgm = np.asarray([radar_lat_aacgm] * len(radar_mlts))

    # Regardless of hemisphere, we will always use the North Pole projection because it gives us a stereo projection
    # with the zero degree line at the bottom of the plot and East is CCW (this is what we need for a clock diagram).
    projection = ccrs.NorthPolarStereo()

    # Since we are using the North Pole projection we need all latitudes to be positive
    df['lat'] = df['lat'].abs()

    print("     Preparing the plot...")
    fig, ax = plt.subplots(figsize=[10, 6], dpi=300, nrows=1, ncols=2, constrained_layout=True,
                           subplot_kw={'projection': projection})
    fig.suptitle(month_string + " " + str(year) + " at " + station.upper() + "; " + beam_string + "; " + freq_string +
                 "\n" + gate_string + "; Data Plotted in AACGM Latitudes and " + time_units.upper() +
                 "\nProduced by " + str(os.path.basename(__file__)), fontsize=18)

    # Apply common subplot formatting
    lat_extreme = 70  # TODO: Adjust lat extreme based on the radar (use the most extreme point in the fan)
    for i in range(ax.size):
        ax[i].set_extent([-180, 180, 90, lat_extreme], crs=ccrs.PlateCarree())
        ax[i].set_boundary(circle, transform=ax[i].transAxes)  # Comment out this line to enable auto labelling

        # Add gridlines
        gl = ax[i].gridlines(draw_labels=True, linestyle='--', linewidth=0.75, zorder=4)
        gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        # Latitude labels won't draw on a circular plot, add them on manually
        add_lat_labels_to_clock_diagram(ax=ax[i], hemisphere=this_radars_info.hemisphere,
                                        lat_extreme=lat_extreme, lat_labels=[85, 80, 75])

        # Print clock numbers
        text_offset_multiplier = 1.03
        ax[i].text(0, text_offset_multiplier * ax[i].get_ylim()[1], "12", ha='center', va='bottom', fontsize=14)
        ax[i].text(0, text_offset_multiplier * ax[i].get_ylim()[0], "00", ha='center', va='top', fontsize=14)
        ax[i].text(text_offset_multiplier * ax[i].get_xlim()[1], 0, "06", ha='left', va='center', fontsize=14)
        ax[i].text(text_offset_multiplier * ax[i].get_xlim()[0], 0, "18", ha='right', va='center', fontsize=14)

        # Print out the time units
        ax[i].text(ax[i].get_xlim()[1], ax[i].get_ylim()[1], time_units.upper(), ha='right', va='top', fontsize=14)

        # Plot radar track
        ax[i].plot(radar_mlts, radar_lats_aacgm, color='k', linewidth=1, linestyle="--", transform=ccrs.Geodetic())

    # Print out echo types
    ax[0].text(ax[0].get_xlim()[0], ax[0].get_ylim()[1], "IS", ha='left', va='top', fontsize=14)
    ax[1].text(ax[1].get_xlim()[0], ax[1].get_ylim()[1], "GS", ha='left', va='top', fontsize=14)

    print("     Computing binned occ rates...")
    # Compute mlt edges
    deg_mlt_per_bin = 2
    n_bins_mlt = int(360 / deg_mlt_per_bin)
    mlt_edges = np.linspace(0, 360, num=(n_bins_mlt + 1))
    delta_mlt = mlt_edges[1] - mlt_edges[0]

    # Compute latitude edges
    n_bins_lat = 90 - lat_extreme  # One bin per degree of latitude
    lat_edges = np.linspace(lat_extreme, 90, num=(n_bins_lat + 1))
    delta_lat = lat_edges[1] - lat_edges[0]

    contour_data_is = np.empty(shape=(n_bins_mlt, n_bins_lat))
    contour_data_is[:] = math.nan

    contour_data_gs = np.empty(shape=(n_bins_mlt, n_bins_lat))
    contour_data_gs[:] = math.nan

    for mlt_idx, start_mlt in enumerate(mlt_edges):
        if start_mlt == mlt_edges[-1]:
            continue  # The last edge is not a start
        end_mlt = start_mlt + delta_mlt
        df_mlt = df[(df['xdata'] >= start_mlt) & (df['xdata'] <= end_mlt)]

        for lat_idx, start_lat in enumerate(lat_edges):
            if start_lat == lat_edges[-1]:
                continue  # The last edge is not a start

            end_lat = start_lat + delta_lat
            df_mlt_lat = df_mlt[(df_mlt['lat'] >= start_lat) & (df_mlt['lat'] <= end_lat)]

            try:
                contour_data_is[mlt_idx][lat_idx] = sum(df_mlt_lat['good_iono_echo']) / len(df_mlt_lat)
                contour_data_gs[mlt_idx][lat_idx] = sum(df_mlt_lat['good_grndscat_echo']) / len(df_mlt_lat)
            except ZeroDivisionError:
                # There are no point in this interval
                contour_data_is[mlt_idx][lat_idx] = math.nan
                contour_data_gs[mlt_idx][lat_idx] = math.nan
            except BaseException as e:
                print("MLT index: " + str(mlt_idx))
                print("LAT index: " + str(lat_idx))
                raise e

    if plot_type == "contour":

        # Compute bin centers
        bin_xwidth = (mlt_edges[1] - mlt_edges[0])
        bin_ywidth = (lat_edges[1] - lat_edges[0])
        bin_xcenters = mlt_edges[1:] - bin_xwidth / 2
        bin_ycenters = lat_edges[1:] - bin_ywidth / 2

        # Adding a cyclic points is required to complete the circle
        # Without adding in this point, there will be one pie-shaped piece missing from the circle
        contour_data_is_cyclic, bin_xcenters_cyclic = add_cyclic_point(contour_data_is.transpose(), coord=bin_xcenters)
        contour_data_gs_cyclic, bin_xcenters_cyclic = add_cyclic_point(contour_data_gs.transpose(), coord=bin_xcenters)

        n_levels = 12
        is_levels = np.linspace(start=0, stop=is_vmax, num=(n_levels + 1))
        gs_levels = np.linspace(start=0, stop=gs_vmax, num=(n_levels + 1))
        cmap = 'jet'
        # cmap = modified_jet(levels=len(levels) - 1)

        if is_vmax < 1:
            plot0 = ax[0].contourf(bin_xcenters_cyclic, bin_ycenters, contour_data_is_cyclic,
                                   cmap=cmap, levels=is_levels, transform=ccrs.PlateCarree(), extend='max')
        else:
            plot0 = ax[0].contourf(bin_xcenters_cyclic, bin_ycenters, contour_data_is_cyclic,
                                   cmap=cmap, levels=is_levels, transform=ccrs.PlateCarree())

        if gs_vmax < 1:
            plot1 = ax[1].contourf(bin_xcenters_cyclic, bin_ycenters, contour_data_gs_cyclic,
                                   cmap=cmap, levels=gs_levels, transform=ccrs.PlateCarree(), extend='max')
        else:
            plot1 = ax[1].contourf(bin_xcenters_cyclic, bin_ycenters, contour_data_gs_cyclic,
                                  cmap=cmap, levels=gs_levels, transform=ccrs.PlateCarree())

    elif plot_type == "pixel":

        cmap = 'jet'
        plot0 = ax[0].pcolormesh(mlt_edges, lat_edges, contour_data_is.transpose(),
                                 transform=ccrs.PlateCarree(), cmap=cmap, vmin=0, vmax=is_vmax)
        plot1 = ax[1].pcolormesh(mlt_edges, lat_edges, contour_data_gs.transpose(),
                                 transform=ccrs.PlateCarree(), cmap=cmap, vmin=0, vmax=gs_vmax)

    else:
        raise Exception("plot_type not recognized")

    if is_vmax < 1:
        cbar0 = fig.colorbar(plot0, ax=ax[0], shrink=0.7, orientation="horizontal", format='%.1f', extend='max')
    else:
        cbar0 = fig.colorbar(plot0, ax=ax[0], shrink=0.7, orientation="horizontal", format='%.1f')
    cbar0.ax.tick_params(labelsize=16, labelrotation=0)

    if gs_vmax < 1:
        cbar1 = fig.colorbar(plot1, ax=ax[1], shrink=0.7, orientation="horizontal", format='%.1f', extend='max')
    else:
        cbar1 = fig.colorbar(plot1, ax=ax[1], shrink=0.7, orientation="horizontal", format='%.1f')
    cbar1.ax.tick_params(labelsize=16, labelrotation=0)

    return df, fig


def add_lat_labels_to_clock_diagram(ax, hemisphere, lat_extreme, lat_labels, color='white', fontsize=12):
    """

    Clock diagram axis are clipped using a matplotlib path object and labels wont draw on these circular axis,
     so this function prints them our manually.
    Also, we always NorthPoleStereo projection and this allows us to print the proper hemisphere (N or S)

    :param ax: matplotlib.pyplot.axis:
            The axis to draw on
    :param hemisphere: A pydarn hemisphere object:
            Which hemisphere is the radar in?
    :param lat_extreme:
    :param lat_labels: list/array of ints:
            The latitude labels to print out
    :param color: Matplotlib named colour as a string: (optional)
            The colour for the labels, default is 'white'
    :param fontsize: int: (optional)
            The fontsize, default is 12
    """

    degree_sign = u'\N{DEGREE SIGN}'
    angle_to_print_on = 135  # deg clockwise starting at the top

    if hemisphere.value == 1:
        hemi_str = "N"
    elif hemisphere.value == -1:
        hemi_str = "S"
    else:
        raise Exception("Error in add_lat_labels_to_clock_diagram(): Hemisphere not recognized")

    lat_range = 90 - lat_extreme  # The range of latitudes that are covered by the plot

    for lat in lat_labels:
        distance_from_edge = lat - lat_extreme
        percent_of_way_to_edge = 1 - distance_from_edge / lat_range

        x_scale = percent_of_way_to_edge * math.sin(math.radians(angle_to_print_on))
        y_scale = percent_of_way_to_edge * math.cos(math.radians(angle_to_print_on))

        ax.text(x_scale * ax.get_xlim()[1], y_scale * ax.get_ylim()[1], str(lat) + degree_sign + hemi_str,
                color=color, fontsize=fontsize, ha='center', va='center', zorder=5)


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    if local_testing:
        station = "dce"

        _, fig = occ_clock_diagram(station=station, year=2011, month_range=(11, 11), day_range=None,
                                   gate_range=(0, 74), beam_range=(6, 7), freq_range=None,
                                   plot_type='pixel', time_units='mlt',
                                   local_testing=local_testing)

        plt.show()


    else:
        station = "dce"
        year = 2019
        month = 3
        freq_range = (9.5, 12.5)
        plot_type = 'pixel'

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        _, fig = occ_clock_diagram(station=station, year=year, month_range=(month, month), day_range=(14, 28),
                                   gate_range=(0, 74), beam_range=None, freq_range=freq_range,
                                   plot_type=plot_type, time_units='mlt',
                                   local_testing=local_testing)

        out_fig = out_dir + "/occ_clock_diagram_" + station + "-" + str(year) + "-" + \
                  str(month) + "_" + str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz" + "_" + plot_type
        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
