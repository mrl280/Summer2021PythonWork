import calendar
import math
import pathlib
import pydarn

from aacgmv2 import get_aacgm_coord, convert_latlon_arr
from cartopy.util import add_cyclic_point
from matplotlib import pyplot as plt
from pydarn import radar_fov, SuperDARNRadars

import cartopy.crs as ccrs
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import numpy as np

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

    Produce a full circle stereographic occurrence plot in either ut or mlt (12 at the top).

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
    :param time_units: str: 'ut' for universal time or 'mlt' for magnetic local time:
            The time units to plot on the circle, 12 is always at the top.  Default is 'mlt'
    :param plot_type: str (optional):
            The type of plot, either 'contour' or 'pixel', default is 'contour'
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :return: matplotlib.pyplot.figure: The figure.
            It can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    time_units = check_time_units(time_units)
    year = check_year(year)
    month_range = check_month_range(month_range)
    freq_range = check_freq_range(freq_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid
    hemisphere = this_radars_info.hemisphere
    radar_lon = this_radars_info.hardware_info.geographic.lon
    radar_lat = this_radars_info.hardware_info.geographic.lat

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    print("     Retrieving data...")
    df = get_data_handler(station, year_range=(year, year), month_range=month_range, day_range=day_range,
                          gate_range=gate_range, beam_range=beam_range, freq_range=freq_range, occ_data=True,
                          local_testing=local_testing)
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

    df = add_mlt_to_df(cell_corners_aacgm_lons=cell_corners_aacgm_lons,
                       cell_corners_aacgm_lats=cell_corners_aacgm_lats, df=df)

    # Get our raw x-data
    if time_units == "mlt":
        df['xdata'] = df['mlt']

    else:
        print("     Computing UTs for " + str(year) + " data...")

        ut_time = []
        for i in range(len(df)):
            ut_time_here = df['datetime'].iat[i].hour + df['datetime'].iat[i].minute / 60 + \
                           df['datetime'].iat[i].second / 3600

            if ut_time_here > 24:
                ut_time_here = ut_time_here - 24
            elif ut_time_here < 0:
                ut_time_here = ut_time_here + 24

            ut_time.append(ut_time_here)

        df['xdata'] = np.asarray(ut_time)

    # Compute a circle in axis coordinates which can be used as a boundary
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5  # 0.5 means middle of the circle
    vertices = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(vertices * radius + center)

    # Convert radar coordinates to aacgm
    radar_lat_aacgm, radar_lon_aacgm, radar_mlt = get_aacgm_coord(radar_lat, radar_lon, 0, date_time_est)
    radar_mlts = np.arange(0, 360, 1)
    radar_lats_aacgm = np.asarray([radar_lat_aacgm] * len(radar_mlts))

    if hemisphere.value == 1:
        projection = ccrs.NorthPolarStereo()
        lat_extreme = int(radar_lat_aacgm - 3)
    elif hemisphere.value == -1:
        projection = ccrs.SouthPolarStereo()
        lat_extreme = int(radar_lat_aacgm + 3)
    else:
        raise Exception("hemisphere not recognized")

    print("     Preparing the plot...")
    fig, ax = plt.subplots(figsize=[10, 6], dpi=300, nrows=1, ncols=2, subplot_kw={'projection': projection})
    # plt.subplots_adjust(left=0, right=1, bottom=0, top=0.85)
    fig.suptitle(month_string + " " + str(year) + " at " + station.upper() +
                 "; Beams " + str(beam_range[0]) + "-" + str(beam_range[1]) +
                 "; Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz", fontsize=18)

    # Apply common subplot formatting
    for i in range(ax.size):
        ax[i].set_extent([-180, 180, hemisphere.value * 90, lat_extreme], crs=ccrs.PlateCarree())
        ax[i].set_boundary(circle, transform=ax[i].transAxes)

        # Add gridlines   # Note: Labels wont draw on a circular axis
        gl = ax[i].gridlines(draw_labels=True, linestyle='--', linewidth=0.5, zorder=5)
        gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        # Print clock numbers
        text_offset_multiplier = 1.03
        ax[i].text(0, text_offset_multiplier * ax[i].get_ylim()[1], "12", ha='center', va='bottom')
        ax[i].text(0, text_offset_multiplier * ax[i].get_ylim()[0], "00", ha='center', va='top')
        ax[i].text(text_offset_multiplier * ax[i].get_xlim()[1], 0, "06", ha='left', va='center')
        ax[i].text(text_offset_multiplier * ax[i].get_xlim()[0], 0, "18", ha='right', va='center')

        # Print out the time units
        ax[i].text(ax[i].get_xlim()[1], ax[i].get_ylim()[1], time_units.upper(), ha='right', va='top')

        # Plot radar track
        ax[i].plot(radar_mlts, radar_lats_aacgm, color='k', linewidth=1, linestyle="--", transform=ccrs.Geodetic())

    # Print out echo types
    ax[0].text(ax[0].get_xlim()[0], ax[0].get_ylim()[1], "IS", ha='left', va='top')
    ax[1].text(ax[1].get_xlim()[0], ax[1].get_ylim()[1], "GS", ha='left', va='top')

    # Right now xdata is in the range 0-24, we need to put it in the range 0-360 for circular plotting
    df['xdata'] = 15 * df['xdata']

    print("     Applying correctional rotation...")
    # The 0 degree aacgm line does not line up with the 0 degree geo line, so we need to rotate everything
    zero_degree_lons_aacgm = np.asarray([0] * len(df['lat']))
    heights = np.asarray([250] * len(df['lat']))  # TODO: Figure out what to do about heights

    # returns are lat_out, lon_out, r_out
    _, zero_degree_lons_geo, _ = convert_latlon_arr(in_lat=df['lat'], in_lon=zero_degree_lons_aacgm, height=heights,
                                                    dtime=date_time_est, method_code="A2G")

    # And rotate, keeping everything in the 0-360 range
    for i in range(len(df)):
        adjusted_xdata_point = df['xdata'].iat[i] - zero_degree_lons_geo[i]

        if adjusted_xdata_point > 360:
            adjusted_xdata_point = adjusted_xdata_point - 360
        elif adjusted_xdata_point < 0:
            adjusted_xdata_point = adjusted_xdata_point + 360

        df['xdata'].iat[i] = adjusted_xdata_point

    print("     Computing binned occ rates...")
    # Compute mlt edges
    deg_mlt_per_bin = 2
    n_bins_mlt = int(360 / deg_mlt_per_bin)
    mlt_edges = np.linspace(0, 360, num=(n_bins_mlt + 1))
    delta_mlt = mlt_edges[1] - mlt_edges[0]

    # Compute latitude edges
    n_bins_lat = 90 - abs(lat_extreme)  # One bin per degree of latitude
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

        levels = 12
        levels = np.linspace(start=0, stop=1, num=(levels + 1))
        cmap = 'jet'
        # cmap = modified_jet(levels=len(levels) - 1)

        plot0 = ax[0].contourf(bin_xcenters_cyclic, bin_ycenters, contour_data_is_cyclic,
                               cmap=cmap, levels=levels, transform=ccrs.PlateCarree())
        plot1 = ax[1].contourf(bin_xcenters_cyclic, bin_ycenters, contour_data_gs_cyclic,
                               cmap=cmap, levels=levels, transform=ccrs.PlateCarree())

    elif plot_type == "pixel":

        cmap = 'jet'
        plot0 = ax[0].pcolormesh(mlt_edges, lat_edges, contour_data_is.transpose(),
                                 transform=ccrs.PlateCarree(), cmap=cmap, vmin=0, vmax=1)
        plot1 = ax[1].pcolormesh(mlt_edges, lat_edges, contour_data_gs.transpose(),
                                 transform=ccrs.PlateCarree(), cmap=cmap, vmin=0, vmax=1)

    else:
        raise Exception("plot_type not recognized")

    cbar0 = fig.colorbar(plot0, ax=ax[0], shrink=0.7, orientation="horizontal", format='%.2f')
    cbar0.ax.tick_params(labelsize=14, labelrotation=30)

    cbar1 = fig.colorbar(plot1, ax=ax[1], shrink=0.7, orientation="horizontal", format='%.2f')
    cbar1.ax.tick_params(labelsize=14, labelrotation=30)

    return df, fig


if __name__ == '__main__':
    """ Testing """

    local_testing = False

    if local_testing:
        station = "rkn"

        _, fig = occ_clock_diagram(station=station, year=2011, month_range=(11, 11), day_range=None,
                                   gate_range=(0, 74), beam_range=(6, 7), freq_range=None,
                                   plot_type='pixel', time_units='mlt',
                                   local_testing=local_testing)

        plt.show()


    else:
        station = "rkn"
        datetime_now = datetime.datetime.now()
        freq_range = (9.5, 12.5)

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        year = 2016

        for month in range(4, 7, 1):
            if year >= datetime_now.year and month > datetime_now.month:
                # No data here yet
                continue

            # Make contour plot
            _, fig = occ_clock_diagram(station=station, year=year, month_range=(month, month), day_range=None,
                                       gate_range=(0, 74), beam_range=None, freq_range=freq_range,
                                       plot_type='contour', time_units='mlt',
                                       local_testing=local_testing)

            out_fig = out_dir + "/occ_full_circle_" + station + "-" + str(year) + "-" + \
                      str(month) + "_" + str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz_contour"
            print("Saving plot as " + out_fig)
            fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)

        # # Then make a pixel plot
        # _, fig = occ_full_circle(station=station, year=year, month_range=month_range, day_range=None,
        #                          gate_range=(0, 74), beam_range=None, freq_range=freq_range,
        #                          plot_type='pixel', time_units='mlt',
        #                          local_testing=local_testing)
        #
        # out_fig = out_dir + "/occ_full_circle_" + station + "-" + str(year) + "-" + \
        #           str(month_range[0]) + "_" + str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz_pixel"
        # print("Saving plot as " + out_fig)
        # fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
