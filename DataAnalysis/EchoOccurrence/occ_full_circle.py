import math
import pathlib
import pydarn

from aacgmv2 import convert_latlon_arr, convert_latlon, get_aacgm_coord
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import radar_fov, SuperDARNRadars
from scipy import stats

import cartopy.crs as ccrs
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import numpy as np

from DataAnalysis.EchoOccurrence.lib.cm.modified_jet import modified_jet
from lib.add_mlt_to_df import add_mlt_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def occ_full_circle(station, year, month_range=None, day_range=None, gate_range=None, beam_range=None,
                    time_units='mlt', plot_type='contour', local_testing=False):
    """

    Produce a full circle stereographic plot in either ut or mlt (12 at the top).
    Can plot a simple echo count, ground scatter count, or average a fitACF parameter over the provided time range.

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Does not distinguish frequency
        - Only considers 45 km data.
            (a warning will be printed if other spatial resolution data is stripped from the dataset)
        - This program uses fitACF 3.0 data.  To change this, modify the source code.

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

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid
    hemisphere = this_radars_info.hemisphere
    radar_lon = this_radars_info.hardware_info.geographic.lon
    radar_lat = this_radars_info.hardware_info.geographic.lat

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    print("Retrieving data...")
    df = get_data_handler(station, year_range=(year, year), month_range=month_range, day_range=day_range,
                          gate_range=gate_range, beam_range=beam_range, occ_data=True,
                          local_testing=local_testing)
    df = only_keep_45km_res_data(df)

    # To compute mlt we need longitudes.. use the middle of the month as magnetic field estimate
    date_time_est, _ = build_datetime_epoch(year, 6, 15, 0)

    print("Computing MLTs for " + str(year) + " data...")
    cell_corners_aacgm_lats, cell_corners_aacgm_lons = radar_fov(stid=radar_id, coords='aacgm', date=date_time_est)

    df = add_mlt_to_df(cell_corners_aacgm_lons=cell_corners_aacgm_lons,
                       cell_corners_aacgm_lats=cell_corners_aacgm_lats, df=df)

    # Get our raw x-data
    if time_units == "mlt":
        df['xdata'] = df['mlt']

    else:
        print("Computing UTs for " + str(year) + " data...")

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

    print("Preparing the plot...")
    fig = plt.figure(figsize=(5, 5), dpi=300)
    lat_extreme = 40 * hemisphere.value  # deg
    if hemisphere.value == 1:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
    elif hemisphere.value == -1:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())
    else:
        raise Exception("hemisphere not recognized")

    ax.set_extent([-180, 180, hemisphere.value * 90, lat_extreme], crs=ccrs.PlateCarree())

    # Compute a circle in axis coordinates which can be used as a boundary
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    vertices = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(vertices * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    # Add gridlines and mlt labels
    text_offset_multiplier = 1.03
    gl = ax.gridlines(draw_labels=True, linestyle='--', zorder=5)
    gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Print clock numbers
    ax.text(0, text_offset_multiplier * ax.get_ylim()[1], "12", ha='center', va='bottom')
    ax.text(0, text_offset_multiplier * ax.get_ylim()[0], "00", ha='center', va='top')
    ax.text(text_offset_multiplier * ax.get_xlim()[1], 0, "06", ha='left', va='center')
    ax.text(text_offset_multiplier * ax.get_xlim()[0], 0, "18", ha='right', va='center')

    # Print time units
    ax.text(text_offset_multiplier * ax.get_xlim()[0],
            text_offset_multiplier * ax.get_ylim()[1], time_units.upper(), ha='left', va='bottom')

    # Convert radar coordinates to aacgm and plot radar track
    radar_lat_aacgm, radar_lon_aacgm, radar_mlt = get_aacgm_coord(radar_lat, radar_lon, 0, date_time_est)
    radar_mlts = np.arange(0, 360, 1)
    radar_lats_aacgm = np.asarray([radar_lat_aacgm] * len(radar_mlts))
    ax.plot(radar_mlts, radar_lats_aacgm, color='k', linewidth=0.5, linestyle="--", transform=ccrs.Geodetic(),
            label="Radar Path")

    # Right now xdata is in the range 0-24, we need to put it in the range 0-360 for circular plotting
    df['xdata'] = 15 * df['xdata']

    print("Computing binned occ rates...")
    # Compute mlt edges
    deg_mlt_per_bin = 2
    n_bins_mlt = int(360 / deg_mlt_per_bin)
    mlt_edges = np.linspace(0, 360, num=(n_bins_mlt + 1))
    delta_mlt = mlt_edges[1] - mlt_edges[0]

    # Compute latitude edges
    n_bins_lat = 90 - abs(lat_extreme)  # One bin per degree of latitude
    lat_edges = np.linspace(lat_extreme, 90, num=(n_bins_lat + 1))
    delta_lat = lat_edges[1] - lat_edges[0]

    contour_data = np.empty(shape=(n_bins_mlt, n_bins_lat))
    contour_data[:] = math.nan

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
                contour_data[mlt_idx][lat_idx] = sum(df_mlt_lat['good_echo']) / len(df_mlt_lat)
            except ZeroDivisionError:
                # There are no point in this interval
                contour_data[mlt_idx, lat_idx] = math.nan
            except BaseException as e:
                print("MLT index: " + str(mlt_idx))
                print("LAT index: " + str(lat_idx))
                raise e

    contour_range = [[0, 360], [hemisphere.value * 37, hemisphere.value * 90]]

    # Compute bin centers
    bin_xwidth = (mlt_edges[1] - mlt_edges[0])
    bin_ywidth = (lat_edges[1] - lat_edges[0])
    bin_xcenters = mlt_edges[1:] - bin_xwidth / 2
    bin_ycenters = lat_edges[1:] - bin_ywidth / 2

    levels = 12
    levels = np.linspace(start=0, stop=1, num=(levels + 1))
    if plot_type == "contour":
        plot = ax.contourf(bin_xcenters, bin_ycenters, contour_data.transpose(), cmap='jet', levels=levels, transform=ccrs.PlateCarree())

    elif plot_type == "pixel":
        plot = ax.imshow(np.flip(contour_data.transpose(), axis=0), aspect='auto', cmap="jet", transform=ccrs.PlateCarree())

    else:
        raise Exception("plot_type not recognized")

    cbar = fig.colorbar(plot, ax=ax, shrink=0.75, orientation="horizontal", format='%.2f')
    cbar.ax.tick_params(labelsize=14, labelrotation=30)

    return df, fig


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    if local_testing:
        station = "rkn"

        df, fig = occ_full_circle(station=station, year=2011, month_range=None, day_range=None, time_units='mlt',
                                  gate_range=(0, 74), beam_range=None, plot_type='contour',
                                  local_testing=local_testing)

        plt.show()


    else:
        station = "rkn"
        year = 2011
        month = 2

        df, fig = occ_full_circle(station=station, year=year, month_range=None, day_range=None, time_units='mlt',
                                  gate_range=(0, 74), beam_range=(13, 15), plot_type='contour',
                                  local_testing=local_testing)

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        out_file = out_dir + "/occ_full_circle" + station + str(year) + str(month)
        print("Saving plot as " + out_file)
        fig.savefig(out_file + ".jpg", format='jpg', dpi=300)

        out_dir = loc_root + "/data"
        out_file = out_dir + "/occ_full_circle_df_" + station + str(year) + str(month) + ".pkl"
        print("Pickling df as " + out_file)
        df.to_pickle(out_file)
