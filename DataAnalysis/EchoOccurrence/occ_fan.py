import math
import pathlib
import statistics
import pydarn

import matplotlib.path as mpath
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib import pyplot as plt
from pydarn import SuperDARNRadars, radar_fov

from lib.z_min_max_defaults import z_min_max_defaults
from lib.cm.modified_viridis import modified_viridis
from lib.data_getters.get_data import get_data
from lib.data_getters.get_local_dummy_data import get_local_dummy_data
from lib.range_checkers import *


def occ_fan(station, year_range, month_range=None, day_range=None, hour_range=None, gate_range=None, beam_range=None,
            local_testing=False, parameter=None, plot_ground_scat=False):
    """

    Produce a fan plot.  Can plot a simple echo count, ground scatter count, or average a fitACF parameter over the
     provided time range.

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Does not distinguish frequency
        - Only considers 45 km data.
            (a warning will be printed if other spatial resolution data is stripped from the dataset)
        - This program uses fitACF 3.0 data.  To change this, modify the source code.
        - All times and dates are assumed UT

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year_range: (<int>, <int>):
            Inclusive. The year range to consider.
    :param month_range: (<int>, <int>) (optional):
            Inclusive. The months of the year to consider.  If omitted (or None), then all days will be considered.
    :param day_range: (<int>, <int>) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param hour_range: (<int>, <int>) (optional):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param gate_range: (<int>, <int>) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :param parameter: str (optional):
            Parameter to be averaged (e.g. 'v' or 'p_l')
            If omitted, then a simple echo count will be plotted.
    :param plot_ground_scat: bool (optional)
            Set this to true if you would like to plot ground scatter counts.  Default is False
            If plot_ground_scat is set to True, then :param parameter is ignored.
    :return: matplotlib.pyplot.figure, 2d np.array: The figure and the scan data plotted.
            The figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    if parameter is not None:
        zmin, zmax = z_min_max_defaults(parameter)

    year_range = check_year_range(year_range)
    month_range = check_month_range(month_range)
    day_range = check_day_range(day_range)
    hour_range = check_hour_range(hour_range)

    if isinstance(station, str):
        hdw_info = pydarn.read_hdw_file(station)  # Get the hardware file, there is lots of good stuff in there
    else:
        raise Exception("Error: Please enter the station as a 3 character string.  e.g. 'rkn'")

    beam_range = check_beam_range(beam_range, hdw_info)
    gate_range = check_gate_range(gate_range, hdw_info)

    radar_lon = hdw_info.geographic.lon
    radar_lat = hdw_info.geographic.lat
    radar_id = hdw_info.stid

    all_radars_info = SuperDARNRadars()
    additional_radar_info = all_radars_info.radars[radar_id]  # Grab some additional radar info
    hemisphere = additional_radar_info.hemisphere

    print("Retrieving data...")
    if local_testing:
        # Just read in some test data
        warnings.warn("Running in local testing mode, just going to use local dummy data", category=Warning)
        # df = get_local_dummy_data(station=station, year=2011, month=9, day=29, start_hour_UT=0, end_hour_UT=23)
        df = get_local_dummy_data(station=station, year=2014, month=3, day=3, start_hour_UT=0, end_hour_UT=12)
        # print(df.keys())
    else:
        df = get_data(station=station, year_range=year_range, month_range=month_range, day_range=day_range,
                      hour_range=hour_range, gate_range=gate_range, beam_range=beam_range)

    df = df.loc[(df['p_l'] >= 3)]  # Restrict to points with at least 3 dB

    if not plot_ground_scat and parameter is not None:
        df = df.loc[(df[parameter] >= zmin) & (df[parameter] <= zmax)]

    df.reset_index(drop=True, inplace=True)

    # Restrict to 45 km mode data, you need to modify the fan if you want to use data of a different spatial resolution
    number_of_records_before_spatial_resolution_check = df.shape[0]
    df = df.loc[(df['frang'] == 180) & (df['rsep'] == 45)]
    df.reset_index(drop=True, inplace=True)
    number_of_records_after_spatial_resolution_check = df.shape[0]
    if number_of_records_before_spatial_resolution_check > number_of_records_after_spatial_resolution_check:
        warnings.warn("Not all data within the specified range is 45 km resolution data.  "
                      "All non-45 km data was discarded.", category=Warning)

    print("Preparing the plot...")
    # Prepare the figure
    fig = plt.figure(figsize=(5, 5), dpi=300)
    if hemisphere.value == 1:
        # Northern hemisphere
        min_lat = 37  # deg
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, 90, min_lat], crs=ccrs.PlateCarree())
    elif hemisphere.value == -1:
        # Southern hemisphere
        max_lat = -34  # deg
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90, max_lat], crs=ccrs.PlateCarree())
    else:
        raise Exception("Error: hemisphere not recognized")

    ax.gridlines()
    # ax.stock_img()
    # ax.add_feature(cfeature.COASTLINE)
    # ax.add_feature(cfeature.BORDERS, linestyle="-", linewidth=0.5)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)

    # Compute a circle in axis coordinates which can be used as a boundary
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    vertices = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(vertices * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    # Plot the radar as a red dot
    plt.plot([radar_lon, radar_lon], [radar_lat, radar_lat], 'ro', markersize=1, transform=ccrs.Geodetic(),
             label=additional_radar_info.name)

    beam_corners_lats, beam_corners_lons = radar_fov(radar_id, coords='geo')  # Get the radar field of view

    print("Computing scan...")

    num_gates = (gate_range[1] + 1) - gate_range[0]
    num_beams = (beam_range[1] + 1) - beam_range[0]

    # Loop through all the gate/beam cells and build the scans
    # First index will be gates, the second will be beams
    scans = np.zeros((num_gates, num_beams))
    grndsct_scans = np.zeros((num_gates, num_beams))
    for gate_idx in range(num_gates):
        for beam_idx in range(num_beams):
            gate = gate_range[0] + gate_idx
            beam = beam_range[0] + beam_idx
            # print("Gate: " + str(gate) + ", beam : " + str(beam))
            cell_df = df[(df['slist'] == gate) & (df['bmnum'] == beam)]
            grndsct_scans[gate_idx, beam_idx] = cell_df[(cell_df['gflg'] == 1)].shape[0]

            if parameter is None:
                # We want a simple echo count
                scans[gate_idx, beam_idx] = cell_df[(cell_df['gflg'] == 0)].shape[0]
            else:
                # Otherwise average the provided parameter
                try:
                    scans[gate_idx, beam_idx] = statistics.mean(cell_df.query('gflg == 0')[parameter])
                except statistics.StatisticsError:
                    # We can't take a median because there are no points
                    scans[gate_idx, beam_idx] = math.nan

    # Build reduced arrays containing only the cells in the specified gate/beam range
    reduced_beam_corners_lons = beam_corners_lons[gate_range[0]: gate_range[1] + 2,
                                                  beam_range[0]: beam_range[1] + 2]
    reduced_beam_corners_lats = beam_corners_lats[gate_range[0]: gate_range[1] + 2,
                                                  beam_range[0]: beam_range[1] + 2]

    print("Plotting Data...")
    if plot_ground_scat:
        cmap = modified_viridis()
        data = ax.pcolormesh(reduced_beam_corners_lons, reduced_beam_corners_lats, grndsct_scans,
                             transform=ccrs.PlateCarree(), cmap=cmap, zorder=3)
    else:
        if parameter == 'v':
            cmap = 'seismic_r'
        else:
            cmap = modified_viridis()
        data = ax.pcolormesh(reduced_beam_corners_lons, reduced_beam_corners_lats, scans,
                             transform=ccrs.PlateCarree(), cmap=cmap, zorder=3)
    fig.colorbar(data, ax=ax)

    # plot all the beam boundary lines
    for beam_line in range(beam_range[0], beam_range[1] + 2):
        plt.plot(beam_corners_lons[gate_range[0]:gate_range[1] + 2, beam_line],
                 beam_corners_lats[gate_range[0]:gate_range[1] + 2, beam_line],
                 color='black', linewidth=0.1, transform=ccrs.Geodetic(), zorder=4)

    # plot the arcs boundary lines
    for range_ in range(gate_range[0], gate_range[1] + 2):
        plt.plot(beam_corners_lons[range_, beam_range[0]:beam_range[1] + 2],
                 beam_corners_lats[range_, beam_range[0]:beam_range[1] + 2],
                 color='black', linewidth=0.1, transform=ccrs.Geodetic(), zorder=4)

    print("Returning the figure and scan...")
    if plot_ground_scat:
        return fig, grndsct_scans
    else:
        return fig, scans


if __name__ == '__main__':
    """ Testing """

    station = "rkn"
    fig, scans = occ_fan(station=station, year_range=(2007, 2009), month_range=(2, 2), day_range=None,
                         gate_range=(0, 74), beam_range=(0, 15), plot_ground_scat=False, parameter='v',
                         local_testing=True)

    loc_root = str((pathlib.Path().parent.absolute()))
    out_dir = loc_root + "/out"
    out_file = out_dir + "/occ_fan_" + station
    print("Saving plot as " + out_file)

    # fig.savefig(out_file + ".jpg", format='jpg', dpi=300)
    plt.show()
