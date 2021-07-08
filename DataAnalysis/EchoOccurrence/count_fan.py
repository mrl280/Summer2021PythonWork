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

from lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.z_min_max_defaults import z_min_max_defaults
from lib.cm.modified_viridis import modified_viridis


def occ_fan(station, year_range, month_range=None, day_range=None, hour_range=None,
            gate_range=None, beam_range=None, freq_range=None,
            time_units='ut', echo_type='is', parameter=None,
            local_testing=False):
    """

    Produce a fan plot.  Can plot a simple echo count, ground scatter count, or average a fitACF parameter over the
     provided time range.

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only considers 45 km data.
            (a warning will be printed if other spatial resolution data is stripped from the dataset)
        - To check which fitACF program is being used, refer to the data readers in lib.data_getters
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
    :param freq_range: (<float>, <float>) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.
    :param time_units: str:
            The time units of the provided hour range.  Default is 'ut'.
                'ut' for universal time
                'mlt' for magnetic local time
                'lt' for local time (based on longitude)
                'lst' for local standard time (based on time zones)
                'ast' for apparent solar time (based on the apparent angular motion of the sun across the sky)
    :param echo_type: str (optional): Default is 'is' (ionospheric scatter)
            The type of echo to consider.  Either ionospheric scatter 'is' or ground scatter 'gs'.
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :param parameter: str (optional):
            Parameter to be averaged (e.g. 'v' or 'p_l')
            If omitted, then a simple echo count will be plotted.

    :return: pandas.DataFrame, matplotlib.pyplot.figure: The dataframe used, and the figure produced.
            The figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    echo_type = check_echo_type(echo_type)
    hour_range = check_hour_range(hour_range)

    print("Retrieving data...")
    df = get_data_handler(station, year_range=year_range, month_range=month_range, day_range=day_range,
                          gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                          local_testing=local_testing)

    print("Getting some hardware info...")
    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    hemisphere = this_radars_info.hemisphere
    radar_lon = this_radars_info.hardware_info.geographic.lon
    radar_lat = this_radars_info.hardware_info.geographic.lat
    radar_id = this_radars_info.hardware_info.stid

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    # Add decimal hour to df in whatever units were requested
    # Use the middle of the mid year as magnetic field estimate
    mid_year = int(year_range[0] + (year_range[1] - year_range[0]) / 2)
    date_time_est, _ = build_datetime_epoch(year=mid_year, month=6, day=15, hour=0)
    df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id, date_time_est=date_time_est)

    print("Filtering data...")
    df = df.loc[(df[time_units] >= hour_range[0]) & (df[time_units] <= hour_range[1])]
    if echo_type == 'is' and parameter is not None:
        # Then we are plotting a parameter, like v or p_l
        zmin, zmax = z_min_max_defaults(parameter)
        df = df.loc[(df[parameter] >= zmin) & (df[parameter] <= zmax)]

    df.reset_index(drop=True, inplace=True)
    df = only_keep_45km_res_data(df)

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
             label=this_radars_info.name)

    cell_corners_lats, cell_corners_lons = radar_fov(stid=radar_id, coords='geo')

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
                    scans[gate_idx, beam_idx] = statistics.median(cell_df.query('gflg == 0')[parameter])
                except statistics.StatisticsError:
                    # We can't take a median because there are no points
                    scans[gate_idx, beam_idx] = math.nan

    # Build reduced arrays containing only the cells in the specified gate/beam range
    reduced_cell_corners_lons = cell_corners_lons[gate_range[0]: gate_range[1] + 2,
                                beam_range[0]: beam_range[1] + 2]
    reduced_cell_corners_lats = cell_corners_lats[gate_range[0]: gate_range[1] + 2,
                                beam_range[0]: beam_range[1] + 2]

    print("Plotting Data...")
    if echo_type == 'gs':
        cmap = modified_viridis()
        data = ax.pcolormesh(reduced_cell_corners_lons, reduced_cell_corners_lats, grndsct_scans,
                             transform=ccrs.PlateCarree(), cmap=cmap, zorder=3)
    else:
        if parameter == 'v':
            cmap = 'seismic_r'
        else:
            cmap = modified_viridis()
        data = ax.pcolormesh(reduced_cell_corners_lons, reduced_cell_corners_lats, scans,
                             transform=ccrs.PlateCarree(), cmap=cmap, zorder=3)
    fig.colorbar(data, ax=ax)

    # plot all the beam boundary lines
    for beam_line in range(beam_range[0], beam_range[1] + 2):
        plt.plot(cell_corners_lons[gate_range[0]:gate_range[1] + 2, beam_line],
                 cell_corners_lats[gate_range[0]:gate_range[1] + 2, beam_line],
                 color='black', linewidth=0.1, transform=ccrs.Geodetic(), zorder=4)

    # plot the arcs boundary lines
    for range_ in range(gate_range[0], gate_range[1] + 2):
        plt.plot(cell_corners_lons[range_, beam_range[0]:beam_range[1] + 2],
                 cell_corners_lats[range_, beam_range[0]:beam_range[1] + 2],
                 color='black', linewidth=0.1, transform=ccrs.Geodetic(), zorder=4)

    print("Returning the figure...")
    return df, fig


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    if local_testing:
        station = "rkn"

        df, fig = occ_fan(station=station, year_range=(2011, 2011), month_range=(11, 11), day_range=None,
                          gate_range=(0, 74), beam_range=None, freq_range=None, echo_type='is', parameter='v',
                          local_testing=local_testing)
        plt.show()


    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        station = "rkn"

        df, fig = occ_fan(station=station, year_range=(2007, 2009), month_range=(2, 2), day_range=None,
                          gate_range=(0, 74), beam_range=(0, 15), freq_range=None, echo_type='is', parameter='v',
                          local_testing=local_testing)

        out_fig = out_dir + "/count_fan_" + station
        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
