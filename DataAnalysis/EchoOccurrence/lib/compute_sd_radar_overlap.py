import math
import pathlib
import pydarn

import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import numpy as np
import cartopy.feature as cfeature
import datetime as datetime

from matplotlib import pyplot as plt
from pydarn import radar_fov, SuperDARNRadars

try:
    # Assume we are importing from one level up
    from lib.add_mlt_to_df import centroid
except ImportError:
    # Needed for testing
    from add_mlt_to_df import centroid


def compute_df_radar_overlap(station1, station1_beam, station1_gate, station2, gate_range=(20, 74), beam_range=(0, 15)):
    """

    Given a SuperDARN radar cell - compute the overlapping cell from a second radar.

    For example, if you wanted to know which MCM radar cell corresponded to the cell defined by beam 10 gate 40 at DCE:
        mcm_beam, mcm_gate = compute_df_radar_overlap(station1=dce, station1_gate=40, station1_beam=10, station2=mcm)

    Results are only valid when both radars are operating in the standard 45 km mode.

    :param station1: str:
            The station whose gate and beam we know, as 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param station1_beam: int:
            The beam of the first station's cell
    :param station1_gate: int:
            The gate of the first station's cell
    :param station2:
            The second station we want to match up with
    :param gate_range: (int, int):
            The range of allowable gates.  Gate matches from outside this range (at either radar) will not be considered
            Often the near gates will not see F region echoes.
            So, you probably want to start around gate 20.
    :param beam_range: (int, int):
            The range of allowable beams.  Beam matches from outside this range (at either radar) will not be considered

    :return station2_beam, station2_gate: int, int:
            The beam and gate of station2 that overlap with the provided station1_beam and station1_gate
            If no valid overlap is found - then None, None is returned
    """

    if station1_gate < gate_range[0] or station1_gate > gate_range[1] or \
            station1_beam < beam_range[0] or station1_beam > beam_range[1]:
        # We are outside of the allowable gate/beam range
        return None, None

    all_radars_info = SuperDARNRadars()
    station1_stid = pydarn.read_hdw_file(station1).stid
    station2_stid = pydarn.read_hdw_file(station2).stid
    station1_info = all_radars_info.radars[station1_stid]
    station2_info = all_radars_info.radars[station2_stid]

    reference_hemisphere = station1_info.hemisphere
    if station2_info.hemisphere != reference_hemisphere:
        # The stations are in different hemispheres - they have no overlap
        return None, None

    # Compute lon and lat for the centroid of the given cell
    station1_cell_corners_lats, station1_cell_corners_lons = radar_fov(stid=station1_stid, coords='geo')
    station1_cell_centroid_lon, station1_cell_centroid_lat = centroid(
        [(station1_cell_corners_lons[station1_gate, station1_beam],
          station1_cell_corners_lats[station1_gate, station1_beam]),
         (station1_cell_corners_lons[station1_gate + 1, station1_beam],
          station1_cell_corners_lats[station1_gate + 1, station1_beam]),
         (station1_cell_corners_lons[station1_gate, station1_beam + 1],
          station1_cell_corners_lats[station1_gate, station1_beam + 1]),
         (station1_cell_corners_lons[station1_gate + 1, station1_beam + 1],
          station1_cell_corners_lats[station1_gate + 1, station1_beam + 1])])

    # We need to compute all of station2's centroids, and then compute the distance from each of those centroids to
    #  station1's centroid
    station2_cell_corners_lats, station2_cell_corners_lons = radar_fov(stid=station2_stid, coords='geo')
    fan_shape = station2_cell_corners_lons.shape

    # Save the distances so we can find the smallest one
    station2_cell_centroid_distances = np.full(shape=(fan_shape[0] - 1, fan_shape[1] - 1), fill_value=np.inf)

    for gate_corner in range(fan_shape[0] - 1):
        for beam_corner in range(fan_shape[1] - 1):
            cent_lon, cent_lat = centroid([(station2_cell_corners_lons[gate_corner, beam_corner],
                                            station2_cell_corners_lats[gate_corner, beam_corner]),
                                           (station2_cell_corners_lons[gate_corner + 1, beam_corner],
                                            station2_cell_corners_lats[gate_corner + 1, beam_corner]),
                                           (station2_cell_corners_lons[gate_corner, beam_corner + 1],
                                            station2_cell_corners_lats[gate_corner, beam_corner + 1]),
                                           (station2_cell_corners_lons[gate_corner + 1, beam_corner + 1],
                                            station2_cell_corners_lats[gate_corner + 1, beam_corner + 1])])

            # Find the Euclidean distance between this station2 centroid and the provided station1's centroid
            distance = math.sqrt((cent_lon - station1_cell_centroid_lon) ** 2
                                 + (cent_lat - station1_cell_centroid_lat) ** 2)

            station2_cell_centroid_distances[gate_corner, beam_corner] = distance

    # Look through the array of distances - we need to find the smallest one
    station2_gate, station2_beam = np.unravel_index(station2_cell_centroid_distances.argmin(),
                                                    station2_cell_centroid_distances.shape)

    if station2_gate < gate_range[0] or station2_gate > gate_range[1] or \
            station2_beam < beam_range[0] or station2_beam > beam_range[1]:
        # We are outside of the allowable gate/beam range
        return None, None

    if station2_cell_centroid_distances[station2_gate, station2_beam] > 1.5:
        # TODO: Optimize the acceptable difference - using degrees for difference might not be acceptable because 1 deg
        #  cells are smaller near the poles
        # We found something, but what we found is unacceptably far away
        return None, None

    else:
        return station2_beam, station2_gate


if __name__ == '__main__':
    """ 
    Testing - Produces a fan map to visually confirm the validity of the results
    """

    SAVE_PLOT = False

    # station1 = "mcm"
    # station1_gate = 60
    # station1_beam = 9
    # station2 = "dcn"
    # gate_min = 20

    # station1 = "sas"
    # station1_gate = 68
    # station1_beam = 0
    # station2 = "cly"
    # gate_min = 20

    station1 = "dcn"
    station1_gate = 40
    station1_beam = 9
    station2 = "mcm"
    gate_min = 20

    station2_beam, station2_gate = compute_df_radar_overlap(station1=station1, station1_gate=station1_gate,
                                                            station1_beam=station1_beam, station2=station2,
                                                            gate_min=gate_min)

    print("Found station2_beam=" + str(station2_beam))
    print("Found station2_gate=" + str(station2_gate))

    all_radars_info = SuperDARNRadars()
    station1_info = all_radars_info.radars[pydarn.read_hdw_file(station1).stid]
    reference_hemisphere = station1_info.hemisphere

    # Produce a fan map to visually confirm the validity of the results
    fig = plt.figure(figsize=(5, 5), dpi=300, constrained_layout=True)
    lat_extreme = 40 * reference_hemisphere.value  # deg
    if reference_hemisphere.value == 1:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
    elif reference_hemisphere.value == -1:
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())
    else:
        raise Exception("hemisphere not recognized")

    ax.set_extent([-180, 180, reference_hemisphere.value * 90, lat_extreme], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)

    # Add gridlines
    gl = ax.gridlines(draw_labels=True, linestyle='--', color='white')
    gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135])

    # Move the longitude labels off the 135 deg line and onto the 45 degree line
    plt.draw()
    for ea in gl.label_artists:
        if ea[1] == False:
            tx = ea[2]
            xy = tx.get_position()
            # print(xy)
            if xy[0] == 135:
                tx.set_position([45, xy[1]])

    gate_range = (gate_min, 74)
    beam_range = (0, 15)

    # Plot the fan for each station
    colours = {station1: 'red', station2: 'blue'}
    for station in [station1, station2]:
        station_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]
        radar_id = station_info.hardware_info.stid
        radar_lon = station_info.hardware_info.geographic.lon
        radar_lat = station_info.hardware_info.geographic.lat
        print("Looking at " + station)
        print("     Radar lat: " + str(radar_lat))
        print("     Radar lon: " + str(radar_lon))

        # Plot the radar as a dot
        plt.plot([radar_lon, radar_lon], [radar_lat, radar_lat], marker="o", markersize=3, transform=ccrs.Geodetic(),
                 label=station_info.name, color=colours[station])

        cell_corners_lats, cell_corners_lons = radar_fov(stid=radar_id, coords='geo', date=datetime.datetime.now())

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

    # Plot the centroid of the first radars cell
    cell_corners_lats, cell_corners_lons = radar_fov(stid=pydarn.read_hdw_file(station1).stid, coords='geo')
    station1_cell_lon, station1_cell_lat = centroid([(cell_corners_lons[station1_gate, station1_beam],
                                                      cell_corners_lats[station1_gate, station1_beam]),
                                                     (cell_corners_lons[station1_gate + 1, station1_beam],
                                                      cell_corners_lats[station1_gate + 1, station1_beam]),
                                                     (cell_corners_lons[station1_gate, station1_beam + 1],
                                                      cell_corners_lats[station1_gate, station1_beam + 1]),
                                                     (cell_corners_lons[station1_gate + 1, station1_beam + 1],
                                                      cell_corners_lats[station1_gate + 1, station1_beam + 1])])

    plt.plot([station1_cell_lon, station1_cell_lon], [station1_cell_lat, station1_cell_lat], marker=".", markersize=1,
             transform=ccrs.Geodetic(), color=colours[station1])

    try:
        # Plot the centroid of the returned matching cell
        cell_corners_lats, cell_corners_lons = radar_fov(stid=pydarn.read_hdw_file(station2).stid, coords='geo')
        station2_cell_lon, station2_cell_lat = centroid([(cell_corners_lons[station2_gate, station2_beam],
                                                          cell_corners_lats[station2_gate, station2_beam]),
                                                         (cell_corners_lons[station2_gate + 1, station2_beam],
                                                          cell_corners_lats[station2_gate + 1, station2_beam]),
                                                         (cell_corners_lons[station2_gate, station2_beam + 1],
                                                          cell_corners_lats[station2_gate, station2_beam + 1]),
                                                         (cell_corners_lons[station2_gate + 1, station2_beam + 1],
                                                          cell_corners_lats[station2_gate + 1, station2_beam + 1])])

        plt.plot([station2_cell_lon, station2_cell_lon], [station2_cell_lat, station2_cell_lat], marker=".",
                 markersize=1, transform=ccrs.Geodetic(), color=colours[station2])
    except BaseException as e:
        print("Unable to plot the found matching cell")
        print(e)
        pass

    plt.legend(loc='best')

    loc_root = str((pathlib.Path().parent.absolute().parent.absolute()))
    out_dir = loc_root + "/out"
    out_file = out_dir + "/testing_SuperDARN_overlap"

    if SAVE_PLOT:
        print("Saving figure as " + out_file)
        fig.savefig(out_file + ".jpg", format='jpg', dpi=300)

    plt.show()
    plt.close(fig)
