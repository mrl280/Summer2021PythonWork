import statistics

import numpy as np
import pydarn

import cartopy.crs as ccrs
import cartopy.feature as cfeature


from matplotlib import pyplot as plt
from pydarn import SuperDARNRadars, radar_fov

from compute_sd_radar_overlap import compute_df_radar_overlap
from get_data_handler import get_data_handler
from only_keep_45km_res_data import only_keep_45km_res_data


def only_keep_overlap(station, df, other_station, gate_min=30):
    """

    Given a SuperDARN dataframe (df) with fit data for station, remove all data for cells that don't overlap with the
    other stations FoV.  The original dataframe is modified and also returned.

    :param station: str:
            The station whose df we are restricting, as 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param df: pandas.DataFrame
            station's dataframe.  Must contain 'gate' and 'beam' columns
    :param other_station: str
            The other station, also as a 3 character string.
    :param gate_min:
            Often the near gates will not see F region echoes.
            So, restrict both stations valid region of overlap to be only gates greater than this gate_min

    :return restricted_df:
            The provided dataframe, expect restricted to only those cells for which there is valid overlap with
             other_station
    """

    #  TODO: Figure out what to do about gate and beam ranges
    #       - right now we assume the station has 75 gates and 16 beams
    gate_range = (gate_min, 74)
    beam_range = (0, 15)

    # Quickest just to kill this stuff straight away, then we don't have to loop though it
    df = df.loc[((df['slist'] >= gate_range[0]) & (df['slist'] <= gate_range[1])) &
                ((df['bmnum'] >= beam_range[0]) & (df['bmnum'] <= beam_range[1]))]

    # Loop though all the station's cells - if there is an overlap we will keep that data - but if there is no overlap
    #  then we will remove it
    gates = np.linspace(start=gate_range[0], stop=gate_range[1], num=(gate_range[1] - gate_range[0] + 1),
                        endpoint=True, dtype=int)
    beams = np.linspace(start=beam_range[0], stop=beam_range[1], num=(beam_range[1] - beam_range[0] + 1),
                        endpoint=True, dtype=int)

    for beam in beams:
        for gate in gates:
            other_station_beam, other_station_gate = compute_df_radar_overlap(station1=station,
                                                                              station1_beam=beam, station1_gate=gate,
                                                                              station2=other_station, gate_min=gate_min)

            if other_station_beam is None or other_station_gate is None:
                # Then this cell is not overlapping, remove all data for this cell
                # print("We are killing " + station + " data at beam=" + str(beam) + " gate=" + str(gate))
                df = df.drop(df[(df['slist'] == gate) & (df['bmnum'] == beam)].index)

    df.reset_index(drop=True, inplace=True)
    return df


if __name__ == "__main__":
    """
    Testing
    
    Produces a quick fan plot to visually verify only_keep_overlap() is working as expected
    """

    # station = "dcn"
    # other_station = "mcm"

    station = "mcm"
    other_station = "dcn"

    year = 2019
    gate_range = (0, 74)
    beam_range = (0, 15)

    print("     Retrieving SuperDARN data...")
    df = get_data_handler("rkn", year_range=(year, year), month_range=None, day_range=None,
                          gate_range=gate_range, beam_range=beam_range, freq_range=None,
                          local_testing=True, even_odd_days=None)
    df = only_keep_45km_res_data(df)
    print(df.keys())
    print("The dataframe is of length " + str(len(df)))

    # print("Here is a list of remaining gates:")
    # print(df['slist'].unique())
    #
    # print("Here is a list of remaining beams:")
    # print(df['bmnum'].unique())

    # Produce a quick fan plot to verify the results
    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    hemisphere = this_radars_info.hemisphere
    radar_lon = this_radars_info.hardware_info.geographic.lon
    radar_lat = this_radars_info.hardware_info.geographic.lat
    radar_id = this_radars_info.hardware_info.stid

    print("Preparing the plot...")
    parameter = 'v'  # Parameter to plot on the fan

    # Prepare the figure
    fig = plt.figure(figsize=(10, 5), constrained_layout=True, dpi=300)
    gs = fig.add_gridspec(ncols=2, nrows=1)

    if hemisphere.value == 1:
        proj = ccrs.NorthPolarStereo()  # Northern hemisphere

        min_lat = 34  # deg
        ax1 = fig.add_subplot(gs[0, 0], projection=proj)
        ax2 = fig.add_subplot(gs[0, 1], projection=proj)

        ax1.set_extent([-180, 180, 90, min_lat], crs=ccrs.PlateCarree())
        ax2.set_extent([-180, 180, 90, min_lat], crs=ccrs.PlateCarree())

    elif hemisphere.value == -1:
        proj = ccrs.SouthPolarStereo()  # Southern hemisphere
        max_lat = -34  # deg

        ax1 = fig.add_subplot(gs[0, 0], projection=proj)
        ax2 = fig.add_subplot(gs[0, 1], projection=proj)

        ax1.set_extent([-180, 180, -90, max_lat], crs=ccrs.PlateCarree())
        ax2.set_extent([-180, 180, -90, max_lat], crs=ccrs.PlateCarree())


    else:
        raise Exception("Error: hemisphere not recognized")

    ax1.set_title("Before only_keep_overlap()")
    ax2.set_title("After only_keep_overlap()")

    for ax in [ax1, ax2]:
        ax.gridlines()
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.LAND)


    """ On the first plot, plot the fan before restriction """
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

           # Average the provided parameter
            try:
                scans[gate_idx, beam_idx] = statistics.median(cell_df.query('gflg == 0')[parameter])
            except statistics.StatisticsError:
                # We can't take a median because there are no points
                scans[gate_idx, beam_idx] = np.nan

    # Build reduced arrays containing only the cells in the specified gate/beam range
    reduced_cell_corners_lons = cell_corners_lons[gate_range[0]: gate_range[1] + 2,
                                beam_range[0]: beam_range[1] + 2]
    reduced_cell_corners_lats = cell_corners_lats[gate_range[0]: gate_range[1] + 2,
                                beam_range[0]: beam_range[1] + 2]

    cmap = 'seismic_r'
    data = ax1.pcolormesh(reduced_cell_corners_lons, reduced_cell_corners_lats, scans,
                          transform=ccrs.PlateCarree(), cmap=cmap, zorder=3)
    # fig.colorbar(data, ax=ax1)


    print("Running only_keep_overlap().. ")
    df = only_keep_overlap(station=station, df=df, other_station=other_station, gate_min=20)
    print("The dataframe is now of length " + str(len(df)))


    """ On the second plot, plot the fan after restriction """
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

            # Average the provided parameter
            try:
                scans[gate_idx, beam_idx] = statistics.median(cell_df.query('gflg == 0')[parameter])
            except statistics.StatisticsError:
                # We can't take a median because there are no points
                scans[gate_idx, beam_idx] = np.nan

    # Build reduced arrays containing only the cells in the specified gate/beam range
    reduced_cell_corners_lons = cell_corners_lons[gate_range[0]: gate_range[1] + 2,
                                beam_range[0]: beam_range[1] + 2]
    reduced_cell_corners_lats = cell_corners_lats[gate_range[0]: gate_range[1] + 2,
                                beam_range[0]: beam_range[1] + 2]

    cmap = 'seismic_r'
    data = ax2.pcolormesh(reduced_cell_corners_lons, reduced_cell_corners_lats, scans,
                          transform=ccrs.PlateCarree(), cmap=cmap, zorder=3)
    # fig.colorbar(data, ax=ax1)


    """ Plot the radars fans """
    colours = {station: 'red', other_station: 'blue'}
    for ax in [ax1, ax2]:
        for station_here in [station, other_station]:
            station_info = all_radars_info.radars[pydarn.read_hdw_file(station_here).stid]
            radar_id = station_info.hardware_info.stid
            radar_lon = station_info.hardware_info.geographic.lon
            radar_lat = station_info.hardware_info.geographic.lat
            print("Plotting the fan for " + station_here)
            print("     Radar lat: " + str(radar_lat))
            print("     Radar lon: " + str(radar_lon))

            # Plot the radar as a dot
            ax.plot([radar_lon, radar_lon], [radar_lat, radar_lat], marker="o", markersize=3, transform=ccrs.Geodetic(),
                     label=station_info.name, color=colours[station_here])

            cell_corners_lats, cell_corners_lons = radar_fov(stid=radar_id, coords='geo')

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

    plt.show()
