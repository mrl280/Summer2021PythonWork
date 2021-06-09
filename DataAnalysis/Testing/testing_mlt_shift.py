import aacgmv2
import pydarn
import datetime

import numpy as np
import cartopy.crs as ccrs
import matplotlib.path as mpath
import matplotlib.ticker as mticker

from matplotlib import pyplot as plt
from pydarn import radar_fov

if __name__ == "__main__":
    """
    Look at difference between shifting the whole fan and computing mlt individually.
    
    Looks like difference is negligible
    """

    # Look at RKN to see if the mlt shift produces the same as computing mlt for each
    station = "han"
    beam_range = (0, 15)
    gate_range = (0, 74)
    date = datetime.datetime.now()

    hdw_info = pydarn.read_hdw_file(station)
    radar_lon = hdw_info.geographic.lon
    radar_lat = hdw_info.geographic.lat
    radar_id = hdw_info.stid

    # Prepare the figure
    min_lat = 37  # deg
    fig = plt.figure(figsize=(5, 5), dpi=600)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180, 180, 90, min_lat], crs=ccrs.PlateCarree())

    # ax.gridlines()
    # ax.stock_img()
    # ax.add_feature(cfeature.COASTLINE)
    # ax.add_feature(cfeature.BORDERS, linestyle="-", linewidth=0.5)
    # ax.add_feature(cfeature.OCEAN)
    # ax.add_feature(cfeature.LAND)

    beam_corners_aacgm_lats, beam_corners_aacgm_lons = \
        radar_fov(stid=hdw_info.stid, coords='aacgm', date=date)
    fan_shape = beam_corners_aacgm_lons.shape

    """ Shift the whole fan - like in pydarn """
    # beam_corners_mlts = np.zeros((fan_shape[0], fan_shape[1]))
    mltshift = beam_corners_aacgm_lons[0, 0] - (aacgmv2.convert_mlt(beam_corners_aacgm_lons[0, 0], date) * 15)
    beam_corners_mlts_1 = beam_corners_aacgm_lons - mltshift

    # plot all the beam boundary lines
    for beam_line in range(beam_range[0], beam_range[1] + 2):
        plt.plot(beam_corners_mlts_1[gate_range[0]:gate_range[1] + 2, beam_line],
                 beam_corners_aacgm_lats[gate_range[0]:gate_range[1] + 2, beam_line],
                 color='black', linewidth=0.5, transform=ccrs.Geodetic(), zorder=4)

    # plot the arcs boundary lines
    for range_ in range(gate_range[0], gate_range[1] + 2):
        plt.plot(beam_corners_mlts_1[range_, beam_range[0]:beam_range[1] + 2],
                 beam_corners_aacgm_lats[range_, beam_range[0]:beam_range[1] + 2],
                 color='black', linewidth=0.5, transform=ccrs.Geodetic(), zorder=4)

    """ Try to compute each point individually """
    beam_corners_mlts_2 = np.zeros((fan_shape[0], fan_shape[1]))
    for i in range(fan_shape[0]):
        for j in range(fan_shape[1]):
            beam_corners_mlts_2[i, j] = aacgmv2.convert_mlt(beam_corners_aacgm_lons[i, j], date) * 15

    # plot all the beam boundary lines
    for beam_line in range(beam_range[0], beam_range[1] + 2):
        plt.plot(beam_corners_mlts_2[gate_range[0]:gate_range[1] + 2, beam_line],
                 beam_corners_aacgm_lats[gate_range[0]:gate_range[1] + 2, beam_line],
                 color='red', linewidth=0.1, transform=ccrs.Geodetic(), zorder=4)

    # plot the arcs boundary lines
    for range_ in range(gate_range[0], gate_range[1] + 2):
        plt.plot(beam_corners_mlts_2[range_, beam_range[0]:beam_range[1] + 2],
                 beam_corners_aacgm_lats[range_, beam_range[0]:beam_range[1] + 2],
                 color='red', linewidth=0.1, transform=ccrs.Geodetic(), zorder=4)

    """ Format Plot """
    # Compute a circle in axis coordinates which can be used as a boundary
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    vertices = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(vertices * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    # Add gridlines and mlt labels
    text_offset_multiplier = 1.03
    gl = ax.gridlines(draw_labels=False)
    gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135])
    ax.text(0, text_offset_multiplier * ax.get_ylim()[1], "12", ha='center', va='bottom')
    ax.text(0, text_offset_multiplier * ax.get_ylim()[0], "00", ha='center', va='top')
    ax.text(text_offset_multiplier * ax.get_xlim()[1], 0, "06", ha='left', va='center')
    ax.text(text_offset_multiplier * ax.get_xlim()[0], 0, "18", ha='right', va='center')

    """ Look at difference """
    difference = beam_corners_mlts_2 - beam_corners_mlts_1
    difference = np.ndarray.flatten(difference)
    # print(np.sort(difference))

    plt.show()
