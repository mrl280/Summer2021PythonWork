import warnings
import calendar

import pydarn
import time

import matplotlib.path as mpath
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import pyplot as plt
from pydarn import SuperDARNRadars, radar_fov


def check_month_range(month_range):
    """
    :param month_range: (<int>, <int>) or None: The month range to check
    :return: (<int>, <int>): A suitable month range
    """
    if month_range is None:
        month_range = (1, 12)  # Assume we want all months
    if month_range[0] < 1:
        month_range = (1, month_range[1])
    if month_range[1] > 12:
        month_range = (month_range[0], 12)
    return month_range


def check_day_range(day_range):
    """
    :param day_range: (<int>, <int>) or None: The day range to check
    :return: (<int>, <int>): A suitable day range
    """
    if day_range is None:
        day_range = (1, 31)  # Assume we want all days
    if day_range[0] < 1:
        day_range = (1, day_range[1])
    if day_range[1] > 31:
        day_range = (day_range[0], 31)
    return day_range


def check_hour_range(hour_range):
    """
    :param hour_range: (<int>, <int>) or None: The hour range to check
    :return: (<int>, <int>): A suitable hour range
    """
    if hour_range is None:
        hour_range = (0, 24)  # Assume we want all hours
    if hour_range[0] < 0:
        hour_range = (0, hour_range[1])
    if hour_range[1] > 24:
        hour_range = (hour_range[0], 24)
    return hour_range


def check_gate_range(gate_range, hdw_info):
    """
    :param gate_range: (<int>, <int>) or None: The gate range to check
    :param hdw_info: _HdwInfo: A hardware info object
    :return: (<int>, <int>): A suitable gate range
    """
    if gate_range is None:
        gate_range = (0, 99)  # This is 100 gates
    if gate_range[0] < 0:
        gate_range = (0, gate_range[1])
    if gate_range[1] > hdw_info.gates:
        gate_range = (gate_range[0], hdw_info.gates - 1)
    return gate_range


def check_beam_range(beam_range, hdw_info):
    """
    :param beam_range: (<int>, <int>) or None: The beam range to check
    :param hdw_info: _HdwInfo: A hardware info object
    :return: (<int>, <int>): A suitable beam range
    """
    if beam_range is None:
        beam_range = (0, 15)  # This is 16 beams
    if beam_range[0] < 0:
        beam_range = (0, beam_range[1])
    if beam_range[1] > hdw_info.beams:
        beam_range = (beam_range[0], hdw_info.beams - 1)
    return beam_range


def build_date_epoch(year, month, day, hour):
    """
    :param year: int: the year to consider
    :param month: int: the month to consider
    :param day: int: the day to consider
    :param hour: int: The hour to consider
    :return: time.struct_time, int: The datetime and epoch
    """
    pattern = '%Y.%m.%d %H:%M:%S'

    if hour == 24:
        datetime = str(year) + "." + str(month) + "." + str(day) \
            + " " + "23:59:59"
    else:
        datetime = str(year) + "." + str(month) + "." + str(day) \
            + " " + str(hour) + ":00:00"
    datetime = time.strptime(datetime, pattern)
    epoch = calendar.timegm(datetime)

    return datetime, epoch


def occ_fan(station, year_range, month_range=None, day_range=None, hour_range=None, gate_range=None, beam_range=None):
    """

    Produce a fan plot

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only plots 45 km data
            (a warning will be printed if other spatial resolution data has been stripped from the dataset)
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
    :return: matplotlib.pyplot.figure: The figure.  It can then be modified, added to, printed out, or saved
            in whichever file format is desired.
    """

    show_plot = True  # TODO: remove this

    month_range = check_month_range(month_range)
    day_range = check_day_range(day_range)
    hour_range = check_hour_range(hour_range)

    start_datetime, start_epoch = build_date_epoch(year_range[0], month_range[0], day_range[0], hour_range[0])
    end_datetime, end_epoch = build_date_epoch(year_range[1], month_range[1], day_range[1], hour_range[1])

    if isinstance(station, str):
        hdw_info = pydarn.read_hdw_file(station)  # Get the hardware file, there is lots of good stuff in there
    else:
        raise Exception("Error: Please enter the station as a 3 character string.  e.g. 'rkn'")

    beam_range = check_beam_range(beam_range, hdw_info)
    gate_range = check_gate_range(gate_range, hdw_info)

    radar_lon = hdw_info.geographic.lon
    radar_lat = hdw_info.geographic.lat
    radar_id = hdw_info.stid

    # We need some additional radar info to determine which hemisphere the radar is in
    all_radars_info = SuperDARNRadars()
    additional_radar_info = all_radars_info.radars[radar_id]
    hemisphere = additional_radar_info.hemisphere





    # print("Station id: " + str(hdw_info.stid))
    print("Gates: " + str(hdw_info.gates))
    # print("Beams: " + str(hdw_info.beams))
    # print("Location: " + str(hdw_info.geographic.lon))
    # print(all_radars_info.radars[radar_id])

    # Prepare the figure
    fig = plt.figure(figsize=(5, 5), dpi=300)
    if hemisphere.value == 1:
        # Northern hemisphere: plot above 37 deg lat
        min_lat = 37  # deg
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, 90, min_lat], crs=ccrs.PlateCarree())
    elif hemisphere.value == -1:
        # Southern hemisphere: plot below 50 deg lat
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
    plt.plot([radar_lon, radar_lon], [radar_lat, radar_lat], 'ro', markersize=1.5, transform=ccrs.Geodetic(),
             label=additional_radar_info.name)

    beam_corners_lats, beam_corners_lons = radar_fov(radar_id, coords='geo')  # Get the radar field of view
    fan_shape = beam_corners_lons.shape

    # plot all the beam boundary lines
    for beam_line in range(beam_range[0], beam_range[1] + 2):
        plt.plot(beam_corners_lons[gate_range[0]:gate_range[1] + 2, beam_line],
                 beam_corners_lats[gate_range[0]:gate_range[1] + 2, beam_line],
                 color='black', linewidth=0.1, transform=ccrs.Geodetic())

    # plot all the arcs
    for range_ in range(gate_range[0], gate_range[1] + 2):
        plt.plot(beam_corners_lons[range_, beam_range[0]:beam_range[1] + 2],
                 beam_corners_lats[range_, beam_range[0]:beam_range[1] + 2],
                 color='black', linewidth=0.1, transform=ccrs.Geodetic())

    if show_plot:
        plt.show()
        warnings.warn("This is a test warning", category=Warning)


if __name__ == '__main__':
    thisTuple = (0, 15)
    occ_fan(station="hal", year_range=(2010, 2011))
