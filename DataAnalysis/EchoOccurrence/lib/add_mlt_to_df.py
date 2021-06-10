import datetime as datetime
import numpy as np

import aacgmv2


def add_mlt_to_df(cell_corners_aacgm_lons, cell_corners_aacgm_lats, df):
    """

    Add magnetic local time ('mlt'), aacgm longitude ('lon'), and aacgm latitude ('lat') columns to a dataframe

    The first date is used to compute the mlt shift, this shift is then applied to the whole dataframe.
     While less accurate, this approach is much faster than computing mlt for each df row independently.

    :param cell_corners_aacgm_lons: 2d numpy.ndarray: Longitudes of the cell corners in aacgm units
    :param cell_corners_aacgm_lats: 2d  numpy.ndarray: Latitudes of the cell corners in aacgm units
    :param df: pandas.DataFrame: The dataframe to which you want to add mlt
    :return: pandas.DataFrame: The input dataframe, except now with 'mlt', 'lon', and 'lat' columns
    """

    if len(df) <= 0:
        df['lon'], df['lat'], df['mlt'] = [], [], []
        return df

    fan_shape = cell_corners_aacgm_lons.shape

    # Compute cell centroids
    cell_centers_aacgm_lons = np.zeros(shape=(fan_shape[0], fan_shape[1]))
    cell_centers_aacgm_lats = np.zeros(shape=(fan_shape[0], fan_shape[1]))

    for gate_corner in range(fan_shape[0] - 1):
        for beam_corner in range(fan_shape[1] - 1):
            cent_lon, cent_lat = centroid([(cell_corners_aacgm_lons[gate_corner, beam_corner],
                                            cell_corners_aacgm_lats[gate_corner, beam_corner]),
                                           (cell_corners_aacgm_lons[gate_corner + 1, beam_corner],
                                            cell_corners_aacgm_lats[gate_corner + 1, beam_corner]),
                                           (cell_corners_aacgm_lons[gate_corner, beam_corner + 1],
                                            cell_corners_aacgm_lats[gate_corner, beam_corner + 1]),
                                           (cell_corners_aacgm_lons[gate_corner + 1, beam_corner + 1],
                                            cell_corners_aacgm_lats[gate_corner + 1, beam_corner + 1])])
            cell_centers_aacgm_lons[gate_corner, beam_corner] = cent_lon
            cell_centers_aacgm_lats[gate_corner, beam_corner] = cent_lat

    aacgm_lons = []
    aacgm_lats = []
    dates = []

    #  Loop through the dataframe, and build up aacgm_lons and dates
    for i in range(len(df)):
        date = datetime.datetime(df['year'][i], df['month'][i], df['day'][i],
                                 df['hour'][i], df['minute'][i], int(df['second'][i]))
        dates.append(date)

        gate = df['slist'][i]
        beam = df['bmnum'][i]

        aacgm_lons.append(cell_centers_aacgm_lons[gate, beam])
        aacgm_lats.append(cell_centers_aacgm_lats[gate, beam])

    df['lon'] = aacgm_lons
    df['lat'] = aacgm_lats
    df['mlt'] = aacgmv2.convert_mlt(arr=aacgm_lons, dtime=dates, m2a=False)

    return df


def centroid(vertexes):
    """
    Compute the centroid of a polygon
    :param vertexes: a tuple with coordinates of the polygon points
    :return: a tuple with the centroid coordinates
    """
    _x_list = [vertex[0] for vertex in vertexes]
    _y_list = [vertex[1] for vertex in vertexes]
    _len = len(vertexes)
    _x = sum(_x_list) / _len
    _y = sum(_y_list) / _len
    return _x, _y
