import datetime as datetime

import aacgmv2


def add_mlt_to_df(beam_corners_aacgm_lons, beam_corners_aacgm_lats, df):
    """

    Add a magnetic local time ('mlt') column to a dataframe

    :param beam_corners_aacgm_lats: numpy.ndarray: Longitudes in aacgm units
    :param beam_corners_aacgm_lons: numpy.ndarray: Latitudes in aacgm units
    :param df: pandas.DataFrame: The dataframe to which you want to add mlt
    :return: pandas.DataFrame: The input dataframe, except now with an 'mlt' column
    """

    aacgm_lons = []
    dates = []

    #  Loop through the dataframe, and build up aacgm_lons and dates
    for i in range(len(df)):
        date = datetime.datetime(df['year'][i], df['month'][i], df['day'][i],
                                 df['hour'][i], df['minute'][i], int(df['second'][i]))
        dates.append(date)

        # TODO: Figure out if you need to get beam_corners_aacgm_lons for every time
        # beam_corners_aacgm_lats, beam_corners_aacgm_lons = \
        #     radar_fov(stid=hdw_info.stid, coords='aacgm', date=date)

        gate_corner = df['slist'][i]
        beam_corner = df['bmnum'][i]

        # estimate the cell with the centroid
        cent_lon, cent_lat = centroid([(beam_corners_aacgm_lons[gate_corner, beam_corner],
                                        beam_corners_aacgm_lats[gate_corner, beam_corner]),
                                       (beam_corners_aacgm_lons[gate_corner + 1, beam_corner],
                                        beam_corners_aacgm_lats[gate_corner + 1, beam_corner]),
                                       (beam_corners_aacgm_lons[gate_corner, beam_corner + 1],
                                        beam_corners_aacgm_lats[gate_corner, beam_corner + 1]),
                                       (beam_corners_aacgm_lons[gate_corner + 1, beam_corner + 1],
                                        beam_corners_aacgm_lats[gate_corner + 1, beam_corner + 1])])
        aacgm_lons.append(cent_lon)

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