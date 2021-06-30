import numpy as np
import pydarn
import pytz

import datetime as datetime

from pydarn import radar_fov, SuperDARNRadars
from timezonefinder import TimezoneFinder

from .add_mlt_to_df import add_mlt_to_df, centroid
from .get_data_handler import get_data_handler


def add_decimal_hour_to_df(df, time_units, stid, date_time_est):
    """
    Add decimal time (hours) column to the dataframe.
    Decimal time is added in the given time_units, and the column is named the same as the given time_units.
    Here is an example of how one might use this function:

        time_units = 'mlt'
        df = add_decimal_hour_to_df(df, time_units, 65, datetime.datetime.now())
        print(df[time_units])  # The added column has the same name as the provided time_units

    :param df: pandas.DataFrame:
            The dataframe to which you want to add some type of decimal hour
    :param time_units: str:
            Decimal time will be added in these units.
            Column will be named with this string
                'ut' for universal time
                'mlt' for magnetic local time
                'lt' for local time (based on longitude)
                'lst' for local standard time (based on time zones).
    :param stid: int:
            The radar station id.
    :param date_time_est: datetime.datetime:
            A datetime object to use for the magnetic field estimate.
    :return: pandas.DataFrame: The input dataframe, except now with a <time_units> columns containing decimal hour time
            in that time_unit
    """

    time_units = time_units.lower()
    if time_units != "ut" and time_units != "mlt" and time_units != "lt" and time_units != "lst":
        raise Exception("add_decimal_hour_to_df(): time_units not recognized.")

    print(" Computing " + time_units.upper() + "s...")

    if time_units != "mlt":
        # All time units except mlt require ut time

        ut = []
        for i in range(len(df)):
            datetime_obj = df['datetime'].iat[i]
            ut.append(datetime_obj.hour + datetime_obj.minute / 60 + datetime_obj.second / 3600)

        if time_units == "ut":
            df['ut'] = np.asarray(ut)

        elif time_units == "lst":
            # Local standard time depends on time zone
            tf = TimezoneFinder()
            cell_corners_lats, cell_corners_lons = radar_fov(stid=stid, coords='geo', date=date_time_est)
            fan_shape = cell_corners_lons.shape
            cell_time_offsets = np.empty(shape=(fan_shape[0] - 1, fan_shape[1] - 1))

            for gate_corner in range(fan_shape[0] - 1):
                for beam_corner in range(fan_shape[1] - 1):
                    cent_lon, cent_lat = centroid([(cell_corners_lons[gate_corner, beam_corner],
                                                    cell_corners_lats[gate_corner, beam_corner]),
                                                   (cell_corners_lons[gate_corner + 1, beam_corner],
                                                    cell_corners_lats[gate_corner + 1, beam_corner]),
                                                   (cell_corners_lons[gate_corner, beam_corner + 1],
                                                    cell_corners_lats[gate_corner, beam_corner + 1]),
                                                   (cell_corners_lons[gate_corner + 1, beam_corner + 1],
                                                    cell_corners_lats[gate_corner + 1, beam_corner + 1])])

                    # Find the time zone offset for this cell, it just depends on the timezone offset
                    cell_timezone_str = tf.timezone_at(lng=cent_lon, lat=cent_lat)
                    cell_timezone = pytz.timezone(cell_timezone_str)
                    dt_now = datetime.datetime.now()
                    time_delta = cell_timezone.utcoffset(dt_now)
                    cell_time_offsets[gate_corner, beam_corner] = time_delta.days * 24 + time_delta.seconds / 3600

            lst = []
            for i in range(len(df)):
                gate = df['slist'].iat[i]
                beam = df['bmnum'].iat[i]

                lst_time = ut[i] + cell_time_offsets[gate, beam]

                # Always keep decimal time in the range 0 to 24
                if lst_time > 24:
                    lst.append(lst_time - 24)
                elif lst_time < 0:
                    lst.append(lst_time + 24)
                else:
                    lst.append(lst_time)

            df['lst'] = np.asarray(lst)

        else:

            # Go ahead and compute lt; local time depends on longitude
            cell_corners_lats, cell_corners_lons = radar_fov(stid=stid, coords='geo', date=date_time_est)
            fan_shape = cell_corners_lons.shape

            cell_time_offsets = np.empty(shape=(fan_shape[0] - 1, fan_shape[1] - 1))

            for gate_corner in range(fan_shape[0] - 1):
                for beam_corner in range(fan_shape[1] - 1):
                    cent_lon, _ = centroid([(cell_corners_lons[gate_corner, beam_corner],
                                             cell_corners_lats[gate_corner, beam_corner]),
                                            (cell_corners_lons[gate_corner + 1, beam_corner],
                                             cell_corners_lats[gate_corner + 1, beam_corner]),
                                            (cell_corners_lons[gate_corner, beam_corner + 1],
                                             cell_corners_lats[gate_corner, beam_corner + 1]),
                                            (cell_corners_lons[gate_corner + 1, beam_corner + 1],
                                             cell_corners_lats[gate_corner + 1, beam_corner + 1])])

                    # Find the time zone offset for this cell
                    # There are 15 degrees in every hour
                    if cent_lon <= 180:
                        # We move East of the Prime Meridian, use positive time offsets
                        time_offset = cent_lon / 15
                    else:
                        # We move West of the Prime Meridian, use negative time offsets
                        cent_lon_deg_W = 360 - cent_lon
                        time_offset = -1 * cent_lon_deg_W / 15

                    cell_time_offsets[gate_corner, beam_corner] = time_offset

            lt = []
            for i in range(len(df)):
                gate = df['slist'].iat[i]
                beam = df['bmnum'].iat[i]

                lt_time = ut[i] + cell_time_offsets[gate, beam]

                # Always keep decimal time in the range 0 to 24
                if lt_time > 24:
                    lt.append(lt_time - 24)
                elif lt_time < 0:
                    lt.append(lt_time + 24)
                else:
                    lt.append(lt_time)

            df['lt'] = np.asarray(lt)

    else:
        # Go ahead and compute mlt
        cell_corners_aacgm_lats, cell_corners_aacgm_lons = radar_fov(stid=stid, coords='aacgm', date=date_time_est)
        df = add_mlt_to_df(cell_corners_aacgm_lons=cell_corners_aacgm_lons,
                           cell_corners_aacgm_lats=cell_corners_aacgm_lats, df=df)

    return df


if __name__ == "__main__":
    """ Testing """

    local_testing = True
    station = "rkn"
    year = 2011
    date_time_est = datetime.datetime.now()

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid

    df = get_data_handler(station, year_range=(year, year), month_range=None, day_range=None,
                          gate_range=(10, 30), beam_range=(6, 8), freq_range=(11, 13), occ_data=True,
                          local_testing=local_testing)

    df = df.head(n=25)

    print("Keys before adding any decimal times:")
    print(df.keys())

    print("Adding ut..")
    df = add_decimal_hour_to_df(df=df, time_units='ut', stid=radar_id, date_time_est=date_time_est)

    print("Adding mlt...")
    df = add_decimal_hour_to_df(df=df, time_units='mlt', stid=radar_id, date_time_est=date_time_est)

    print("Adding lt...")
    df = add_decimal_hour_to_df(df=df, time_units='lt', stid=radar_id, date_time_est=date_time_est)

    print("Adding lst...")
    df = add_decimal_hour_to_df(df=df, time_units='lst', stid=radar_id, date_time_est=date_time_est)

    print(df[['slist', 'bmnum', 'datetime', 'ut', 'mlt', 'lt', 'lst']])
