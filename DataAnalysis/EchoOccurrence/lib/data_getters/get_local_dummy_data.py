import bz2
import pathlib
import time
import calendar

import _pickle as cPickle


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


def get_local_dummy_data(station, year, month, day, start_hour_UT, end_hour_UT, occ_data=False):
    """
    Get some local dummy data, this program is just a helper for local testing
    On Maxwell data will come from get_data.py, so this program renames some fields to match

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year: int: Y
    :param month: int: m
    :param day: int: d
    :param start_hour_UT: int: Start hour
    :param end_hour_UT: int: End hour
    :param occ_data: bool (optional):
            Set this to true if you need echo occurrence data.
            If False, you will get normal data (basically a reduced fitACF datafile),  Default if False.
    :return:  A dataframe with some dummy local data
    """
    if month < 10:
        month = "0" + str(month)
    else:
        month = str(month)
    if day < 10:
        day = "0" + str(day)
    else:
        day = str(day)

    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).absolute())
    in_dir = loc_root + "/DataReading/SD/data/" + station + "/" + station + str(year) + month + day
    if occ_data:
        in_file = in_dir + "/" + station + str(year) + month + day + "_occ.pbz2"
    else:
        in_file = in_dir + "/" + station + str(year) + month + day + ".pbz2"

    try:
        data_stream = bz2.BZ2File(in_file, "rb")
        df = cPickle.load(data_stream)
    except FileNotFoundError as e:
        # We might be calling from one level down, recompute the path and try again
        loc_root = str((((pathlib.Path().parent.absolute()).parent.absolute()).absolute()).parent.absolute())
        in_dir = loc_root + "/DataReading/SD/data/" + station + "/" + station + str(year) + month + day

        if occ_data:
            in_file = in_dir + "/" + station + str(year) + month + day + "_occ.pbz2"
        else:
            in_file = in_dir + "/" + station + str(year) + month + day + ".pbz2"

        data_stream = bz2.BZ2File(in_file, "rb")
        df = cPickle.load(data_stream)

    test_start_datetime, test_start_epoch = build_date_epoch(year, int(month), int(day), start_hour_UT)
    test_end_datetime, test_end_epoch = build_date_epoch(year, int(month), int(day), end_hour_UT)

    df = df.loc[(df['epoch'] >= test_start_epoch) & (df['epoch'] <= test_end_epoch)]
    df.reset_index(drop=True, inplace=True)

    # Fix up some of the column names so they match what is coming in from the occurrence packages' data getters
    if not occ_data:
        df['slist'] = df['gate']
        df['v'] = df['vel']
        df['w_l'] = df['wdt']
        df['p_l'] = df['pwr']
        df['phi0'] = df['phase']
        df['tfreq'] = df['transFreq']
        df['frang'] = df['firstRang']
        df['rsep'] = df['rangeSep']

    return df
