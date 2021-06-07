import pathlib
import pandas as pd
import time
import calendar


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


def get_local_dummy_data(station, year, month, day, start_hour_UT, end_hour_UT):
    """
    Get some local dummy data, this program is just a helper for local testing
    On Maxwell data will come from get_data.py, so this program renames some fields to match
    Needs to be called from the same level as the lib folder
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
    in_file = in_dir + "/" + station + str(year) + month + day + ".pkl"
    df = pd.read_pickle(in_file)

    test_start_datetime, test_start_epoch = build_date_epoch(year, int(month), int(day), start_hour_UT)
    test_end_datetime, test_end_epoch = build_date_epoch(year, int(month), int(day), end_hour_UT)

    df = df.loc[(df['epoch'] >= test_start_epoch) & (df['epoch'] <= test_end_epoch)]
    df.reset_index(drop=True, inplace=True)

    # Fix up some of the column names so they match what is coming in from get_data()
    df['slist'] = df['gate']
    df['v'] = df['vel']
    df['w_l'] = df['wdt']
    df['p_l'] = df['pwr']
    df['phi0'] = df['phase']
    df['tfreq'] = df['transFreq']
    df['frang'] = df['firstRang']
    df['rsep'] = df['rangeSep']

    return df
