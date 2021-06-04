import pathlib
import pandas as pd

from DataAnalysis.EchoOccurrence.lib.build_date_epoch import build_date_epoch


def get_local_dummy_data(station, year, month, day, start_hour_UT, end_hour_UT):
    """
    Get some local dummy data, this program is just a helper for local testing
    On Maxwell data will come from get_data.py, so this program renames some fields to match
    Needs to be called from the same level as the lib folder
    :return:  A dataframe with some dummy local data
    """
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).absolute())
    in_dir = loc_root + "/DataReading/SD/data/" + station + "/" + station + str(year) + str(month) + str(day)
    in_file = in_dir + "/" + station + str(year) + str(month) + str(day) + ".pkl"
    df = pd.read_pickle(in_file)

    test_start_datetime, test_start_epoch = build_date_epoch(year, month, day, start_hour_UT)
    test_end_datetime, test_end_epoch = build_date_epoch(year, month, day, end_hour_UT)

    df = df.loc[(df['epoch'] >= test_start_epoch) & (df['epoch'] <= test_end_epoch)]

    df['v'] = df['vel']
    df['w_l'] = df['wdt']
    df['p_l'] = df['pwr']
    df['phi0'] = df['phase']
    df['tfreq'] = df['transFreq']

    return df
