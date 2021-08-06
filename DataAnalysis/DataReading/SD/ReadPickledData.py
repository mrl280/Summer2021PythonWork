import bz2
import time

import numpy as np

import pandas as pd

if __name__ == '__main__':
    """
    Read in and look at pickled SuperDARN data
    Doesn't produce anything, just for looking at the file structure
    """

    station = "rkn"
    date = "20111112"

    in_dir = "data/" + station + "/" + station + date
    in_file = in_dir + "/" + station + date + ".pbz2"

    data_stream = bz2.BZ2File(in_file, "rb")
    df = pd.read_pickle(data_stream)

    # pd.set_option('display.max_columns', None)
    print(df.head())
    df['phaseInDg'] = np.degrees(np.asarray(df['phase']))
    print(df[['phase', 'phaseInDg', 'elv']])

    # Grab start and end times
    pattern = '%Y-%m-%d %H:%M:%S'
    print("Data start time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[0])))
    print("Data end time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[df.shape[0] - 1])))

    print(type(df['year'][0]))
