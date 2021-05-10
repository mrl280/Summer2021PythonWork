import time

import pandas as pd

if __name__ == '__main__':
    """
    Read in and look at pickled SuperDARN data
    Doesn't produce anything, just for looking at the file structure
    """

    station = "rkn"
    date = "20190318"

    in_dir = "data/" + station + "/" + station + date
    in_file = in_dir + "/" + station + date + ".pkl"
    df = pd.read_pickle(in_file)

    pd.set_option('display.max_columns', None)
    print(df.head())
    # print(df.keys())

    # Grab start and end times
    pattern = '%Y-%m-%d %H:%M:%S'
    print("Data start time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[0])))
    print("Data end time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[df.shape[0] - 1])))

