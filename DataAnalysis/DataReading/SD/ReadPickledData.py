import time

import pandas as pd

if __name__ == '__main__':
    """
    Read in pickled SuperDARN data
    """

    station = "rkn"
    date = "20190318"
    start_time = "00"

    in_dir = "data/" + station + date + "/"
    in_file = in_dir + date + "." + start_time + "." + station + ".pkl"
    df = pd.read_pickle(in_file)
    print(df.head())

    # Grab start and end times
    pattern = '%Y-%m-%d %H:%M:%S'
    print("Data start time: " + time.strftime(pattern, time.gmtime(df['time'].iloc[0])))
    print("Data end time: " + time.strftime(pattern, time.gmtime(df['time'].iloc[df.shape[0] - 1])))

