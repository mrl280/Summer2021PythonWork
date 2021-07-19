import bz2

import statistics as stat
import time
import pandas as pd
import _pickle as cPickle

if __name__ == '__main__':
    """
    Read in pickled RISR Data (data comes from text file)
    """

    station = "ran"
    date = "20101105"

    in_dir = "data/" + station + "/" + station + date + "/"
    in_file = in_dir + station + date + ".1min.pbz2"

    data_stream = bz2.BZ2File(in_file, "rb")
    df = cPickle.load(data_stream)

    print(df.head())

    # Grab start and end times
    pattern = '%Y-%m-%d %H:%M:%S'
    date_time = time.strftime(pattern, time.gmtime(df['epoch'].iloc[0]))
    print("\nData start time: " + date_time)
    print("Data end time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[df.shape[0] - 1])))

    print(df['datetime'].iloc[0])
    print(df['datetime'].iloc[df.shape[0] - 1])
    print("\n")

    df = df.loc[df['wdBmnum'] == 2]
    print(df['elv'].unique())
    print(min(df['aspect'].unique()))
    print(max(df['aspect'].unique()))
    print(stat.mean(df['aspect'].unique()))

    # df = df.drop_duplicates(subset=['dateTime'])
    # print(df['minute'].iloc[1] - df['minute'].iloc[0])
