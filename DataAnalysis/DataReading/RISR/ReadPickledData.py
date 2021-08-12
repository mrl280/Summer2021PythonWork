import bz2
import time

import numpy as np
import pandas as pd

if __name__ == '__main__':
    """
    Read in pickled RISR Data (data comes from text file)
    """

    station = "ran"
    date = "20091015"
    minute_res = 3

    in_dir = "data/" + station + "/" + station + date + "/"
    in_file = in_dir + station + date + "." + str(minute_res) + "min.pbz2"

    data_stream = bz2.BZ2File(in_file, "rb")
    df = pd.read_pickle(data_stream)

    print(df.head())

    # Grab start and end times
    pattern = '%Y-%m-%d %H:%M:%S'
    date_time = time.strftime(pattern, time.gmtime(df['epoch'].iloc[0]))
    print("\nData start time: " + date_time)
    print("Data end time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[df.shape[0] - 1])))

    print("")
    print("Restricting to world day beam 2...")
    df = df.loc[df['wdBmnum'] == 2]
    print("Here is a list of unique elevations here: " + str(df['elv'].unique()))
    print("Min aspect: " + str(min(df['aspect'].unique())))
    print("Max aspect: " + str(max(df['aspect'].unique())))
    print("Mean aspect: " + str(np.mean(df['aspect'].unique())))

    # df = df.drop_duplicates(subset=['dateTime'])
    # print(df['minute'].iloc[1] - df['minute'].iloc[0])
