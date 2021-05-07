import time

import pandas as pd

if __name__ == '__main__':
    """
    Read in pickled SuperDARN data
    """

    station = "ran"
    date = "201432"
    start_time = "0001"

    in_dir = "data/" + station + date + "/"
    in_file = in_dir + station + date + ".LongPulse.pkl"
    df = pd.read_pickle(in_file)
    print(df.head())

    # Grab start and end times
    pattern = '%Y-%m-%d %H:%M:%S'
    print("\nData start time: " + time.strftime(pattern, time.gmtime(df['time'].iloc[0])))
    print("Data end time: " + time.strftime(pattern, time.gmtime(df['time'].iloc[df.shape[0] - 1])))

    print("\nColumns: " + str(list(df.columns)))
    print(df.size)

    print(df['wd_bmnum'].unique())
    print("Max range: " + str(df['range'].unique().max()))
    print("Min range: " + str(df['range'].unique().min()))
    print("Max magnetic lat: " + str(df['cgm_lat'].unique().max()))
    print("Min magnetic lat: " + str(df['cgm_lat'].unique().min()))

