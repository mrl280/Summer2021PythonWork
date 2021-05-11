import time

import pandas as pd

if __name__ == '__main__':
    """
    Read in pickled SuperDARN data
    """

    station = "ras"
    date = "2016923"
    start_time = "0001"

    in_dir = "data/" + station + "/" + station + date + "/"
    in_file = in_dir + station + date + ".LongPulse.pkl"
    df = pd.read_pickle(in_file)
    print(df.head())

    # Grab start and end times
    pattern = '%Y-%m-%d %H:%M:%S'
    print("\nData start time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[0])))
    print("Data end time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[df.shape[0] - 1])))

    print("\nColumns: " + str(list(df.columns)))
    print(df.size)

    print(df['bmId'].unique())
    print(df['wdBmnum'].unique())
    print("Max range: " + str(df['range'].unique().max()))
    print("Min range: " + str(df['range'].unique().min()))
    print("Max magnetic lat: " + str(df['cgmLat'].unique().max()))
    print("Min magnetic lat: " + str(df['cgmLat'].unique().min()))

