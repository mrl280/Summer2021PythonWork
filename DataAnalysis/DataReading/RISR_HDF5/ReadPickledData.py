import time
import bz2
import deprecation

import pandas as pd


@deprecation.deprecated(details="I no longer use RISR data from HDF5 files because they don't include all derived "
                                "parameters.")
def ReadPickleData():
    """
    Read in pickled HDF5 RISR Data
    """

    station = "ras"
    date = "2016923"
    time_interval = '5min'

    in_dir = "data/" + station + "/" + station + date + "/"
    in_file = in_dir + station + date + "." + time_interval + ".LongPulse.pbz2"
    data_stream = bz2.BZ2File(in_file, "rb")
    df = pd.read_pickle(data_stream)
    print(df.head())

    # Grab start and end times
    pattern = '%Y-%m-%d %H:%M:%S'
    date_time = time.strftime(pattern, time.gmtime(df['epoch'].iloc[0]))
    print("\nData start time: " + date_time)
    print("Data end time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[df.shape[0] - 1])))

    year = int(time.strftime(pattern, time.gmtime(df['epoch'].iloc[0]))[:4])
    month = int(time.strftime(pattern, time.gmtime(df['epoch'].iloc[0]))[5:7])
    day = int(time.strftime(pattern, time.gmtime(df['epoch'].iloc[0]))[8:10])

    hour = int(time.strftime(pattern, time.gmtime(df['epoch'].iloc[0]))[11:13])
    minute = int(time.strftime(pattern, time.gmtime(df['epoch'].iloc[0]))[14:16])
    second = int(time.strftime(pattern, time.gmtime(df['epoch'].iloc[0]))[17:19])

    print("Year: " + str(year))
    print("Month: " + str(month))
    print("Day: " + str(day))

    print("Hour: " + str(hour))
    print("Minute: " + str(minute))
    print("Second: " + str(second))

    # print("\nColumns: " + str(list(df.columns)))
    # print(df.size)
    #
    # print(df['bmId'].unique())
    # print(df['wdBmnum'].unique())
    # print("Max range: " + str(df['range'].unique().max()))
    # print("Min range: " + str(df['range'].unique().min()))
    # print("Max magnetic lat: " + str(df['cgmLat'].unique().max()))
    # print("Min magnetic lat: " + str(df['cgmLat'].unique().min()))


if __name__ == '__main__':
    ReadPickleData()
