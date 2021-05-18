import time
import pandas as pd

if __name__ == '__main__':
    """
    Read in pickled RISR Data (data comes from text file)
    """

    station = "ran"
    date = "20140302"

    in_dir = "data/" + station + "/" + station + date + "/"
    in_file = in_dir + station + date + ".pkl"
    df = pd.read_pickle(in_file)
    print(df.head())

    # Grab start and end times
    pattern = '%Y-%m-%d %H:%M:%S'
    date_time = time.strftime(pattern, time.gmtime(df['epoch'].iloc[0]))
    print("\nData start time: " + date_time)
    print("Data end time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[df.shape[0] - 1])))

    print(df['dateTime'].iloc[0])
    print(df['dateTime'].iloc[df.shape[0] - 1])
    print("\n")

    df = df.drop_duplicates(subset=['dateTime'])

    print(df['minute'].iloc[1] - df['minute'].iloc[0])
