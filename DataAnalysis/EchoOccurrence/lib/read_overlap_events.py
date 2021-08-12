import bz2
import glob
import pathlib

import deprecation
import pandas as pd
import datetime as datetime


@deprecation.deprecated(details="This program is no longer necessary because the new version of "
                                "find_high_velocity_overlap_events() directly ads in datetime info, filters the data, "
                                "and saves the event dataframe as an easy-view .csv file")
def read_overlap_events():
    """

    Loop thorough all of the pickled overlap event files, and read them in.
    The dataframe can then be sorted and the required information saved to csv.

    See find_high_vel_overlap_events() for more on how the pickled files are produced.

    """

    # We need to go one higher than usual since this is a library file
    loc_root = str((pathlib.Path().parent.absolute().parent.absolute()))
    in_dir = loc_root + "/data/overlap_events"

    for in_file in glob.iglob(in_dir + "/list_of_overlap_events*.pbz2"):
        print("Reading in data from " + in_file)

        data_stream = bz2.BZ2File(in_file, "rb")
        df = pd.read_pickle(data_stream)

        df.sort_values(by=['station1_mean_vel', 'station2_mean_vel'], ascending=[False, False], inplace=True)

        df = df.loc[(df['station1_mean_vel'] >= 300) & (df['station2_mean_vel'] >= 300)]
        df.reset_index(drop=True, inplace=True)

        print("Here are the dataframe's keys:")
        print(df.keys())

        # Compute datetime info
        starting_datetimes = []
        ending_datetimes = []
        for i in range(len(df)):
            starting_datetimes.append(datetime.datetime.utcfromtimestamp(df['start_epoch'].iat[i]))
            ending_datetimes.append(datetime.datetime.utcfromtimestamp(df['end_epoch'].iat[i]))

        df['start'] = starting_datetimes
        df['end'] = ending_datetimes

        print("")
        print(df[['station1', 'station2', 'start_epoch', 'station1_mean_vel', 'station2_mean_vel']])

        df.to_csv(in_file[:-5] + '.csv')


if __name__ == "__main__":

    read_overlap_events()
