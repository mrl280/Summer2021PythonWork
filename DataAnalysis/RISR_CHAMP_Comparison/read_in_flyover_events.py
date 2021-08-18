import pathlib

import pandas as pd


def read_in_flyover_events(in_file):
    """

    Read in the list of flyover events for which there exists both RISR and CHAMP data.

    :param in_file: str:
            The name of the .csv file containing the list of flyovers.
            Should contain the following columns:

                Time:               datetime of flyover
                Distance:           great circle distance between CHAMP and RISR at time of flyover
                champ_data:         is there CHAMP data for this event?  YES or NO
                risr_data:          is there RISR data for this event? YES or NO
                risr_start_day:     Often a single RISR file will contain data for several days. So, this is the day
                                     used in the RISR file name.

    :return event_df: pandas.DataFrame:
            A dataframe with one row for every flyover event
    """

    # Read in the data from csv
    root = str((pathlib.Path().parent.absolute()))
    in_dir = root + "/data"
    in_file = in_dir + "/" + in_file
    print("We are reading in flyover events from: " + in_file)
    df = pd.read_csv(in_file)

    df = df.loc[(df['champ_data'] == "YES") & (df['risr_data'] == "YES")]
    df.reset_index(drop=True, inplace=True)

    return df


if __name__ == "__main__":
    """ """

    in_file = "risr_champ_500km_conjunctions.csv"
    df = read_in_flyover_events(in_file)

    print(df.head())

    print("\nNumber of events: " + str(len(df)))
