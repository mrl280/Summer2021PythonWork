import bz2
import glob

import numpy as np
import pandas as pd
import datetime as datetime
import _pickle

try:
    from .build_datetime_epoch import build_datetime_epoch
except ImportError:
    # We are performing local testing
    from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch


def pickle_imf(start_year=2010, end_year=2021):
    """
    Take IMF (interplanetary magnetic file) ascii files from OMNI, turn them into a reduced Pandas DataFrame,
     and pickle it for later

    OMNI ascii files (.lst) files need to be uncompressed and in the <EchoOccurrence/data/omni> directory.

    More on OMNI data here: https://omniweb.gsfc.nasa.gov/html/omni_min_data.html
    OMNI data files were created and downloaded from here: https://omniweb.gsfc.nasa.gov/form/omni_min.html

    Here is the file format for the OMNI IMF data:

              FORMAT OF THE SUBSETTED FILE

            ITEMS                      FORMAT

         1 Year                          I4
         2 Day                           I4
         3 Hour                          I3
         4 Minute                        I3
         5 ID for IMF spacecraft         I3
         6 # of points in IMF averages   I4
         7 Timeshift                     I7
         8 RMS, Timeshift                I7
         9 Time btwn observations,sec    I7
        10 Field magnitude average, nT   F8.2
        11 BX, nT (GSE, GSM)             F8.2
        12 BY, nT (GSE)                  F8.2
        13 BZ, nT (GSE)                  F8.2
        14 BY, nT (GSM)                  F8.2
        15 BZ, nT (GSM)                  F8.2
        16 RMS SD B scalar, nT           F8.2
        17 RMS SD field vector, nT       F8.2

    Note that selection of different parameters will necessitate modifying this function.

    :param start_year: int (optional; default is 2010):
            The starting year to consider
    :param end_year: int (optional; default is 2021):
            The final year to consider (inclusive - if you say 2020, you get all data up to Jan 1, 2021)
    """

    column_names = ["year", "day", "hour", "minute",
                    "imf_spacecraft_id", "n_imf_averages", "time_shift", "rms_timeshift", "time_btwn_observations_s",
                    "B_field_avg_nT", "Bx_nT", "By_nT_GSE", "Bz_nT_GSE", "By_nT_GSM", "Bz_nT_GSM",
                    "rms_stddev_B_nT_scalar", "rms_stddev_B_nT_vector"]
    df = pd.DataFrame(columns=column_names)

    # Loop through all the files any add data to the dataframe
    in_dir = "../data/omni"
    for in_file in glob.iglob(in_dir + "/omni_imf*.lst"):
        year = int(in_file[22:26])

        if start_year <= year <= end_year:
            # Then we are good to go ahead and use this file

            print("     Adding in data for " + str(year))
            try:
                # Read in the data for this year, and append it to the dataframe
                temp_df = pd.read_csv(in_file, names=column_names, delim_whitespace=True)
                df = pd.concat([df, temp_df])

            except BaseException as e:
                print("Something went wrong while trying to read " + in_file + ".  Skipping it.")
                print(e)
                pass

        else:
            continue

    # Sometimes the data can spill over the desired year range by a few hours
    df = df.loc[(df['year'] >= start_year) & (df['year'] <= end_year)]
    df.reset_index(drop=True, inplace=True)

    print("     Computing datetime and epoch..")
    date_time, epoch = [], []  # Add datetime and epoch to file
    month, day = [], []  # While we are at it we might as well hold on to this information too
    for i in range(len(df)):
        year = df['year'].iat[i]
        day_of_the_year = df['day'].iat[i]

        # Use a helper date to convert from day of the year to month/day
        helper_date = datetime.datetime.strptime('{} {}'.format(day_of_the_year, year), '%j %Y')
        month.append(helper_date.month)
        day.append(helper_date.day)

        date_time_here, epoch_here = build_datetime_epoch(year=year, month=helper_date.month, day=helper_date.day,
                                                          hour=df['hour'].iat[i], minute=df['minute'].iat[i],
                                                          second=None)

        date_time.append(date_time_here)
        epoch.append(epoch_here)

    df['datetime'] = date_time
    df['epoch'] = epoch
    df['month'] = month
    df['day'] = day

    # OMNI uses 9's as fill values - swap them for nans
    df.replace(to_replace={9999.99: np.nan, 999999: np.nan, 999: np.nan, 99.99: np.nan}, inplace=True)

    # Save to file
    out_file = in_dir + "/omni_imf_" + str(start_year) + "-" + str(end_year) + ".pbz2"
    print("     Pickling as " + out_file + "...")
    with bz2.BZ2File(out_file, "w") as file:
        _pickle.dump(df, file, protocol=4)


if __name__ == '__main__':
    """
    Handler to call pickle_imf()
    """

    start_year = 2013
    end_year = 2021

    print("Pickling OMNI IMF data from " + str(start_year) + " to " + str(end_year) + "...")
    pickle_imf(start_year=start_year, end_year=end_year)
