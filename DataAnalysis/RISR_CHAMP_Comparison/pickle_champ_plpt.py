import bz2
import glob
import os

import pandas as pd
import dill as pickle

from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch


def pickle_champ_plpt(year, month, day):
    """

    Take PLPT ascii files from CHAMP, turn them into a reduced Pandas DataFrame, and pickle it for later.

    Runs on a single day's data.  ASCII files (with the .dat extension) must be placed in the data/champ directory.

    From CHAMP:
        Product Identifier:     CH-ME-2-PLPT
        Level:                  2 (Level 2 data are calibrated and corrected time series, but more limited on one
                                    dedicated instrument)
        Source:                 PLP
        Comment:                Electron Density and Temperature

    Instructions on how to access CHAMP data can be found here:
     https://isdc.gfz-potsdam.de/champ-isdc/access-to-the-champ-data/

    :param year: int:
            The year of the event to pickle
    :param month: int:
            The month of the event to pickle
    :param day: int
            The day of the event to pickle
    """

    in_dir = "data/champ"
    for in_file in glob.iglob(in_dir + "/*PLPT*.dat"):
        file_name = str(os.path.basename(in_file))

        # TODO: This check is insufficient, we will usually just pickle all data
        # # Check to see if the file is the one we are looking for
        #
        # if str(year) in file_name and str(month) in file_name and str(day) in file_name:
        #     print("We are going to use the following .dat file: " + file_name)
        #     print(file_name)

        # Build a temp file all the unnecessary data stripped out
        temp_file = in_dir + "/temp_" + file_name
        with open(in_file, 'r') as infile, open(temp_file, 'w+') as outfile:
            for line in infile:
                stripped_line = line.strip()
                if not stripped_line or stripped_line[0] == "#":
                    continue  # The line is blank or starts with a #, either way we don't want it

                outfile.write(line.lstrip())  # non-empty line. Write it to output

        # Read in the data from the temp file
        column_names = ["GPS_s", "year", "month", "day",
                        "hour", "minute", "second", "radius_km",
                        "gdlat", "gdlon", "Ne_cm-3", "eTemp_K"]
        df = pd.read_csv(temp_file, names=column_names, delim_whitespace=True, lineterminator='\n')

        print("     Computing datetime and epoch..")
        date_time, epoch = [], []  # Add datetime and epoch to file
        for i in range(len(df)):
            date_time_here, epoch_here = build_datetime_epoch(year=df['year'].iat[i], month=df['month'].iat[i],
                                                              day=df['day'].iat[i], hour=df['hour'].iat[i],
                                                              minute=df['minute'].iat[i], second=df['second'].iat[i])

            date_time.append(date_time_here)
            epoch.append(epoch_here)

        df['datetime'] = date_time
        df['epoch'] = epoch

        print(df.head())

        # Compress and save
        out_file = in_file[:-4] + ".pbz2"
        print("     Pickling as " + out_file + "...")
        with bz2.BZ2File(out_file, "w") as file:
            pickle.dump(df, file)

        # Remove the temporary temp file
        for file in glob.iglob(in_dir + "/temp*.dat"):
            os.remove(file)


if __name__ == '__main__':
    """
    Handler to call pickle_champ_plpt()
    """

    year = 2009
    month = 9
    day = 15

    print("Pickling CHAMP PLPT data from " + str(day) + "/" + str(month) + "/" + str(year) + "...")
    pickle_champ_plpt(year=year, month=month, day=day)
