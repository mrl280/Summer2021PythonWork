import math
import os

import pandas as pd
import numpy as np
import glob
import h5py

from DataAnalysis.DataReading.RISR_HDF5.wd_beam_num import wd_beam_num

if __name__ == '__main__':
    station = "ran"
    date = "20161012"

    # Loop through all the files for this station/date and pickle all Long Pulse files
    in_dir = "data/" + station + "/" + station + date
    for in_file in glob.iglob(in_dir + "/*.h5"):

        try:
            file = h5py.File(in_file, "r")
        except:
            continue

        # Look at the file type, we are only interested in Long Pulse Files here
        exp_notes = file['/Metadata/Experiment Notes']
        if station == "ran":
            file_type = str(exp_notes[25][0])
        elif station == "ras":
            file_type = str(exp_notes[39][0])
        else:
            raise Exception("Station " + station + " not recognized")
        file_type = file_type[10:len(file_type) - 1].rstrip()  # Strip out everything but the file type name
        if file_type != "Long Pulse":
            print("\n" + in_file + " is a " + file_type + " file, skipping it...")
        else:
            # The file is a long pulse file, we are good to go ahead
            print("\n" + in_file + " is a " + file_type + " file, reading in data...")

            keys = [key for key in file['/Data/Array Layout/'].keys()]  # These will be the list of beams
            try:
                data_time = file['/Data/Array Layout/' + keys[0] + '/timestamps']
            except:
                raise Exception("Error: no keys found, unable to obtain data resolution")
            resolution = int((data_time[1] - data_time[0]) / 60)
            print("Data time resolution: " + str(resolution) + " minute data")
            # Strip the "x.h5" and replace with data resolution and "LongPulse.pkl"
            out_file = in_file[:len(in_file) - 4] + str(resolution) + "min.LongPulse.pkl"

            data = pd.HDFStore(in_file, mode='r')
            # print(type(data))
            print(data.keys())
