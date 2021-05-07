import h5py
import glob
import numpy as np
import time

if __name__ == '__main__':
    """
    Try and read in multiple RISR h5 files and look at the timestamps
    Doesn't produce anything, just for looking at the file structure
    """
    # in_dir = "ran2020728"
    in_dir = "ran2017515"
    exp_times = []

    # Loop through every file in the directory and look at the time array
    for filepath in glob.iglob("data/" + in_dir + "/*"):
        try:
            file = h5py.File(filepath, "r")
            data_time = file['/Data/Array Layout/timestamps']
            exp_times += data_time
            print("\n" + filepath)
            print("Here is the start of the time data: " + str(np.asarray(data_time)[0:10]))
            print("Length of time array: " + str(len(data_time)))
            print("First time stamp: " + time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(data_time[0])))
            print("Last time stamp: " + time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(data_time[len(data_time) - 1])))

            # Print out integration time
            data_1Dparam_inttms = file['/Data/Array Layout/1D Parameters/inttms']
            print("Integration time: " + str(data_1Dparam_inttms[()][0]))
            file.close()
        except:
            print("Unable to read: " + filepath)

    # # Look at the time data
    # print(len(exp_times))
    # print(np.asarray(exp_times))
    # print("First time stamp: " + time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(exp_times[0])))
    # print("Last time stamp: " + time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(exp_times[len(exp_times) - 1])))

