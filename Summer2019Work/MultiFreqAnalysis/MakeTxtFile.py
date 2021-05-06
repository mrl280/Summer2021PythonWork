# Michael Luciuk
# Aug 29, 2019

# Last updated Sept 14, 2019

# program to make a text file from an IDL .sav file.  It takes idlsave a long time to read in an IDL .sav file
#   (for good reasons as everything must be converted to Python types)
# The data from the IDL file is used to filter the data so that only relevant data is left in the text file - this
#   relevant data can then be read in quickly by other programs for quick analysis

import idlsave
import csv
import numpy as np


def MakeTxtFile(in_file=None):
    """
    Purpose:
        To read in an IDL.sav file, filter for only the data of interest, and put the important info into a text file,
        The text file can then be read quickly by other analysis programs.
    Pre-conditions:
        :param: in_file: the IDL .sav file that you would like to create a text file from
        If no in_file is given, then Sept 5, 2016 is selected by default.  The given IDL file must be in the
        Data folder of this project
    Post-conditions:
        A text file is created in the same location as the IDL .sav file was read in from.
        The created file has the same name except with the time it was filtered for attached
        Currnet form of output text file is:

        Gate    Ut decimal time     frequency       velocity
        xxx     xxx                 xx              xx
        xxx     xxx                 xx              xx

    Return:
        none
    """
    START_UT = 1.50
    END_UT = 8.00
    PWR_MIN = 3     # minimum power in dB
    N_DECIMAL_PLACES = 5    # number of decimal points to keep when printing to txt file

    if not in_file:
        in_file = "20160905rkn"

    # Read in data from IDL .sav FILE
    print("Reading in IDL.sav file...")
    s = idlsave.read("Data/" + in_file + ".sav")

    # fit data is read in multidimensional arrays of size # of measurements by gates
    # i.e. vel[measurement_idx][gate_idx]

    # prm data is all different
    # Pick out everything that is important from prm
    print("Getting prm data...")
    time = makeTimeArray(s.prm['time'])
    freq = s.prm['tfreq']/1000.
    gate = np.arange(0, 75, 1)

    # duplicate the prm things so they will line up with the flattened things from fit
    NumOfMeasurements = len(s.fit['v'])
    print("# of measurements: " + str(NumOfMeasurements))

    time = np.repeat(time, len(gate))
    freq = np.repeat(freq, len(gate))
    gate = np.tile(gate, NumOfMeasurements)

    # Pick out the important things from fit and flatten them in row-major (C-style) order
    print("Getting and flattening fit data...")
    qflg = flattenFitArrays(s.fit['qflg'])
    gflg = flattenFitArrays(s.fit['gflg'])
    vel = flattenFitArrays(s.fit['v'])
    pwr = flattenFitArrays(s.fit['pwr0'])

    # filter for low quality points or those flagged as ground scatter
    valid = (time >= START_UT) & (time <= END_UT) & (qflg == 1) & (gflg == 0) & (pwr >= PWR_MIN)

    # Write to text file
    print("Writing data to file...")
    OUT_FILE = in_file + "_" + str(START_UT) + "-" + str(END_UT) + "UT"
    full_file_name = "Data/" + OUT_FILE + ".txt"
    with open(full_file_name, "w+") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(gate[valid], np.round(time[valid], N_DECIMAL_PLACES), freq[valid],
                             np.round(vel[valid], N_DECIMAL_PLACES)))

    # remove blank lines from file
    print("Cleaning up file...")
    with open(full_file_name) as filehandle:
        lines = filehandle.readlines()
    with open(full_file_name, 'w') as filehandle:
        lines = filter(lambda x: x.strip(), lines)
        filehandle.writelines(lines)

    print("Program Complete")


def makeTimeArray(time_from_prm):
    # loop thorugh all records in the time array and make a 1D array containing UT time as a decimal
    # return the UT time array
    time = []  # preallocate
    for i in range(len(time_from_prm)):
        # time = np.append(time, time_from_prm[i]['hr'][0] + time_from_prm[i]['mt'][0] / 60.
        #                  + time_from_prm[i]['sc'][0] / 3600. + time_from_prm[i]['us'][0] / 3600. / 1e6)
        time = np.append(time, time_from_prm[i]['hr'][0] + time_from_prm[i]['mt'][0] / 60.
                         + time_from_prm[i]['sc'][0] / 3600.)
    return time


def flattenFitArrays(fit_param):
    # loop through each measurement and add all the gates from that measurement to the array
    # return the row-major (C-style) flattened array
    param = []  # preallocate
    for i in range(len(fit_param)):
        param = np.append(param, fit_param[i])
    return param


if __name__ == "__main__":
    MakeTxtFile()
