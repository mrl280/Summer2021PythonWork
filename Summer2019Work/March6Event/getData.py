# Michael Luciuk
# Aug 27, 2019

# Read in the data from Koustav's text file

import os


def getData(in_file=None):
    """
    Purpose:
        Read in data from one of Koustav's data files
    Pre-conditions:
        The data file must exist
        At the time of writing, txt files from Koustav are in the /Data/koustav_txt_files folder
        Current data files include
            03_10_12_diff_RKN_15_21_ut_all_data.txt
            03_10_12_diff_RKN_all_data.txt
            06_06_10_12_diff_RKN.txt
            06_06_10_12_diff_RKN_14_23_ut.txt
    Post-conditions:
        none
    Return:
        :return: data as a multidimensional list
    """

    if not in_file:
        in_file = "03_10_12_diff_RKN_15_21_ut_all_data.txt"

    # Read in data from FILE
    cur_path = os.path.dirname(__file__)  # where we are
    file = open(cur_path + "/Data/koustav_txt_files/" + in_file, "r")
    lines = file.readlines()
    month = []      # assuming 2019
    day = []
    gate_min = []
    gate_max = []
    beam_min = []
    beam_max = []

    start_UT_10 = []    # 10 MHz data
    numPts10 = []
    vel10 = []
    stdd_vel10 = []       # standard deviation in 10 MHz data

    start_UT_12 = []    # 12 MHz data
    numPts12 = []
    vel12 = []
    stdd_vel12 = []  # standard deviation in 10 MHz data

    # Go through each line and add data to arrays:
    for line in lines:
        month.append(float(line.split()[0]))
        day.append(float(line.split()[1]))
        gate_min.append(float(line.split()[10]))
        gate_max.append(float(line.split()[11]))
        beam_min.append(float(line.split()[12]))
        beam_max.append(float(line.split()[13]))

        numPts10.append(float(line.split()[2]))
        start_UT_10.append(float(line.split()[3]))
        vel10.append(float(line.split()[4]))
        stdd_vel10.append(float(line.split()[5]))

        numPts12.append(float(line.split()[6]))
        start_UT_12.append(float(line.split()[7]))
        vel12.append(float(line.split()[8]))
        stdd_vel12.append(float(line.split()[9]))

    file.close()

    return [month, day, gate_min, gate_max, beam_min, beam_max,
            start_UT_10, numPts10, vel10, stdd_vel10,
            start_UT_12, numPts12, vel12, stdd_vel12]


if __name__ == "__main__":
    data = getData()
