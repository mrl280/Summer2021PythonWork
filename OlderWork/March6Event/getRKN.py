# Michael Luciuk
# Aug 12, 2019

# Program to read in RKN data from text file

import os
import sys


def getRKN():
    """
    Purpose:
        Read in RKN data from text file
    Pre-conditions:
        none
    Post-conditions:
        Data file must exist
    Return:
        :return: RKN data as a multidimensional list
    """
    cur_path = os.path.dirname(__file__)    # where we are

    # Read in RKN data
    file = open(cur_path + "/Data/RKN_vel_pwr_wdt_wGrouping.txt", "r")
    lines = file.readlines()
    start_UT = []   # start time
    end_UT = []     # end time
    gate = []       # range gate
    n10 = []        # number of points at 10 MHz
    med_vel10 = []  # median velocity at 10 MHz
    std_vel10 = []  # standard deviation in 10 MHz velocities
    n12 = []        # number of points at 12 MHz
    med_vel12 = []  # median velocity at 12 MHz
    stddev_vel12 = []  # standard deviation in 12 MHz velocities

    # Also include the groupings
    group1 = []
    group2 = []
    group3 = []
    group4 = []

    # Go through each line and add data to arrays:
    for i, line in enumerate(lines):
        if i == 0:  # first element is a string stating what the column holds
            start_UT.append(line.split()[0])
            end_UT.append(line.split()[1])
            gate.append(line.split()[2])
            n10.append(line.split()[3])
            med_vel10.append(line.split()[4])
            std_vel10.append(line.split()[5])
            n12.append(line.split()[10])
            med_vel12.append(line.split()[11])
            stddev_vel12.append(line.split()[12])

            group1.append(line.split()[17])
            group2.append(line.split()[18])
            group3.append(line.split()[19])
            group4.append(line.split()[20])
        else:   # the rest is numerical data
            start_UT.append(float(line.split()[0]))
            end_UT.append(float(line.split()[1]))
            gate.append(float(line.split()[2]))
            n10.append(float(line.split()[3]))
            med_vel10.append(float(line.split()[4]))
            std_vel10.append(float(line.split()[5]))
            n12.append(float(line.split()[10]))
            med_vel12.append(float(line.split()[11]))
            stddev_vel12.append(float(line.split()[12]))

            group1.append(float(line.split()[17]))
            group2.append(float(line.split()[18]))
            group3.append(float(line.split()[19]))
            group4.append(float(line.split()[20]))
    file.close()

    # Check to make sure we are reading in what we think we are
    if start_UT[0] != "start" or end_UT[0] != "end" or gate[0] != "gate" or n10[0] != "n10" or \
            med_vel10[0] != "med_vel10" or std_vel10[0] != "std_vel10" or n12[0] != "n12" or med_vel12[0] != \
            "med_vel12" or med_vel12[0] != "med_vel12" or stddev_vel12[0] != "stddev_vel12" or group1[0] != \
            "group1" or group2[0] != "group2" or group3[0] != "group3" or group4[0] != "group4":
        print("RKN data not read in correctly")
        sys.exit()
    else:
        return [start_UT, end_UT, gate, n10, med_vel10, std_vel10, n12, med_vel12, stddev_vel12,
                group1, group2, group3, group4]


if __name__ == "__main__":
    RKN_data = getRKN()
    print(RKN_data[0])
