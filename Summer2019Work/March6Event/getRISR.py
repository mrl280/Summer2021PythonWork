# Michael Luciuk
# Aug 12, 2019

# Program to read in RISR data from text file

import os
import sys


def getRISR():
    """
    Purpose:
        Read in RISR data from text file
    Pre-conditions:
        none
    Post-conditions:
        Data file must exist
    Return:
        :return: RISR data as a multidimensional list
    """
    cur_path = os.path.dirname(__file__)    # where we are

    # Read in RISR data
    file = open(cur_path + "/Data/risrc_vel_20160306.txt", "r")
    lines = file.readlines()
    start_UT = []  # start time
    end_UT = []    # end time
    vel = []       # velocity
    dvel = []      # error in velocity
    lat = []       # geographical latitude

    # Go through each line and add data to arrays:
    for i, line in enumerate(lines):
        if i == 0:  # first element is a string stating what the column holds
            start_UT.append(line.split()[0])
            end_UT.append(line.split()[1])
            lat.append(line.split()[6])
            vel.append(line.split()[8])
            dvel.append(line.split()[9])
        else:   # the rest is numerical data
            start_UT.append(float(line.split()[0]))
            end_UT.append(float(line.split()[1]))
            lat.append(float(line.split()[6]))
            vel.append(float(line.split()[8]))
            dvel.append(float(line.split()[9]))
    file.close()

    if start_UT[0] != "START_UT" or end_UT[0] != "END_UT" or vel[0] != "VEL" or dvel[0] != "D_VEL" or \
            lat[0] != "LATITUDE":
        print("RISR data not read in correctly")
        sys.exit()
    else:
        return [start_UT, end_UT, lat, vel, dvel]


if __name__ == "__main__":
    RISR_data = getRISR()
