# Michael Luciuk
# Aug 12, 2019

# Functions to read in Gate to geographical latitude data from a text file and convert geographical latitude to
# corresponding SuperDARN gate

import os
import sys
import numpy as np
import matplotlib.pyplot as plt


def getGateToGlat():
    """
    Purpose:
        Read in gate to geographical latitude conversion data from text file
    Pre-conditions:
        none
    Post-conditions:
        Data file must exist
    Return:
        :return: a multidimensional list containing the gate and corresponding geographical latitude
    """
    cur_path = os.path.dirname(__file__)    # where we are

    # Read in RKN data
    file = open(cur_path + "/Data/gateTOglat_110km.txt", "r")
    lines = file.readlines()
    gate = []   # gate
    glat = []   # corresponding geographical latitude

    # Go through each line and add data to arrays:
    for i, line in enumerate(lines):
        if i == 0:  # first element is a string stating what the column holds
            gate.append(line.split()[0])
            glat.append(line.split()[3])
        else:   # the rest is numerical data
            gate.append(float(line.split()[0]))
            glat.append(float(line.split()[3]))
    file.close()

    if gate[0] != "GATE" or glat[0] != "GEOLAT":
        print("gate to glat data not read in correctly")
        sys.exit()
    else:
        return [gate, glat]


def glatToGate(glat_arr):
        """
        Purpose:
            Take an array of geolat data and turn it into an array containing the corresponding RKN gate
            :param: glat_arr: A numpy array of geo_lat data
        Pre-conditions:
            none
        Post-conditions:
            none
        Return:
            :return: CorrGateArr: The corresponding gate array as numpy array
        """

        gateToGlat = getGateToGlat()    # get conversion data

        # remove the string descriptions and create numpy arrays
        gate = np.asarray(gateToGlat[0]) if gateToGlat[0].pop(0) == "GATE" else print("Error: Gate data")
        geolat = np.asarray(gateToGlat[1]) if gateToGlat[1].pop(0) == "GEOLAT" else print("Error: Gate data")

        # Do a quadratic fit
        conversion_coeffs = np.polyfit(geolat, gate, 2)
        CorrGateArr = conversion_coeffs[1] * glat_arr + \
                      conversion_coeffs[0] * (np.power(glat_arr, 2)) + conversion_coeffs[2]

        # plt.plot(geolat, gate, 'r-', linewidth=1.0)
        # plt.plot(glat_arr, CorrGateArr, 'ko')
        # plt.show()

        return CorrGateArr


if __name__ == "__main__":
    CorrespondingGateArr = glatToGate(np.asarray([66, 67, 70, 75]))
