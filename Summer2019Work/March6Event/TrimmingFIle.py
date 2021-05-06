# Michael Luciuk
# Aug 28, 2019

# Create a histograms for the 10/12 ratios

import os
import matplotlib.pyplot as plt
import numpy as np
from getData import getData


def TrimmingFile():
    """
    Purpose:
        Trim Koustav's Data file to get a list of dates
    Pre-conditions:
        Data files and program to read in data files must exist
        Current data files include
            03_10_12_diff_RKN_15_21_ut_all_data.txt     ; many different days
    Post-conditions:
        Output a file with just the days
    Return:
        none
    """

    FILE = "03_10_12_diff_RKN_15_21_ut_all_data.txt"
    SAVE_PLOTS = False

# < -- GET DATA -- >
    my_data = getData(FILE)
    month = np.asarray(my_data[0])
    day = np.asarray(my_data[1])

    # file = open(“testfile.txt”, ”w”)
    # file.write(“Hello World”)
    # file.close()


if __name__ == "__main__":
    TrimmingFile()
