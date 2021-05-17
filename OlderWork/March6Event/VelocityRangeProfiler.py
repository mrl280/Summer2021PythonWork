# Michael Luciuk
# Aug 15, 2019

# Program to plot RKN and RISR velocities as a function of gate range

from matplotlib import patches
from matplotlib.ticker import MultipleLocator

from getGateToGlat import glatToGate
from getRISR import getRISR
from getRKN import getRKN
import matplotlib.pyplot as plt
import numpy as np
import math
import os


def VelocityRangeProfiler():
    """
    Purpose: Plot RKN and RISR velocities as a function of range gate.  Plots +/- given error in RISR velocities and
        +/- 1 standard deviation in RKN velocities.
        Additionally, boxes are drawn around given groups of RKN velocities, we have the following groups,
            Group 1 (Coral): Low velocity echoes. Two types -> those of the same polarity and those of opposite
            polarities.

            Group 2 (Pink): Higher velocity echoes. Two types -> those of the same polarity and those of opposite
            polarities (Opposite polarity indicates a malfunction).

            Group 3 (Green): High velocity echoes in good agreement.  Both frequencies are detecting ExB.

            Group 4 (Purple): Medium-High velocity echoes in disagreement.  One frequency (usually 10) is detecting
            F while the other (usually 12) is detecting E.
    Pre-conditions:
        Data files and the programs to read in these files must exist
    Post-conditions:
        Velocity Range Profile plot is created
        Optionally an outfile is created with the given plot
    Return:
        :return: ax: Plot axes
    """

    OUT_FILE_NAME = "06032016_Velocity_Range_Profile"
    SAVE_PLOTS = False  # Save the plot in the Processed data folder
    PLOT_GROUPINGS = False  # Plot the morphological groups as described above
    START_UT = 19.10
    END_UT = 19.19
    MIN_NUM_RAW = 3  # minimum number of data points required for a RKN median to be plotted
    TEXT_SIZE = 17  # size of text in plot
    RKN_PT_OFFSET = 0.1  # how offset from the actual position RKN data points are, done so the points are not right \
    # on top of each other
    GROUP_PT_OFFSET = 0.3  # how offset from the actual gate the RKN grouping boxes are
    FIG_SIZE = (10, 6)  # size of the figure created
    PLT_RANGE = [-0.5, 30.5, -1250, 500]    # plot range [xmin, xmax, ymin, ymax]
    OUT_FILE_NAME = "06032016_Velocity_Range_Profile" + str(START_UT) + "-" + str(END_UT) + "UT"

    RKN_data = getRKN()  # get RKN data

    # remove the string descriptions and create numpy arrays
    RKN_start_UT = np.asarray(RKN_data[0]) if RKN_data[0].pop(0) == "start" else print("Error: RKN start time data")
    # RKN_end_UT = np.asarray(RKN_data[1]) if RKN_data[1].pop(0) == "end" else print("Error getting end time data")
    RKN_gate = np.asarray(RKN_data[2]) if RKN_data[2].pop(0) == "gate" else print("Error: RKN gate data")
    RKN_n10 = np.asarray(RKN_data[3]) if RKN_data[3].pop(0) == "n10" else print("Error: RKN number of 10 MHz points")
    RKN_med_vel10 = np.asarray(RKN_data[4]) if RKN_data[4].pop(0) == "med_vel10" else print("Error: RKN 10MHz vel data")
    RKN_std_vel10 = np.asarray(RKN_data[5]) if RKN_data[5].pop(0) == "std_vel10" else print("Error: RKN std_vel_10")
    RKN_n12 = np.asarray(RKN_data[6]) if RKN_data[6].pop(0) == "n12" else print("Error: RKN number of 12 MHz points")
    RKN_med_vel12 = np.asarray(RKN_data[7]) if RKN_data[7].pop(0) == "med_vel12" else print("Error: RKN 12MHz vel data")
    RKN_std_vel12 = np.asarray(RKN_data[8]) if RKN_data[8].pop(0) == "stddev_vel12" else print("Error: RKN std_vel_12")
    RKN_group1 = np.asarray(RKN_data[9]) if RKN_data[9].pop(0) == "group1" else print("Error: group1 data")
    RKN_group2 = np.asarray(RKN_data[10]) if RKN_data[10].pop(0) == "group2" else print("Error: group2 data")
    RKN_group3 = np.asarray(RKN_data[11]) if RKN_data[11].pop(0) == "group3" else print("Error: group3 data")
    RKN_group4 = np.asarray(RKN_data[12]) if RKN_data[12].pop(0) == "group4" else print("Error: group4 data")

    # Restrict to a specific time period and gate range
    valid_RKN = (RKN_start_UT == START_UT) & (RKN_gate <= 30)

    _RKN_gate = RKN_gate[valid_RKN]
    _RKN_vel10 = RKN_med_vel10[valid_RKN]
    _RKN_vel12 = RKN_med_vel12[valid_RKN]
    _RKN_std_vel10 = RKN_std_vel10[valid_RKN]
    _RKN_std_vel12 = RKN_std_vel12[valid_RKN]
    _RKN_n10 = RKN_n10[valid_RKN]
    _RKN_n12 = RKN_n12[valid_RKN]
    _RKN_group1 = RKN_group1[valid_RKN]
    _RKN_group2 = RKN_group2[valid_RKN]
    _RKN_group3 = RKN_group3[valid_RKN]
    _RKN_group4 = RKN_group4[valid_RKN]

    # remove velocity points where there are less than MIN_NUM_RAW data points used in calculation of the mean
    # also remove any group value that may have been assigned here
    for i in range(len(_RKN_vel10)):
        if _RKN_n10[i] < MIN_NUM_RAW:
            _RKN_vel10[i] = math.nan
            _RKN_group1[i] = 0.
            _RKN_group2[i] = 0.
            _RKN_group3[i] = 0.
            _RKN_group4[i] = 0.
        if _RKN_n12[i] < MIN_NUM_RAW:
            _RKN_vel12[i] = math.nan
            _RKN_group1[i] = 0.
            _RKN_group2[i] = 0.
            _RKN_group3[i] = 0.
            _RKN_group4[i] = 0.

    # Plot velocities as a function of range gate
    # Set up plot
    myPlot = plt.figure(figsize=FIG_SIZE)
    plt.rc('font', size=TEXT_SIZE)
    plt.plot([-1, 31], [0, 0], 'k-', linewidth='0.6')
    plt.xlabel('Range Gate')
    plt.ylabel('Velocity [m/s]')
    plt.title('6 March 2016, ' + str(int(START_UT)) + ':' + str(int((START_UT*60) % 60)).zfill(2) + '-'
              + str(int(END_UT)) + ':' + str(int((END_UT*60) % 60)).zfill(2) + 'UT')
    plt.axis(PLT_RANGE)
    ax = plt.gca()
    ax.minorticks_on()
    ax.tick_params(axis='y', which='minor', left=False)

    # plot error in RKN data
    for i in range(len(_RKN_vel10)):
        if not math.isnan(_RKN_vel10[i]):
            plt.plot([_RKN_gate[i] - RKN_PT_OFFSET, _RKN_gate[i] - RKN_PT_OFFSET],
                     [_RKN_vel10[i] + _RKN_std_vel10[i], _RKN_vel10[i] - _RKN_std_vel10[i]], 'r-', linewidth=1.25)
    for i in range(len(_RKN_vel12)):
        if not math.isnan(_RKN_vel12[i]):
            plt.plot([_RKN_gate[i] + RKN_PT_OFFSET, _RKN_gate[i] + RKN_PT_OFFSET],
                     [_RKN_vel12[i] + _RKN_std_vel12[i], _RKN_vel12[i] - _RKN_std_vel12[i]], 'b-', linewidth=1.25)

    if PLOT_GROUPINGS:
        # Look through all groups and add rectangles to the plot
        for group in range(4):
            if group == 0:  # go through each group
                # Low velocity echoes
                # Two types -> those of the same polarity and those of opposite polarity
                group_i = _RKN_group1
                COL = 'coral'
                BASE_HEIGHT = -225
                HEIGHT = 275
            elif group == 1:
                # Higher velocity echoes, both frequencies detecting E region echoes Two types -> those of the same
                # polarity and those of opposite polarities, opposite polarities indicate malfunction
                group_i = _RKN_group2
                COL = 'deeppink'
                BASE_HEIGHT = -525
                HEIGHT = 500
            elif group == 2:
                # High velocity echoes in good agreement, both frequencies detecting ExB
                group_i = _RKN_group3
                COL = 'forestgreen'
                BASE_HEIGHT = -1025
                HEIGHT = 875
            else:
                # High velocity echoes in disagreement
                # One frequency (usually 10) is detecting F while the other (usually 12) is detecting E
                group_i = _RKN_group4
                COL = 'm'
                BASE_HEIGHT = -1025
                HEIGHT = 1250

            i = 0
            while i < len(group_i):
                if group_i[i]:
                    j = i + 1
                    if j == len(group_i):  # we are at the end of the array
                        group1_rect = patches.Rectangle((_RKN_gate[i] - GROUP_PT_OFFSET, BASE_HEIGHT),
                                                        2 * GROUP_PT_OFFSET, HEIGHT,  # One one spot, always 0.6 wide
                                                        linewidth=1, edgecolor=COL, facecolor='none')
                        plt.gca().add_patch(group1_rect)  # have to add right away encase there is another group2 set
                    else:
                        while group_i[j]:
                            j += 1
                            if j == len(group_i):  # we are at the end of the array
                                break
                        group1_rect = patches.Rectangle((_RKN_gate[i] - GROUP_PT_OFFSET, BASE_HEIGHT),
                                                        (j - i - 1 + 2 * GROUP_PT_OFFSET), HEIGHT,
                                                        linewidth=1, edgecolor=COL, facecolor='none')
                        plt.gca().add_patch(group1_rect)  # have to add right away encase there is another group2 set
                        i = j + 1
                i += 1

    # plot RKN data, 10 MHz slightly offset to the left, and 12 MHz slightly offset to the right
    plt.plot(_RKN_gate - RKN_PT_OFFSET, _RKN_vel10, 'rd', label='10.4 MHz', markersize=6.25)
    plt.plot(_RKN_gate + RKN_PT_OFFSET, _RKN_vel12, 'bd', label='12.3 MHz', markersize=6.25)

    RISR_data = getRISR()  # get RISR data

    RISR_start_UT = np.asarray(RISR_data[0]) if RISR_data[0].pop(0) == "START_UT" else \
        print("Error: RISR start time data")
    # RISR_end_UT = np.asarray(RISR_data[1]) if RISR_data[1].pop(0) == "END_UT" else \
    #     print("Error: RISR end time data")
    RISR_geolat = np.asarray(RISR_data[2]) if RISR_data[2].pop(0) == "LATITUDE" else \
        print("Error: RISR geolat data")
    RISR_vel = np.asarray(RISR_data[3]) if RISR_data[3].pop(0) == "VEL" else \
        print("Error: RISR velocity data")
    RISR_dvel = np.asarray(RISR_data[4]) if RISR_data[4].pop(0) == "D_VEL" else \
        print("Error: RISR velocity data")

    # Restrict to a specific time period (19.10-19.19)
    valid_times_RISR = RISR_start_UT == 19.10  # Get the bool array of all elements with start time 19.10
    _RISR_geolat = RISR_geolat[valid_times_RISR]
    _RISR_vel = RISR_vel[valid_times_RISR]
    _RISR_dvel = RISR_dvel[valid_times_RISR]

    _RISR_gate = glatToGate(_RISR_geolat)  # convert geolat to gate

    # since we are connecting with a line we have to remove all NaN values so there are no breaks in the line
    not_nan = []
    for i in range(len(_RISR_vel)):
        if math.isnan(_RISR_vel[i]):
            not_nan.append(False)
        else:
            not_nan.append(True)

    __RISR_gate = _RISR_gate[not_nan]
    __RISR_vel = _RISR_vel[not_nan]
    __RISR_dvel = _RISR_dvel[not_nan]

    # plot the data
    plt.plot(__RISR_gate, __RISR_vel, 'ks', label='RISR', markersize=6)
    plt.plot(__RISR_gate, __RISR_vel, 'k-', linewidth=1.25)

    # Plot error
    for i in range(len(__RISR_gate)):
        plt.plot([_RISR_gate[i], _RISR_gate[i]],
                 [_RISR_vel[i] + _RISR_dvel[i], _RISR_vel[i] - _RISR_dvel[i]], 'k-', linewidth=1.25)

    # fix up plot and show
    plt.legend(loc='lower left')
    plt.grid(axis='y', linestyle='--')

    # Look at size of plot
    # fig_size = plt.rcParams["figure.figsize"]   # default [6.4, 4.8]
    # print(fig_size)

    # plt.show()

    if SAVE_PLOTS:
        cur_path = os.path.dirname(__file__)  # where we are
        myPlot.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.ps', format='ps', orientation='landscape')
        # myPlot.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.svg', format='svg', dpi=1200)

    return plt.gca()


if __name__ == "__main__":
    returned_ax = VelocityRangeProfiler()
    plt.show()
