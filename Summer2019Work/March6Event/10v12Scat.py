# Michael Luciuk
# Aug 27, 2019

# Create a scatter plot of 12 MHz as a function of 10 MHz

import os
import matplotlib.pyplot as plt
import numpy as np
from getData import getData


def Scat10v12():
    """
    Purpose:
        Create a scatter plot of 12 MHz as a function of 10 MHz
        Create plots for the gate ranges specified below
    Pre-conditions:
        Data files and program to read in data files must exist
        Current data files (made by Koustav) include:
            03_10_12_diff_RKN_15_21_ut_all_data.txt
            03_10_12_diff_RKN_all_data.txt
            06_06_10_12_diff_RKN.txt
            06_06_10_12_diff_RKN_14_23_ut.txt
    Post-conditions:
        Create the plot
        Optionally saves it
    Return:
        none
    """

    FILE = "06_06_10_12_diff_RKN_14_23_ut.txt"
    START_UT = 18.0
    END_UT = 21.0
    MONTH = 3
    DAY = 6
    SAVE_PLOTS = True

    MIN_NUM_RAW = 3  # minimum number of data points required
    OUT_FILE_NAME = "10v12Scat_" + str(START_UT) + "-" + str(END_UT) + " UT"
    FIG_SIZE = (8, 8)
    PLT_RANGE = [-1200, 600, -1200, 600]
    LABEL_SIZE = 16
    TICK_LABEL_SIZE = 17
    HELPER_LINEWIDTH = '1.0'
    PLOT_LETTER_LABELS = ['a', 'b', 'c', 'd']
    PLOT_GATE_LABELS = ['gg: 5-10', 'gg: 11-16', 'gg: 17-22', 'gg: 23-28']
    PLOT_LETTER_LABELS_SIZE = 25
    PLOT_GATE_LABELS_SIZE = 14
    COUNT_SIZE = 16
    MARKER_SIZE = 40
    MARKER_WDT = 1

# < -- SET UP PLOT -- >

    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=FIG_SIZE)

    for i in range(2):
        axs[1, i].tick_params(axis='x', rotation=45)
        for j in range(2):
            axs[i, j].set_xlim([PLT_RANGE[0], PLT_RANGE[1]])
            axs[i, j].set_ylim([PLT_RANGE[2], PLT_RANGE[3]])
            axs[i, j].set_xticks([400, 0, -400, -800, -1200])
            axs[i, j].set_yticks([400, 0, -400, -800, -1200])
            axs[i, j].tick_params(axis='both', which='major', labelsize=TICK_LABEL_SIZE)
            axs[i, j].set_aspect('equal')
            axs[i, j].plot([PLT_RANGE[0], PLT_RANGE[1]],
                           [PLT_RANGE[2], PLT_RANGE[3]], 'k-',
                           linewidth=HELPER_LINEWIDTH, color='blue')  # bisector
            axs[i, j].plot([PLT_RANGE[0], PLT_RANGE[1]], [0, 0], 'k--', linewidth=HELPER_LINEWIDTH)
            axs[i, j].plot([0, 0], [PLT_RANGE[2], PLT_RANGE[3]], 'k--', linewidth=HELPER_LINEWIDTH)
            axs[i, j].plot([-1000, -200], [-400, -400], 'r-', linewidth='0.75')
            axs[i, j].plot([-400, -400], [-1000, -200], 'r-', linewidth='0.75')
            axs[i, j].text(PLT_RANGE[0] + 0.89 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                           PLT_RANGE[2] + 0.045 * abs(PLT_RANGE[3] - PLT_RANGE[2]), PLOT_LETTER_LABELS[2 * i + j],
                           fontsize=PLOT_LETTER_LABELS_SIZE)
            axs[i, j].text(PLT_RANGE[0] + 0.072 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                           PLT_RANGE[2] + 0.83 * abs(PLT_RANGE[3] - PLT_RANGE[2]), PLOT_GATE_LABELS[2 * i + j],
                           fontsize=PLOT_GATE_LABELS_SIZE)

    axs[0, 0].set_ylabel('LOS Velocity (12 MHz) [m/s]', size=LABEL_SIZE)
    axs[1, 1].set_xlabel('LOS Velocity (10 MHz) [m/s]', size=LABEL_SIZE)
    axs[1, 0].set_xlabel('LOS Velocity (10 MHz) [m/s]', size=LABEL_SIZE)
    axs[1, 0].set_ylabel('LOS Velocity (12 MHz) [m/s]', size=LABEL_SIZE)

    # fig2.suptitle("March 6, 2016, " + str(int(START_UT)) + ':' + str(int((START_UT * 60) % 60)).zfill(2) + '-'
    #              + str(int(END_UT)) + ':' + str(int((END_UT * 60) % 60)).zfill(2) + 'UT'
    #              + ", Data from: " + FILE, fontsize=10)
    fig.subplots_adjust(wspace=0.15, hspace=0.15)

    # < -- GET DATA -- >
    my_data = getData(FILE)
    month = np.asarray(my_data[0])
    day = np.asarray(my_data[1])
    gate_min = np.asarray(my_data[2])
    gate_max = np.asarray(my_data[3])
    beam_min = np.asarray(my_data[4])
    beam_max = np.asarray(my_data[5])

    UT_10 = np.asarray(my_data[6])
    numPts10 = np.asarray(my_data[7])
    vel10 = np.asarray(my_data[8])
    stdd_vel10 = np.asarray(my_data[9])
    UT_12 = np.asarray(my_data[10])
    numPts12 = np.asarray(my_data[11])
    vel12 = np.asarray(my_data[12])
    stdd_vel12 = np.asarray(my_data[13])

    # Loop through each plot
    starting_gates = [5, 11, 17, 23]
    for i, start_gate in enumerate(starting_gates):
        valid_here = (UT_10 >= START_UT) & (UT_10 <= END_UT) \
                     & (gate_min >= starting_gates[i]) & \
                     (gate_max <= starting_gates[i] + 5) & \
                     (abs(vel10) > 10) & (abs(vel12) > 10) & \
                     (numPts10 >= MIN_NUM_RAW) & (numPts12 >= MIN_NUM_RAW) & \
                     (month == MONTH) & (day == DAY)
        if i < 2:
            axs[0, i].scatter(vel10[valid_here], vel12[valid_here], s=MARKER_SIZE, facecolors='none', edgecolors='k',
                              linewidths=MARKER_WDT)
            axs[0, i].text(PLT_RANGE[0] + 0.075 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                           PLT_RANGE[2] + 0.905 * abs(PLT_RANGE[3] - PLT_RANGE[2]), "n=" + str(len(vel10[valid_here])),
                           fontsize=COUNT_SIZE)
        else:
            axs[1, i - 2].scatter(vel10[valid_here], vel12[valid_here], s=MARKER_SIZE, facecolors='none',
                                  edgecolors='k', linewidths=MARKER_WDT)
            axs[1, i - 2].text(PLT_RANGE[0] + 0.075 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                               PLT_RANGE[2] + 0.905 * abs(PLT_RANGE[3] - PLT_RANGE[2]), "n=" +
                               str(len(vel10[valid_here])), fontsize=COUNT_SIZE)

    plt.show()

    if SAVE_PLOTS:
        cur_path = os.path.dirname(__file__)  # where we are
        fig.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.eps', format='eps', bbox_inches='tight')
        fig.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.svg', format='svg', dpi=1200)


if __name__ == "__main__":
    Scat10v12()
