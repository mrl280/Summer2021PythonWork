# Michael Luciuk
# Aug 28, 2019

# Create a histograms for the 10/12 ratios

import os
import matplotlib.pyplot as plt
import numpy as np
from getData import getData


def Scat10v12Ratio():
    """
    Purpose:
        Plot a histogram for the 10/12 velocity ratios.
        The percentage of the ratios between 2 and 4 are noted
    Pre-conditions:
        Data files and program to read in data files must exist
        Current data files include
            03_10_12_diff_RKN_15_21_ut_all_data.txt     ; many different days
            03_10_12_diff_RKN_all_data.txt
            06_06_10_12_diff_RKN.txt
            06_06_10_12_diff_RKN_14_23_ut.txt
    Post-conditions:
        Create the plot
        Optionally saves it
    Return:
        none
    """
    #FILE = "03_10_12_diff_RKN_all_data.txt"
    FILE = "03_10_12_diff_RKN_15_21_ut_all_data.txt"
    START_UT = 16.0
    END_UT = 21.0
    MONTH = 3
    DAY = 6
    SAVE_PLOTS = False

    MIN10 = 300     # minimum threshold for 10 MHz data
    MIN12 = 100
    MAX12 = 1000    # maximum threshold for 12 MHz data
    MIN_NUM_RAW = 3  # minimum number of data points required
    OUT_FILE_NAME = "10v12ScatRatio_" + str(START_UT) + "-" + str(END_UT) + " UT"
    FIG_SIZE = (8, 8)
    PLT_RANGE = [-2, 5, 0, 1.2]
    LABEL_SIZE = 16
    TICK_LABEL_SIZE = 17
    HELPER_LINEWIDTH = '1.5'
    PLOT_LETTER_LABELS = ['a', 'b', 'c', 'd']
    PLOT_GATE_LABELS = ['gg: 5-10', 'gg: 11-16', 'gg: 17-22', 'gg: 23-28']
    PLOT_LETTER_LABELS_SIZE = 25
    PLOT_GATE_LABELS_SIZE = 14
    COUNT_SIZE = 16
    PERC_SIZE = 15
    BIN_WDT = 0.2
    INTERVAL_MIN = 2
    INTERVAL_MAX = 4

# < -- SET UP PLOT -- >

    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=FIG_SIZE)

    for i in range(2):
        for j in range(2):
            axs[i, j].set_xlim([PLT_RANGE[0], PLT_RANGE[1]])
            axs[i, j].set_ylim([PLT_RANGE[2], PLT_RANGE[3]])
            axs[i, j].set_xticks(np.arange(-2, 6, 1))
            axs[i, j].set_yticks(np.arange(0.0, 1.4, 0.2))
            axs[i, j].tick_params(axis='both', which='major', labelsize=TICK_LABEL_SIZE)
            axs[i, j].plot([1, 1], [PLT_RANGE[2], PLT_RANGE[3]], 'r-', linewidth=HELPER_LINEWIDTH)
            # axs[i, j].plot([2, 2], [0, 0.5], 'r-', linewidth=HELPER_LINEWIDTH)
            # axs[i, j].plot([4, 4], [0, 0.5], 'r-', linewidth=HELPER_LINEWIDTH)
            axs[i, j].text(PLT_RANGE[0] + 0.89 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                           PLT_RANGE[2] + 0.88 * abs(PLT_RANGE[3] - PLT_RANGE[2]), PLOT_LETTER_LABELS[2 * i + j],
                           fontsize=PLOT_LETTER_LABELS_SIZE)
            axs[i, j].text(PLT_RANGE[0] + 0.075 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                           PLT_RANGE[2] + 0.82 * abs(PLT_RANGE[3] - PLT_RANGE[2]), PLOT_GATE_LABELS[2 * i + j],
                           fontsize=PLOT_GATE_LABELS_SIZE)
            axs[i, j].minorticks_on()

            # axs[i, j].tick_params(axis='y', which='minor', left=False)

    axs[0, 0].set_ylabel('Occurrence', size=LABEL_SIZE)
    axs[1, 1].set_xlabel('10/12 Velocity Ratio', size=LABEL_SIZE)
    axs[1, 0].set_xlabel('10/12 Velocity Ratio', size=LABEL_SIZE)
    axs[1, 0].set_ylabel('Occurrence', size=LABEL_SIZE)

    # fig.suptitle(str(int(START_UT)) + ':' + str(int((START_UT * 60) % 60)).zfill(2) + '-'
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
                     (abs(vel10) >= MIN10) & (abs(vel12) >= MIN12) & (abs(vel12) <= MAX12) & \
                     (numPts10 >= MIN_NUM_RAW) & (numPts12 >= MIN_NUM_RAW)
                     # & (month == MONTH) & (day == DAY)

        ratio = (vel10[valid_here]/vel12[valid_here])

        valid_ratio = (ratio >= -10)
        ratio = ratio[valid_ratio]

        n = len(ratio)
        hist, bin_edges = np.histogram(ratio, bins=np.arange(PLT_RANGE[0], PLT_RANGE[1]+BIN_WDT, BIN_WDT))
        # remove the extra element that has the end of the last bin
        bin_edges = np.delete(bin_edges, len(bin_edges) - 1)

        restricted_ratios = (ratio >= INTERVAL_MIN) & (ratio <= INTERVAL_MAX)
        interval = (bin_edges >= INTERVAL_MIN-BIN_WDT/2) & (bin_edges <= INTERVAL_MAX-BIN_WDT)
        normalized_ratio = hist/max(hist)

        if i < 2:
            axs[0, i].bar(bin_edges, normalized_ratio, width=BIN_WDT, align='edge')
            axs[0, i].bar(bin_edges[interval], normalized_ratio[interval], width=BIN_WDT, align='edge', color='gold')
            axs[0, i].text(PLT_RANGE[0] + 0.075 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                           PLT_RANGE[2] + 0.90 * abs(PLT_RANGE[3] - PLT_RANGE[2]), "n=" + str(n),
                           fontsize=COUNT_SIZE)
            axs[0, i].text(PLT_RANGE[0] + 0.62 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                           PLT_RANGE[2] + 0.3 * abs(PLT_RANGE[3] - PLT_RANGE[2]),
                           str(round(len(ratio[restricted_ratios])/n*100., 1)) + "%",
                           fontsize=PERC_SIZE)
        else:
            axs[1, i - 2].bar(bin_edges, normalized_ratio, width=BIN_WDT, align='edge')
            axs[1, i - 2].bar(bin_edges[interval], normalized_ratio[interval], width=BIN_WDT, align='edge', color='gold')
            axs[1, i - 2].text(PLT_RANGE[0] + 0.075 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                               PLT_RANGE[2] + 0.90 * abs(PLT_RANGE[3] - PLT_RANGE[2]), "n=" +
                               str(n), fontsize=COUNT_SIZE)
            axs[1, i - 2].text(PLT_RANGE[0] + 0.62 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                               PLT_RANGE[2] + 0.3 * abs(PLT_RANGE[3] - PLT_RANGE[2]),
                               str(round(len(ratio[restricted_ratios])/n*100., 1)) + "%",
                               fontsize=PERC_SIZE)

    plt.show()

    if SAVE_PLOTS:
        cur_path = os.path.dirname(__file__)  # where we are
        fig.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.eps', format='eps', bbox_inches='tight')
        fig.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.svg', format='svg', dpi=1200)


if __name__ == "__main__":
    Scat10v12Ratio()
