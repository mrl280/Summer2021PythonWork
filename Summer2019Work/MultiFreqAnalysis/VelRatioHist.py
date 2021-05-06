# Michael Luciuk
# Aug 29, 2019

# Create a histograms for the velocity ratios from 15 km mode observations

# Measurements are made every 3 seconds for 2 minutes and then the radar resets and repeats.
# 10 scans at each frequency every 2 minutes
# Since only 5 decimal places are kept in the text files need to be careful when selecting
#   start times and time conditions.

import matplotlib.pyplot as plt
import numpy as np
import os
from getData import getData


def VelRatioHist():
    """
    Purpose:
        Create histograms for the velocity ratios from 15 km mode observations
        Highlight those ratios of an area of interest
    Pre-conditions:
        Text files and program to read in data from those text files must exist
    Post-conditions:
        Histogram plot is plotted
        Optionally, plot is printed to file in Processed data directory
    Return:
        none
    """

    IN_FILE = "20160905rkn_1.5-8.0UT"
    # make sure the correct time interval is actually in the IN_FILE
    START_UT = 1.53333  # important to pick a start time that is equal to the first measurement after the radar resets
    END_UT = 8.0
    SAVE_PLOTS = True
    OUT_FILE_NAME = "VelRatioHist_" + IN_FILE.split('_')[0] + "_" + str(START_UT) + "-" + str(END_UT) + " UT"

    starting_gates = [2, 20, 38, 56]
    starting_gates_fine = [2, 11, 20, 29, 38, 47, 56, 65]

    MIN_NUM_RAW = 3  # minimum number of raw data points required in a grouping to plot the median
    MIN10 = 300  # minimum threshold for 10 MHz data
    MIN12 = 100
    MAX12 = 1000  # maximum threshold for 12 MHz data
    TIME_INTERVAL = 2/60    # time interval in hours (1/60 = 1 minute)
    FIG_SIZE = (8, 8)
    PLT_RANGE = [-2, 5, 0, 1.2]
    LABEL_SIZE = 16
    TICK_LABEL_SIZE = 17
    HELPER_LINEWIDTH = '1.5'
    PLOT_LETTER_LABELS = ['a', 'b', 'c', 'd']
    PLOT_GATE_LABELS = ['gg: 2-19', 'gg: 20-37', 'gg: 38-55', 'gg: 56-73']
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

# < -- Get Data -->

    my_data = getData(IN_FILE)
    gate = np.asarray(my_data[0])
    time = np.asarray(my_data[1])
    freq = np.asarray(my_data[2])
    vel = np.asarray(my_data[3])

    # restrict to time period of intrest
    valid_time = (time >= START_UT) & (time <= END_UT)
    time = time[valid_time]
    freq = freq[valid_time]
    vel = vel[valid_time]
    gate = gate[valid_time]

    # separate out 10 and 12
    hereIs10 = (freq > 9.5) & (freq < 10.5)
    hereIs12 = (freq > 11.5) & (freq < 12.5)

    gate10 = gate[hereIs10]
    vel10 = vel[hereIs10]
    time10 = time[hereIs10]

    gate12 = gate[hereIs12]
    vel12 = vel[hereIs12]
    time12 = time[hereIs12]

    # apply arbitrary restrictions
    valid10 = (vel10 >= MIN10)
    gate10 = gate10[valid10]
    vel10 = vel10[valid10]
    time10 = time10[valid10]
    valid12 = (vel12 >= MIN12) & (vel12 <= MAX12)
    gate12 = gate12[valid12]
    vel12 = vel12[valid12]
    time12 = time12[valid12]

# < -- Analyse data -->

    # combine adjacent gates within TIME_INTERVAL intervals to make median data
    starting_times = np.arange(START_UT, END_UT, TIME_INTERVAL)

    velMedians10 = []   # preallocate
    velMedians12 = []
    start_gate_of_gate_range = []

    num_of_points_tracker = []

    # loop through and make medians
    for start_time in starting_times:
        timeFor10 = (time10 >= start_time) & (time10 < start_time + TIME_INTERVAL)
        timeFor12 = (time12 >= start_time) & (time12 < start_time + TIME_INTERVAL)
        vel10_ = vel10[timeFor10]
        vel12_ = vel12[timeFor12]
        gate10_ = gate10[timeFor10]
        gate12_ = gate12[timeFor12]

        for start_gate in starting_gates_fine:
            gatesFor10 = (gate10_ >= start_gate) & (gate10_ < start_gate + 8)
            gatesFor12 = (gate12_ >= start_gate) & (gate12_ < start_gate + 8)
            vel10__ = vel10_[gatesFor10]
            vel12__ = vel12_[gatesFor12]
            num_of_points_tracker = np.append(num_of_points_tracker, len(vel10__))
            # if len(vel10__) > 90:     # in 2 minute interval with 9 gates
            #     print("Start time: " + str(start_time))
            #     print("End time: " + str(start_time + TIME_INTERVAL))
            #     print("Start gate: " + str(start_gate))
            #     print("End gate: " + str(start_gate+8))
            #     print(len(vel10__))
            #     return
            if len(vel10__) >= MIN_NUM_RAW and len(vel12__) >= MIN_NUM_RAW:
                velMedians10 = np.append(velMedians10, np.median(vel10__))
                velMedians12 = np.append(velMedians12, np.median(vel12__))
                start_gate_of_gate_range = np.append(start_gate_of_gate_range, start_gate)

    ratio = velMedians10/velMedians12
    valid_ratio = (ratio >= -10)
    ratio = ratio[valid_ratio]
    start_gate_of_gate_range = start_gate_of_gate_range[valid_ratio]

    # Now we have the following:
    #   ratio
    #   start_gate_of_gate_range    # The starting gate of the gate range that the corresponding ratio point is for

# < -- Plot the Data -->
    # Loop through each plot and create a histogram of the data
    for i, start_gate in enumerate(starting_gates):
        RatiosHere = (start_gate_of_gate_range >= start_gate) & (start_gate_of_gate_range < start_gate+18)
        ratio_ = ratio[RatiosHere]
        n = len(ratio_)
        hist, bin_edges = np.histogram(ratio_, bins=np.arange(PLT_RANGE[0], PLT_RANGE[1] + BIN_WDT, BIN_WDT))
        # remove the extra element that has the end of the last bin
        bin_edges = np.delete(bin_edges, len(bin_edges) - 1)

        restricted_ratios = (ratio_ >= INTERVAL_MIN) & (ratio_ <= INTERVAL_MAX)
        interval = (bin_edges >= INTERVAL_MIN - BIN_WDT / 2) & (bin_edges <= INTERVAL_MAX - BIN_WDT)
        normalized_ratio = hist / max(hist)

        if i < 2:
            axs[0, i].bar(bin_edges, normalized_ratio, width=BIN_WDT, align='edge')
            axs[0, i].bar(bin_edges[interval], normalized_ratio[interval], width=BIN_WDT, align='edge', color='gold')
            axs[0, i].text(PLT_RANGE[0] + 0.075 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                           PLT_RANGE[2] + 0.90 * abs(PLT_RANGE[3] - PLT_RANGE[2]), "n=" + str(n),
                           fontsize=COUNT_SIZE)
            axs[0, i].text(PLT_RANGE[0] + 0.62 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                           PLT_RANGE[2] + 0.3 * abs(PLT_RANGE[3] - PLT_RANGE[2]),
                           str(round(len(ratio_[restricted_ratios])/n*100., 1)) + "%",
                           fontsize=PERC_SIZE)
        else:
            axs[1, i - 2].bar(bin_edges, normalized_ratio, width=BIN_WDT, align='edge')
            axs[1, i - 2].bar(bin_edges[interval], normalized_ratio[interval], width=BIN_WDT, align='edge', color='gold')
            axs[1, i - 2].text(PLT_RANGE[0] + 0.075 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                               PLT_RANGE[2] + 0.90 * abs(PLT_RANGE[3] - PLT_RANGE[2]), "n=" +
                               str(n), fontsize=COUNT_SIZE)
            axs[1, i - 2].text(PLT_RANGE[0] + 0.62 * abs(PLT_RANGE[1] - PLT_RANGE[0]),
                               PLT_RANGE[2] + 0.3 * abs(PLT_RANGE[3] - PLT_RANGE[2]),
                               str(round(len(ratio_[restricted_ratios])/n*100., 1)) + "%",
                               fontsize=PERC_SIZE)

    plt.show()

    if SAVE_PLOTS:
        cur_path = os.path.dirname(__file__)  # where we are
        fig.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.eps', format='eps', bbox_inches='tight')
        fig.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.svg', format='svg', dpi=1200)


if __name__ == "__main__":
    VelRatioHist()
