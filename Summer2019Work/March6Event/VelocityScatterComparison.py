# Michael Luciuk
# Aug 15, 2019

# Program to plot a scatter plot comparing RKN 10 and RKN 12 MHz velocities to RISR velocities (ExB)

from getGateToGlat import glatToGate
from getRISR import getRISR
from getRKN import getRKN
from statistics import mean
import matplotlib.pyplot as plt
import numpy as np
import math
import os


def VelocityScatterComparison():
    """
    Purpose:
        Program to plot a scatter plot comparing RKN 10 MHz and RKN 12 MHz velocities to RISR velocities
    Pre-conditions:
        Data files and the programs to read in these files must exist
    Post-conditions:
        Velocity comparison plot is created
        If SAVE_PLOTS then the plot will be saved in the processed data folder
    Return:
        :return: ax: Plot axes
    """

    SAVE_PLOTS = False  # Save the plot in the Processed data folder
    START_UT = 19.0
    END_UT = 20.5
    MIN_NUM_RAW = 3  # minimum number of data points required for a RKN median to be plotted
    TEXT_SIZE = 19  # size of text in plot
    FIG_SIZE = (8, 7)  # size of the PNG figure created
    PLT_RANGE = [-1400, 400, -1200, 400]  # plot range [xmin, xmax, ymin, ymax]
    OUT_FILE_NAME = "06032016_Velocity_Comparison_" + str(START_UT) + "-" + str(END_UT) + "UT"

    RKN_data = getRKN()  # get RKN data
    # remove the string descriptions and create numpy arrays
    RKN_start_UT = np.asarray(RKN_data[0]) if RKN_data[0].pop(0) == "start" else print("Error: RKN start time data")
    RKN_end_UT = np.asarray(RKN_data[1]) if RKN_data[1].pop(0) == "end" else print("Error getting end time data")
    RKN_gate = np.asarray(RKN_data[2]) if RKN_data[2].pop(0) == "gate" else print("Error: RKN gate data")
    RKN_n10 = np.asarray(RKN_data[3]) if RKN_data[3].pop(0) == "n10" else print("Error: RKN number of 10 MHz points")
    RKN_med_vel10 = np.asarray(RKN_data[4]) if RKN_data[4].pop(0) == "med_vel10" else print("Error: RKN 10MHz vel data")
    # RKN_std_vel10 = np.asarray(RKN_data[5]) if RKN_data[5].pop(0) == "std_vel10" else print("Error: RKN std_vel_10")
    RKN_n12 = np.asarray(RKN_data[6]) if RKN_data[6].pop(0) == "n12" else print("Error: RKN number of 12 MHz points")
    RKN_med_vel12 = np.asarray(RKN_data[7]) if RKN_data[7].pop(0) == "med_vel12" else print("Error: RKN 12MHz vel data")
    # RKN_std_vel12 = np.asarray(RKN_data[8]) if RKN_data[8].pop(0) == "stddev_vel12" else print("Err: RKN std_vel_12")

    # Restrict to a specific time period
    valid_times_RKN = (RKN_start_UT >= START_UT) & (RKN_end_UT <= END_UT)
    _RKN_start_UT = RKN_start_UT[valid_times_RKN]
    _RKN_gate = RKN_gate[valid_times_RKN]
    _RKN_vel10 = RKN_med_vel10[valid_times_RKN]
    _RKN_vel12 = RKN_med_vel12[valid_times_RKN]
    _RKN_n10 = RKN_n10[valid_times_RKN]
    _RKN_n12 = RKN_n12[valid_times_RKN]

    # remove velocity points where there are less than MIN_NUM_RAW data points used in calculation of the median
    for i in range(len(_RKN_vel10)):
        if _RKN_n10[i] < MIN_NUM_RAW:
            _RKN_vel10[i] = math.nan
        if _RKN_n12[i] < MIN_NUM_RAW:
            _RKN_vel12[i] = math.nan

    RISR_data = getRISR()  # get RISR data

    RISR_start_UT = np.asarray(RISR_data[0]) if RISR_data[0].pop(0) == "START_UT" else \
        print("Error: RISR start time data")
    RISR_end_UT = np.asarray(RISR_data[1]) if RISR_data[1].pop(0) == "END_UT" else \
        print("Error: RISR end time data")
    RISR_geolat = np.asarray(RISR_data[2]) if RISR_data[2].pop(0) == "LATITUDE" else \
        print("Error: RISR geolat data")
    RISR_vel = np.asarray(RISR_data[3]) if RISR_data[3].pop(0) == "VEL" else \
        print("Error: RISR velocity data")
    # RISR_dvel = np.asarray(RISR_data[4]) if RISR_data[4].pop(0) == "D_VEL" else \
    #     print("Error: RISR velocity data")

    # Restrict RISR data to a specific time period (START_UT-END_UT)
    valid_times_RISR = (RISR_start_UT >= START_UT) & (RISR_end_UT <= END_UT) & (RISR_start_UT < RISR_end_UT)
    _RISR_start_UT = RISR_start_UT[valid_times_RISR]
    _RISR_geolat = RISR_geolat[valid_times_RISR]
    _RISR_vel = RISR_vel[valid_times_RISR]

    _RISR_gate = glatToGate(_RISR_geolat)  # convert geolat to gate

    # Remove all RKN data that is not within the RISR gate range
    valid_gates_RKN = (_RKN_gate >= math.floor(min(_RISR_gate))) & (_RKN_gate <= math.ceil(max(_RISR_gate)))
    __RKN_start_UT = _RKN_start_UT[valid_gates_RKN]
    __RKN_gate = _RKN_gate[valid_gates_RKN]
    __RKN_vel10 = _RKN_vel10[valid_gates_RKN]
    __RKN_vel12 = _RKN_vel12[valid_gates_RKN]

    # Round all RISR gate data to whole numbers so they can be easily picked out and
    # compared with the corresponding RKN gate
    _RISR_gate_rounded = np.around(_RISR_gate)  # includes gates 9 through 22 (14 gates)

    new_RISR_vel = np.zeros(len(__RKN_gate))    # pre-allocate
    # Loop through the RKN data and make points by averaging the corresponding RISR data
    for i in range(len(__RKN_gate)):
        RISR_idexs_that_match = (_RISR_gate_rounded == __RKN_gate[i]) & (_RISR_start_UT == __RKN_start_UT[i])
        new_RISR_vel[i] = mean(_RISR_vel[RISR_idexs_that_match])    # between 1 and 3 points will be averaged

    # Plot RKN velocities as a function of RISR velocities (RISR velocities should be representative of ExB)
    # Set up plot
    myPlot = plt.figure(figsize=FIG_SIZE)
    plt.rc('font', size=TEXT_SIZE)
    plt.plot([PLT_RANGE[0], PLT_RANGE[1]], [PLT_RANGE[0], PLT_RANGE[1]], 'k-', linewidth='0.6')  # bisector
    plt.plot([PLT_RANGE[0], PLT_RANGE[1]], [0, 0], 'k-', linewidth='0.6')
    plt.plot([0, 0], [PLT_RANGE[2], PLT_RANGE[3]], 'k-', linewidth='0.6')
    plt.xlabel('Velocity (RISR) [m/s]')
    plt.ylabel('Velocity (RKN) [m/s]')
    plt.title('6 March 2016, ' + str(int(START_UT)) + ':' + str(int((START_UT*60) % 60)).zfill(2) + '-'
              + str(int(END_UT)) + ':' + str(int((END_UT*60) % 60)).zfill(2) + ' UT')
    plt.axis(PLT_RANGE)
    plt.grid(linestyle='--')
    ax = plt.gca()
    ax.tick_params(axis='x', rotation=30)
    ax.set_xticks([400, 200, 0, -200, -400, -600, -800, -1000, -1200, -1400])

    # Multiply RISR velocities by 1.2, idk why?
    new_RISR_vel = 1.2 * new_RISR_vel

    # plot the data
    plt.scatter(new_RISR_vel, __RKN_vel10, s=40, facecolors='none', edgecolors='r', label='10.4 MHz', linewidths=1.5)
    plt.scatter(new_RISR_vel, __RKN_vel12, s=40, facecolors='none', edgecolors='b', label='12.3 MHz', linewidths=1.5)

    # Count the number of non-nan data points
    count_10 = 0
    count_12 = 0
    for i in range(len(__RKN_vel10)):
        if not np.isnan(__RKN_vel10[i]):
            count_10 = count_10 + 1
    for i in range(len(__RKN_vel12)):
        if not np.isnan(__RKN_vel12[i]):
            count_12 = count_12 + 1

    print("Count 10 = " + str(count_10))
    print("Count 12 = " + str(count_12))

    # Calculate Cs versus ExB drift curves acording to previous study
    # Fit_99 = 309.0 + 0.0000448 * drift^2
    # Fit_102 = 319.0 + 0.0000674 * drift^2

    # nielsen = 300. + 0.00011 * drift ^ 2
    # koustov = 390 + 0.00015 * drift ^ 2
    # fit_105 = 333.4 + 0.0001045 * drift ^ 2
    # fit_108 = 352.6 + 0.000139 * drift ^ 2
    # fit_111 = 379.6 + 0.000156 * drift ^ 2
    # fit_114 = 410.1 + 0.000141 * drift ^ 2
    # fit_117 = 452.5 + 0.000114 * drift ^ 2
    # fit_120 = 479.8 + 0.000119 * drift ^ 2

    drift = [-300, -600, -800, -1000, -1200, -1400]
    fit_99 = 309.0 + 0.0000448 * np.power(drift, 2)
    fit_102 = 319.0 + 0.0000674 * np.power(drift, 2)
    nielsen = 300. + 0.00011 * np.power(drift, 2)
    koustov = 390 + 0.00015 * np.power(drift, 2)
    fit_105 = 333.4 + 0.0001045 * np.power(drift, 2)
    fit_108 = 352.6 + 0.000139 * np.power(drift, 2)
    fit_111 = 379.6 + 0.000156 * np.power(drift, 2)
    fit_114 = 410.1 + 0.000141 * np.power(drift, 2)
    fit_117 = 452.5 + 0.000114 * np.power(drift, 2)
    fit_120 = 479.8 + 0.000119 * np.power(drift, 2)

    Line_Width = 2
    # Plot the Cs versus ExB drift curves
    plt.plot(drift, -fit_99, 'k-', linewidth=Line_Width)
    plt.plot(drift, -fit_102, 'm-', linewidth=Line_Width)
    plt.plot(drift, -nielsen, 'g-', linewidth=Line_Width)
    # plt.plot(drift, -koustov, 'm-', linewidth=Line_Width)
    # plt.plot(drift, -fit_105, 'c-', linewidth=Line_Width)
    # plt.plot(drift, -fit_108, 'k-', linewidth=Line_Width)
    # plt.plot(drift, -fit_111, 'r-', linewidth=Line_Width)
    # plt.plot(drift, -fit_114, 'g-', linewidth=Line_Width)
    # plt.plot(drift, -fit_117, 'm-', linewidth=Line_Width)
    # plt.plot(drift, -fit_120, 'c-', linewidth=Line_Width)

    # fix up plot and show
    plt.legend(loc='lower right')
    plt.show()

    if SAVE_PLOTS:
        cur_path = os.path.dirname(__file__)  # where we are
        myPlot.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.eps', format='eps', dpi=2400)
        myPlot.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.svg', format='svg')
        # myPlot.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.jpeg', format='jpeg')
        myPlot.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.pdf', format='pdf', dpi=500)
        myPlot.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.jpg', format='jpg', dpi=500)
        # myPlot.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.svg', format='svg', dpi=1200)

    return plt.gca()


if __name__ == "__main__":
    returned_ax = VelocityScatterComparison()
    # plt.show()
