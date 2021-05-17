# Michael Luciuk
# Aug 15, 2019

# Plot a velocity scatter comparision of RKN 12 MHz and RISR velocities.  Bin and fit some functions.

from getGateToGlat import glatToGate
from getRISR import getRISR
from getRKN import getRKN
from statistics import mean
import matplotlib.pyplot as plt
import numpy as np
import math
import os


def VelScatCompWFitting():
    """
    Purpose:
        Plot a velocity scatter comparison of RKN 12 MHz and RISR velocities.
        Bin and fit linear, quadratic, rational and logarithmic functions among points that fit certain
        arbitrary criteria.
    Pre-conditions:
        Data files and the programs to read in these files must exist
    Post-conditions:
        Creates a plot of RKN 12 MHz velocity as a function of RISR velocity
        Optionally create an out file with the plot
    Return:
        none
    """

    SAVE_PLOTS = False  # Save the plot in the Processed data folder
    START_UT = 19.0
    END_UT = 20.5
    MIN_NUM_RAW = 3  # minimum number of data points required for a RKN median to be plotted
    TEXT_SIZE = 15  # size of text in plot
    FIG_SIZE = (12, 6)  # size of the PNG figure created
    PLT_RANGE = [-1200, 400, -1200, 400]  # plot range [xmin, xmax, ymin, ymax]
    OUT_FILE_NAME = "06032016_VelScatCompWFitting_" + str(START_UT) + "-" + str(END_UT) + "UT"

# < -- GET RKN DATA -- >

    RKN_data = getRKN()  # get RKN data
    # remove the string descriptions and create numpy arrays
    RKN_start_UT = np.asarray(RKN_data[0]) if RKN_data[0].pop(0) == "start" else print("Error: RKN start time data")
    RKN_end_UT = np.asarray(RKN_data[1]) if RKN_data[1].pop(0) == "end" else print("Error getting end time data")
    RKN_gate = np.asarray(RKN_data[2]) if RKN_data[2].pop(0) == "gate" else print("Error: RKN gate data")
    RKN_n10 = np.asarray(RKN_data[3]) if RKN_data[3].pop(0) == "n10" else print("Error: RKN number of 10 MHz points")
    RKN_med_vel10 = np.asarray(RKN_data[4]) if RKN_data[4].pop(0) == "med_vel10" else print("Error: RKN 10MHz vel data")
    RKN_n12 = np.asarray(RKN_data[6]) if RKN_data[6].pop(0) == "n12" else print("Error: RKN number of 12 MHz points")
    RKN_med_vel12 = np.asarray(RKN_data[7]) if RKN_data[7].pop(0) == "med_vel12" else print("Error: RKN 12MHz vel data")

    # Restrict to a specific time period
    valid_times_RKN = (RKN_start_UT >= START_UT) & (RKN_end_UT <= END_UT)
    _RKN_start_UT = RKN_start_UT[valid_times_RKN]
    _RKN_gate = RKN_gate[valid_times_RKN]
    _RKN_vel10 = RKN_med_vel10[valid_times_RKN]
    _RKN_vel12 = RKN_med_vel12[valid_times_RKN]
    _RKN_n10 = RKN_n10[valid_times_RKN]
    _RKN_n12 = RKN_n12[valid_times_RKN]

    # remove velocity points where there are less than MIN_NUM_RAW data points used in calculation of the mean
    for i in range(len(_RKN_vel10)):
        if _RKN_n10[i] < MIN_NUM_RAW:
            _RKN_vel10[i] = math.nan
        if _RKN_n12[i] < MIN_NUM_RAW:
            _RKN_vel12[i] = math.nan

# < -- GET RISR DATA -- >

    RISR_data = getRISR()  # get RISR data
    RISR_start_UT = np.asarray(RISR_data[0]) if RISR_data[0].pop(0) == "START_UT" else \
        print("Error: RISR start time data")
    RISR_end_UT = np.asarray(RISR_data[1]) if RISR_data[1].pop(0) == "END_UT" else \
        print("Error: RISR end time data")
    RISR_geolat = np.asarray(RISR_data[2]) if RISR_data[2].pop(0) == "LATITUDE" else \
        print("Error: RISR geolat data")
    RISR_vel = np.asarray(RISR_data[3]) if RISR_data[3].pop(0) == "VEL" else \
        print("Error: RISR velocity data")

    # Restrict RISR data to a specific time period (START_UT-END_UT)
    valid_times_RISR = (RISR_start_UT >= START_UT) & (RISR_end_UT <= END_UT) & (RISR_start_UT < RISR_end_UT)
    _RISR_start_UT = RISR_start_UT[valid_times_RISR]
    _RISR_geolat = RISR_geolat[valid_times_RISR]
    _RISR_vel = RISR_vel[valid_times_RISR]

    _RISR_gate = glatToGate(_RISR_geolat)  # convert geolat to gate

# < -- MATCH UP RKN AND RISR DATA -- >

    # Remove all RKN data that is not within the RISR gate range
    valid_gates_RKN = (_RKN_gate >= math.floor(min(_RISR_gate))) & (_RKN_gate <= math.ceil(max(_RISR_gate)))
    __RKN_start_UT = _RKN_start_UT[valid_gates_RKN]
    __RKN_gate = _RKN_gate[valid_gates_RKN]
    __RKN_vel10 = _RKN_vel10[valid_gates_RKN]
    __RKN_vel12 = _RKN_vel12[valid_gates_RKN]

    # Round all RISR gate data to whole numbers so they can be easily picked out and
    # compared with the corresponding RKN gate
    _RISR_gate_rounded = np.around(_RISR_gate)

    new_RISR_vel = np.zeros(len(__RKN_gate))  # pre-allocate
    # Loop through the RKN data and make points by averaging the corresponding RISR data
    for i in range(len(__RKN_gate)):
        RISR_idexs_that_match = (_RISR_gate_rounded == __RKN_gate[i]) & (_RISR_start_UT == __RKN_start_UT[i])
        new_RISR_vel[i] = mean(_RISR_vel[RISR_idexs_that_match])  # between 1 and 3 points will be averaged

# < -- PLOT THE DATA -->

    # Plot RKN velocities as a function of RISR velocities (RISR velocities should be representative of ExB)
    myPlot = plt.figure(figsize=FIG_SIZE)
    plt.rc('font', size=TEXT_SIZE)
    plt.plot([PLT_RANGE[0], PLT_RANGE[1]], [PLT_RANGE[0], PLT_RANGE[1]], 'k-', linewidth='0.6')  # bisector
    plt.plot([PLT_RANGE[0], PLT_RANGE[1]], [0, 0], 'k-', linewidth='0.6')
    plt.plot([0, 0], [PLT_RANGE[2], PLT_RANGE[3]], 'k-', linewidth='0.6')
    plt.xlabel('Velocity (RISR) [m/s]')
    plt.ylabel('12 MHz Velocity (RKN) [m/s]')
    plt.title('6 March 2016, ' + str(int(START_UT)) + ':' + str(int((START_UT * 60) % 60)).zfill(2) + '-'
              + str(int(END_UT)) + ':' + str(int((END_UT * 60) % 60)).zfill(2) + 'UT')
    plt.axis(PLT_RANGE)
    plt.grid(linestyle='--')
    ax = plt.gca()
    ax.tick_params(axis='x', rotation=30)
    ax.set_xticks([400, 200, 0, -200, -400, -600, -800, -1000, -1200])

    # We have to remove all NaN values to make it easier down the line
    not_nan = []
    for i in range(len(__RKN_vel12)):
        if math.isnan(new_RISR_vel[i]) or math.isnan(__RKN_vel12[i]):
            not_nan.append(False)
        else:
            not_nan.append(True)

    new_RISR_vel = new_RISR_vel[not_nan]
    __RKN_vel12 = __RKN_vel12[not_nan]
    # __RKN_vel10 = __RKN_vel10[not_nan]

    # plot the data
    # plt.scatter(new_RISR_vel, __RKN_vel10, s=40, facecolors='none', edgecolors='r', label='10.4 MHz', linewidths=1.5)
    plt.scatter(new_RISR_vel, __RKN_vel12, c='blue', marker="o", label='All 12 MHz pts', s=1)

# < -- FIT SOME OF THE DATA -- >

    # Restrict based on what is currently arbitrary criteria
    # - both have to be of the same sign
    # - points have to be 20% above the bisector
    ForFit = (((__RKN_vel12 < 0) & (new_RISR_vel < 0)) | ((__RKN_vel12 > 0) & (new_RISR_vel > 0))) & \
             (1.3 * __RKN_vel12 > new_RISR_vel)

    ForFit_RKN = __RKN_vel12[ForFit]
    ForFit_RISR = new_RISR_vel[ForFit]
    plt.scatter(ForFit_RISR, ForFit_RKN, s=20, marker="+", c='blue', label='Pts considered', linewidths=0.3)

    # Run binning
    x, y = x_binner(ForFit_RISR, ForFit_RKN, 200, PLT_RANGE[0])
    plt.plot(x, y, c='black', marker='s', label='Binned points', markersize=7)

    # Run rational fit on raw data (y=A+B/x)
    fit = np.polyfit(1 / ForFit_RISR, ForFit_RKN, 1)
    dummy_x = np.arange(PLT_RANGE[0], 0, 1)
    plt.plot(dummy_x, fit[0] / dummy_x + fit[1], label='Rational fit to raw points \ny = '
             + str(round(fit[0], 2)) + '/x ' + str(round(fit[1], 2)), c='m')

    # Run rational fit on the results from binnning (y=A+B/x)
    fit = np.polyfit(1 / x, y, 1)
    plt.plot(dummy_x, fit[0] / dummy_x + fit[1], label='Rational fit to binned results \ny = '
             + str(round(fit[0], 2)) + '/x ' + str(round(fit[1], 2)), c='r')

    # Run quadratic fit on raw data (y=A+Bx+Cx^2)
    fit = np.polyfit(ForFit_RISR, ForFit_RKN, 2)
    plt.plot(dummy_x, fit[0]*np.power(dummy_x, 2) + fit[1]*dummy_x + fit[2],
             label='Quadratic fit to raw points \ny = '
             + str(round(fit[0], 6)) + 'x^2 + ' + str(round(fit[1], 2)) + 'x +' + str(round(fit[2], 2)), c='g')

    # Run quadratic fit on binned data (y=A+Bx+Cx^2)
    fit = np.polyfit(x, y, 2)
    plt.plot(dummy_x, fit[0] * np.power(dummy_x, 2) + fit[1] * dummy_x + fit[2],
             label='Quadratic fit to binned points \ny = '
                   + str(round(fit[0], 6)) + 'x^2 + ' + str(round(fit[1], 2)) + 'x +' + str(round(fit[2], 2)),
                   c='coral')

    # # Run logarithmic fit (y = A + B log(x))
    # fit = np.polyfit(np.log(np.negative(ForFit_RISR)), ForFit_RKN, 1)
    # dummy_x = np.arange(PLT_RANGE[0], 0, 1)
    # plt.plot(dummy_x, fit[0] * np.log(np.negative(dummy_x)) + fit[1], label='Logarithmic fit', c='gold')

    # # Run y = a + bx^2 fit - doesn't do a good job because you need the linear term
    # dummy_y = np.arange(PLT_RANGE[0], 0, 1)
    # fit_coeffs = np.polyfit(np.power(ForFit_RKN, 2), ForFit_RISR, 1)
    # plt.plot(fit_coeffs[0] * np.power(dummy_y, 2) + fit_coeffs[1], dummy_y, label='y = a + bx^2', c='gold')

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()

    if SAVE_PLOTS:
        cur_path = os.path.dirname(__file__)  # where we are
        myPlot.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.eps', format='eps', bbox_inches='tight')
        myPlot.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.svg', format='svg', dpi=1200)


def x_binner(x_data, y_data, interval, start_x):
    """
    Purpose:
        Bin data along the x direction
        Will put a bin every interval from start_x to the end of x_data
    Pre-conditions:
        :param: x_data: numpy array of x data
        :param: y_data: numpy array of y data
        :param: interval: how often to put the data points
        :param: start_x: where to start binning
    Post-conditions:
        none
    Return:
        none
    """
    end_x = math.ceil(max(x_data) / interval) * interval

    # Create x data
    x = np.arange(start_x, end_x, interval)
    y = np.empty(len(x)-1)

    for i in range(len(x)-1):
        x_here = (x_data >= x[i]) & (x_data < x[i+1])
        y[i] = np.median(y_data[x_here])

    x = np.delete(x, len(x)-1)
    x = x + interval/2  # Make sure the x values are in the middle of the interval

    return x, y


if __name__ == "__main__":
    VelScatCompWFitting()
