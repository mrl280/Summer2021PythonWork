# Just a file to learn how to plot with python
from getGateToGlat import glatToGate
from getRISR import getRISR
from getRKN import getRKN
import matplotlib.pyplot as plt
import numpy as np
import math


def main():
    MIN_NUM_RAW = 3
    TEXT_SIZE = 13

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

    # Restrict to a specific time period (19.10-19.19)
    valid_times_RKN = RKN_start_UT == 19.10  # Get the bool array of all elements with start time 19.10
    _RKN_gate = RKN_gate[valid_times_RKN]
    _RKN_vel10 = RKN_med_vel10[valid_times_RKN]
    _RKN_vel12 = RKN_med_vel12[valid_times_RKN]
    _RKN_n10 = RKN_n10[valid_times_RKN]
    _RKN_n12 = RKN_n12[valid_times_RKN]
    _RKN_std_vel10 = RKN_std_vel10[valid_times_RKN]
    _RKN_std_vel12 = RKN_std_vel12[valid_times_RKN]

    # Restrict to range gates <30
    valid_gates_RKN = _RKN_gate <= 30  # Get the bool array of all elements with gates <30
    __RKN_gate = _RKN_gate[valid_gates_RKN]
    __RKN_vel10 = _RKN_vel10[valid_gates_RKN]
    __RKN_vel12 = _RKN_vel12[valid_gates_RKN]
    __RKN_n10 = _RKN_n10[valid_gates_RKN]
    __RKN_n12 = _RKN_n12[valid_gates_RKN]
    __RKN_std_vel10 = _RKN_std_vel10[valid_gates_RKN]
    __RKN_std_vel12 = _RKN_std_vel12[valid_gates_RKN]

    # remove velocity points where there are less than MIN_NUM_RAW data points used in calculation of the mean
    for i in range(len(__RKN_vel10)):
        if __RKN_n10[i] < MIN_NUM_RAW:
            __RKN_vel10[i] = math.nan
        if __RKN_n12[i] < MIN_NUM_RAW:
            __RKN_n12[i] = math.nan

    # Plot velocities as a function of range gate
    plt.rc('font', size=TEXT_SIZE)

    # plot RKN data
    plt.plot(__RKN_gate - 0.1, __RKN_vel10, 'rd', label='10.4 MHz', markersize=5)
    plt.plot(__RKN_gate + 0.1, __RKN_vel12, 'bd', label='12.3 MHz', markersize=5)

    # plot error in RKN data
    for i in range(len(__RKN_vel10)):
        if not math.isnan(__RKN_n10[i]):
            plt.plot([__RKN_gate[i] - 0.1, __RKN_gate[i] - 0.1],
                     [__RKN_vel10[i] + __RKN_std_vel10[i], __RKN_vel10[i] - __RKN_std_vel10[i]], 'r-', linewidth=0.5)
    for i in range(len(__RKN_vel12)):
        if not math.isnan(__RKN_n12[i]):
            plt.plot([__RKN_gate[i] + 0.1, __RKN_gate[i] + 0.1],
                     [__RKN_vel12[i] + __RKN_std_vel12[i], __RKN_vel12[i] - __RKN_std_vel12[i]], 'b-', linewidth=0.5)

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

    _RISR_gate = glatToGate(_RISR_geolat)   # convert geolat to gate

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
    plt.plot(__RISR_gate, __RISR_vel, 'k2', label='RISR', markersize=6)
    plt.plot(__RISR_gate, __RISR_vel, 'k-', linewidth=0.7)

    # Plot error
    for i in range(len(__RISR_gate)):
        plt.plot([_RISR_gate[i], _RISR_gate[i]],
             [_RISR_vel[i] + _RISR_dvel[i], _RISR_vel[i] - _RISR_dvel[i]], 'k-', linewidth=0.5)

    # fix up plot and show
    plt.plot([-1, 31], [0, 0], 'k-', linewidth='0.6')
    plt.xlabel('Range Gate')
    plt.ylabel('Velocity [m/s]')
    plt.title('Velocity Range Profile')
    plt.axis([-0.5, max(__RKN_gate) + 0.5, -1250, 500])
    plt.legend(loc='upper right')
    plt.grid(axis='y', linestyle='--')
    plt.show()


if __name__ == "__main__":
    main()
