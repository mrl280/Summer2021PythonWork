# Michael Luciuk
# Aug 22, 2019

# Try to combine Velocity Range Profile and Velocity Scatter Comparison

import numpy as np
import matplotlib.pyplot as plt
import math
import os
from statistics import mean
from getGateToGlat import glatToGate
from getRISR import getRISR
from getRKN import getRKN
from matplotlib.gridspec import GridSpec


def VelocityScatterComparison():
    """
    Purpose:
        Plot a velocity range profile and a velocity scatter comparision together in one figure
        Does not plot morphological groupings
    Pre-conditions:
        Data files and the programs to read in these files must exist
    Post-conditions:
        Figure is created with 2 subplots:
        - Velocity Range Profile (RKN 10 and 12 MHz data as a function of range gate)
        - Velocity Scatter Comparision (RKN 10 and 12 MHz velocity as a function of RISR_HDF5 velocity (ExB))

        Optionally creates an eps file of the figure
    Return:
        none
    """

    START_UT = 19.10
    END_UT = 19.19
    MIN_NUM_RAW = 3  # minimum number of data points required for a RKN median to be plotted
    SAVE_PLOTS = False
    RKN_PT_OFFSET = 0.1  # how offset from the actual position RKN data points are, done so the points are not right \
    # on top of each other
    OUT_FILE_NAME = "06032016_VelRng_VelScat_" + str(START_UT) + "-" + str(END_UT) + "UT"
    FIG_SIZE = (10, 12)  # size of the figure created
    VEL_RNG_PLT_RANGE = [-0.5, 30.5, -1200, 400]
    VEL_SCAT_PLT_RANGE = [-1200, 400, -1200, 400]
    LABEL_SIZE = 25
    TICK_LABEL_SIZE = 20

    fig = plt.figure(figsize=FIG_SIZE)
    gs = GridSpec(50, 50)  # 50 rows, 50 columns
    fig.suptitle('6 March 2016, ' + str(int(START_UT)) + ':' + str(int((START_UT * 60) % 60)).zfill(2) + '-'
                 + str(int(END_UT)) + ':' + str(int((END_UT * 60) % 60)).zfill(2) + 'UT', fontsize=35)

    VelRngAx = fig.add_subplot(gs[0:20, 2:48])
    VelScatAx = fig.add_subplot(gs[26:45, 12:37])

    # get RKN data
    RKN_data = getRKN()
    # remove the string descriptions and create numpy arrays
    RKN_start_UT = np.asarray(RKN_data[0]) if RKN_data[0].pop(0) == "start" else print("Error: RKN start time data")
    RKN_end_UT = np.asarray(RKN_data[1]) if RKN_data[1].pop(0) == "end" else print("Error getting end time data")
    RKN_gate = np.asarray(RKN_data[2]) if RKN_data[2].pop(0) == "gate" else print("Error: RKN gate data")
    RKN_n10 = np.asarray(RKN_data[3]) if RKN_data[3].pop(0) == "n10" else print("Error: RKN number of 10 MHz points")
    RKN_med_vel10 = np.asarray(RKN_data[4]) if RKN_data[4].pop(0) == "med_vel10" else print("Error: RKN 10MHz vel data")
    RKN_std_vel10 = np.asarray(RKN_data[5]) if RKN_data[5].pop(0) == "std_vel10" else print("Error: RKN std_vel_10")
    RKN_n12 = np.asarray(RKN_data[6]) if RKN_data[6].pop(0) == "n12" else print("Error: RKN number of 12 MHz points")
    RKN_med_vel12 = np.asarray(RKN_data[7]) if RKN_data[7].pop(0) == "med_vel12" else print("Error: RKN 12MHz vel data")
    RKN_std_vel12 = np.asarray(RKN_data[8]) if RKN_data[8].pop(0) == "stddev_vel12" else print("Error: RKN std_vel_12")

    # Restrict to a specific time period and gate range
    valid_RKN = (RKN_start_UT == START_UT) & (RKN_end_UT <= END_UT) & (RKN_gate <= 30)
    _RKN_start_UT = RKN_start_UT[valid_RKN]
    _RKN_gate = RKN_gate[valid_RKN]
    _RKN_vel10 = RKN_med_vel10[valid_RKN]
    _RKN_vel12 = RKN_med_vel12[valid_RKN]
    _RKN_std_vel10 = RKN_std_vel10[valid_RKN]
    _RKN_std_vel12 = RKN_std_vel12[valid_RKN]
    _RKN_n10 = RKN_n10[valid_RKN]
    _RKN_n12 = RKN_n12[valid_RKN]

    # remove velocity points where there are less than MIN_NUM_RAW data points used in calculation of the mean
    # also remove any group value that may have been assigned here if necessary
    for i in range(len(_RKN_vel10)):
        if _RKN_n10[i] < MIN_NUM_RAW:
            _RKN_vel10[i] = math.nan
        if _RKN_n12[i] < MIN_NUM_RAW:
            _RKN_vel12[i] = math.nan

    # Set up velocity as a function of range gate plot (on top)
    VelRngAx.set_xlabel('Range Gate', size=LABEL_SIZE)
    VelRngAx.set_ylabel('Velocity [m/s]', size=LABEL_SIZE)
    VelRngAx.minorticks_on()
    VelRngAx.set_yticks([400, 200, 0, -200, -400, -600, -800, -1000, -1200])
    VelRngAx.tick_params(axis='y', which='minor', left=False)
    VelRngAx.tick_params(axis='both', which='major', labelsize=TICK_LABEL_SIZE)
    VelRngAx.set_xlim(left=VEL_RNG_PLT_RANGE[0], right=VEL_RNG_PLT_RANGE[1])
    VelRngAx.set_ylim(bottom=VEL_RNG_PLT_RANGE[2], top=VEL_RNG_PLT_RANGE[3])
    VelRngAx.grid(axis='y', which='major', linestyle='--')
    VelRngAx.plot([-1, 31], [0, 0], 'k-', linewidth='0.6')
    VelRngAx.text(VEL_RNG_PLT_RANGE[0] + 0.94 * abs(VEL_RNG_PLT_RANGE[1] - VEL_RNG_PLT_RANGE[0]),
                  VEL_RNG_PLT_RANGE[2] + 0.04 * abs(VEL_RNG_PLT_RANGE[3] - VEL_RNG_PLT_RANGE[2]), 'a', fontsize=37)

    # plot RKN data, 10 MHz slightly offset to the left, and 12 MHz slightly offset to the right
    VelRngAx.plot(_RKN_gate - RKN_PT_OFFSET, _RKN_vel10, 'rd', label='10.4 MHz', markersize=6.25)
    VelRngAx.plot(_RKN_gate + RKN_PT_OFFSET, _RKN_vel12, 'bd', label='12.3 MHz', markersize=6.25)

    # plot error in RKN data
    for i in range(len(_RKN_vel10)):
        if not math.isnan(_RKN_vel10[i]):
            VelRngAx.plot([_RKN_gate[i] - RKN_PT_OFFSET, _RKN_gate[i] - RKN_PT_OFFSET],
                          [_RKN_vel10[i] + _RKN_std_vel10[i], _RKN_vel10[i] - _RKN_std_vel10[i]], 'r-', linewidth=1.25)
    for i in range(len(_RKN_vel12)):
        if not math.isnan(_RKN_vel12[i]):
            VelRngAx.plot([_RKN_gate[i] + RKN_PT_OFFSET, _RKN_gate[i] + RKN_PT_OFFSET],
                          [_RKN_vel12[i] + _RKN_std_vel12[i], _RKN_vel12[i] - _RKN_std_vel12[i]], 'b-', linewidth=1.25)

    # get RISR_HDF5 data
    RISR_data = getRISR()
    RISR_start_UT = np.asarray(RISR_data[0]) if RISR_data[0].pop(0) == "START_UT" else \
        print("Error: RISR_HDF5 start time data")
    RISR_end_UT = np.asarray(RISR_data[1]) if RISR_data[1].pop(0) == "END_UT" else \
        print("Error: RISR_HDF5 end time data")
    RISR_geolat = np.asarray(RISR_data[2]) if RISR_data[2].pop(0) == "LATITUDE" else \
        print("Error: RISR_HDF5 geolat data")
    RISR_vel = np.asarray(RISR_data[3]) if RISR_data[3].pop(0) == "VEL" else \
        print("Error: RISR_HDF5 velocity data")
    RISR_dvel = np.asarray(RISR_data[4]) if RISR_data[4].pop(0) == "D_VEL" else \
        print("Error: RISR_HDF5 velocity data")

    # Restrict to a specific time period (START_UT-END_UT)
    valid_times_RISR = (RISR_start_UT >= START_UT) & (RISR_end_UT <= END_UT) & (RISR_start_UT < RISR_end_UT)
    # for i in range(len(valid_times_RISR)):
    #     if valid_times_RISR[i]:
    #         print(True)
    # print(valid_times_RISR)
    _RISR_start_UT = RISR_start_UT[valid_times_RISR]
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
    VelRngAx.plot(__RISR_gate, __RISR_vel, 'ks', label='RISR_HDF5', markersize=6)
    VelRngAx.plot(__RISR_gate, __RISR_vel, 'k-', linewidth=1.25)
    for i in range(len(__RISR_gate)):  # Plot error
        VelRngAx.plot([_RISR_gate[i], _RISR_gate[i]],
                      [_RISR_vel[i] + _RISR_dvel[i], _RISR_vel[i] - _RISR_dvel[i]], 'k-', linewidth=1.25)

    # add legend
    VelRngAx.legend(loc='lower left', fontsize=17)

    # Remove all RKN data that is not within the RISR_HDF5 gate range
    valid_gates_RKN = (_RKN_gate >= math.floor(min(_RISR_gate))) & (_RKN_gate <= math.ceil(max(_RISR_gate)))
    __RKN_start_UT = _RKN_start_UT[valid_gates_RKN]
    __RKN_gate = _RKN_gate[valid_gates_RKN]
    __RKN_vel10 = _RKN_vel10[valid_gates_RKN]
    __RKN_vel12 = _RKN_vel12[valid_gates_RKN]

    # Round all RISR_HDF5 gate data to whole numbers so they can be easily picked out and
    # compared with the corresponding RKN gate
    _RISR_gate_rounded = np.around(_RISR_gate)

    new_RISR_vel = np.zeros(len(__RKN_gate))  # pre-allocate
    # Loop through the RKN data and make points by averaging the corresponding RISR_HDF5 data
    for i in range(len(__RKN_gate)):
        RISR_idexs_that_match = (_RISR_gate_rounded == __RKN_gate[i]) & (_RISR_start_UT == __RKN_start_UT[i])
        new_RISR_vel[i] = mean(_RISR_vel[RISR_idexs_that_match])  # between 1 and 3 points will be averaged

    # Set up RKN velocity as a function of RISR_HDF5 velocity plot (on bottom)
    VelScatAx.set_xlim(left=VEL_SCAT_PLT_RANGE[0], right=VEL_SCAT_PLT_RANGE[1])
    VelScatAx.set_ylim(bottom=VEL_SCAT_PLT_RANGE[2], top=VEL_SCAT_PLT_RANGE[3])
    VelScatAx.grid(axis='both', which='major', linestyle='--')
    VelScatAx.set_xlabel('Velocity (RISR_HDF5) [m/s]', size=LABEL_SIZE)
    VelScatAx.set_ylabel('Velocity (RKN) [m/s]', size=LABEL_SIZE)
    VelScatAx.set_xticks([400, 200, 0, -200, -400, -600, -800, -1000, -1200])
    VelScatAx.set_yticks([400, 200, 0, -200, -400, -600, -800, -1000, -1200])
    VelScatAx.plot([VEL_SCAT_PLT_RANGE[0], VEL_SCAT_PLT_RANGE[1]], [VEL_SCAT_PLT_RANGE[0], VEL_SCAT_PLT_RANGE[1]], 'k-',
                   linewidth='0.6')  # bisector
    VelScatAx.plot([VEL_SCAT_PLT_RANGE[0], VEL_SCAT_PLT_RANGE[1]], [0, 0], 'k-', linewidth='0.6')
    VelScatAx.plot([0, 0], [VEL_SCAT_PLT_RANGE[2], VEL_SCAT_PLT_RANGE[3]], 'k-', linewidth='0.6')
    VelScatAx.text(VEL_SCAT_PLT_RANGE[0] + 0.89 * abs(VEL_SCAT_PLT_RANGE[1] - VEL_SCAT_PLT_RANGE[0]),
                   VEL_SCAT_PLT_RANGE[2] + 0.045 * abs(VEL_SCAT_PLT_RANGE[3] - VEL_SCAT_PLT_RANGE[2]), 'b', fontsize=37)
    VelScatAx.tick_params(axis='both', which='major', labelsize=TICK_LABEL_SIZE)
    VelScatAx.tick_params(axis='x', rotation=45)

    # plot the data
    VelScatAx.scatter(new_RISR_vel, __RKN_vel10, s=45, facecolors='none', edgecolors='r', label='10.4 MHz',
                      linewidths=1.5)
    VelScatAx.scatter(new_RISR_vel, __RKN_vel12, s=45, facecolors='none', edgecolors='b', label='12.3 MHz',
                      linewidths=1.5)

    VelScatAx.legend(loc=(0.465, 0.765), fontsize=17)

    if SAVE_PLOTS:
        cur_path = os.path.dirname(__file__)  # where we are
        fig.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.eps', format='eps', bbox_inches='tight')
        fig.savefig(cur_path + '/Processed_data/' + OUT_FILE_NAME + '.svg', format='svg', dpi=1200)

    plt.show()


if __name__ == "__main__":
    VelocityScatterComparison()
