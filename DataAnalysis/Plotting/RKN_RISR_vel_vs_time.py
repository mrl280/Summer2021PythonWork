import calendar
import time
import os
import pathlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import pandas as pd

if __name__ == '__main__':
    """
    Plot SuperDARN and RISR velocity vs time
    """

    # TODO: Update so it will work with multiple SuperDARN data chunks

    SAVE_PLOT = False
    pattern = '%Y.%m.%d %H:%M:%S'

    year = "2014"
    month = "03"
    day = "02"
    start_time = "20:00:00"  # %H:%M:%S
    end_time = "22:00:00"  # %H:%M:%S

    SD_station = "rkn"
    SD_time_chunk = start_time[0:2]
    SD_beam_range = [5, 5]
    SD_gate_range = [30, 40]

    RISR_station = "ran"
    RISR_wd_beam_range = [2, 2]
    RISR_range = [600, 900]  # Not sure what this is

    # Compute start and end epochs
    start_date_time = year + "." + month + "." + day + " " + start_time
    end_date_time = year + "." + month + "." + day + " " + end_time
    start_epoch = calendar.timegm(time.strptime(start_date_time, pattern))
    end_epoch = calendar.timegm(time.strptime(end_date_time, pattern))

    # Read in SuperDARN data
    loc_root = str((pathlib.Path().parent.absolute()).parent.absolute())
    SD_in_dir = loc_root + "/DataReading/SD/data/" + SD_station + year + month + day + "/"
    SD_in_file = SD_in_dir + year + month + day + "." + SD_time_chunk + "." + SD_station + ".pkl"
    SD_df = pd.read_pickle(SD_in_file)

    # Read in RISR data
    RISR_month = month[1:] if month.startswith('0') else month  # RISR dates don't have the start 0
    RISR_day = day[1:] if day.startswith('0') else day
    RISR_in_dir = loc_root + "/DataReading/RISR/data/" + RISR_station + year + RISR_month + RISR_day + "/"
    RISR_in_file = RISR_in_dir + RISR_station + year + RISR_month + RISR_day + ".LongPulse.pkl"
    RISR_df = pd.read_pickle(RISR_in_file)

    # Filter SuperDARN data
    SD_df = SD_df.loc[(SD_df['time'] >= start_epoch) & (SD_df['time'] <= end_epoch)  # filer times
                      & (SD_df['gate'] >= SD_gate_range[0]) & (SD_df['gate'] <= SD_gate_range[1])  # filter gates
                      & (SD_df['bmnum'] >= SD_beam_range[0]) & (SD_df['bmnum'] <= SD_beam_range[1])]  # filter beams
    SD_df.reset_index(drop=True, inplace=True)

    RISR_df = RISR_df.loc[(RISR_df['time'] >= start_epoch) & (RISR_df['time'] <= end_epoch)  # filer times
                          & (RISR_df['wd_bmnum'] >= RISR_wd_beam_range[0])
                          & (RISR_df['wd_bmnum'] <= RISR_wd_beam_range[1])  # filter beams
                          & (RISR_df['range'] >= RISR_range[0])
                          & (RISR_df['range'] <= RISR_range[1])]  # filter range
    RISR_df.reset_index(drop=True, inplace=True)

    # Recover decimal time from epoch
    SD_decimal_time = []
    for t in range(len(SD_df['time'])):
        hour = str.split(time.strftime(pattern, time.gmtime(SD_df['time'][t])))[1][0:2]
        min = str.split(time.strftime(pattern, time.gmtime(SD_df['time'][t])))[1][3:5]
        sec = str.split(time.strftime(pattern, time.gmtime(SD_df['time'][t])))[1][6:8]
        SD_decimal_time.append(float(hour) + float(min) / 60.0 + float(sec) / 3600.0)
    SD_df['decimal_time'] = SD_decimal_time

    RISR_decimal_time = []
    for t in range(len(RISR_df['time'])):
        hour = str.split(time.strftime(pattern, time.gmtime(RISR_df['time'][t])))[1][0:2]
        min = str.split(time.strftime(pattern, time.gmtime(RISR_df['time'][t])))[1][3:5]
        sec = str.split(time.strftime(pattern, time.gmtime(RISR_df['time'][t])))[1][6:8]
        RISR_decimal_time.append(float(hour) + float(min) / 60.0 + float(sec) / 3600.0)
    RISR_df['decimal_time'] = RISR_decimal_time

    # Set up the plot
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300, nrows=2, ncols=1)
    plt.subplots_adjust(hspace=0.4)
    fig.suptitle("RKN and RISR LOS Velocity Evolution; as produced by " + str(os.path.basename(__file__)),
                 fontsize=12.5)

    # Plot SuperDARN data on the first set of axis
    if SD_beam_range[0] == SD_beam_range[1]:
        SD_beam_string = "Beam " + str(SD_beam_range[0])
    else:
        SD_beam_string = "Beams " + str(SD_beam_range[0]) + "-" + str(SD_beam_range[1])
    if SD_gate_range[0] == SD_gate_range[1]:
        SD_gate_string = "Gate " + str(SD_gate_range[0])
    else:
        SD_gate_string = "Gates " + str(SD_gate_range[0]) + "-" + str(SD_gate_range[1])
    ax[0].title.set_text(SD_station + "; " + SD_beam_string + "; " + SD_gate_string)
    ax[0].set_xlabel('Time [UT]')
    ax[0].set_ylabel('LOS Velocity [m/s]')
    ax[0].set_ylim([-600, 600])
    ax[0].set_xlim([float(start_time[0:2]), float(end_time[0:2])])
    ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0].grid(which='major', axis='y', linestyle='--', linewidth=0.5)
    ax[0].plot(ax[0].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
    ax[0].plot(SD_df['decimal_time'], SD_df['vel'], 'r.', markersize=5)
    ax[0].errorbar(SD_df['decimal_time'], SD_df['vel'], yerr=SD_df['vel_err'],
                   fmt='none', color='black', linewidth=0.75)

    # Plot RISR data on the second plot
    if RISR_wd_beam_range[0] == RISR_wd_beam_range[1]:
        RISR_beam_string = "WD Beam " + str(RISR_wd_beam_range[0])
    else:
        RISR_beam_string = "WD Beams " + str(RISR_wd_beam_range[0]) + "-" + str(RISR_wd_beam_range[1])
    if RISR_range[0] == RISR_range[1]:
        RISR_range_string = "Range " + str(RISR_range[0]) + " km"
    else:
        RISR_range_string = "Range " + str(RISR_range[0]) + "-" + str(RISR_range[1]) + " km"
    ax[1].title.set_text(RISR_station + "; " + RISR_beam_string + "; " + RISR_range_string)
    ax[1].set_xlabel('Time [UT]')
    ax[1].set_ylabel('LOS Ion Velocity [m/s]')
    ax[1].set_ylim([-600, 600])
    ax[1].set_xlim([float(start_time[0:2]), float(end_time[0:2])])
    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[1].grid(which='major', axis='y', linestyle='--', linewidth=0.5)
    ax[1].plot(ax[1].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
    ax[1].plot(RISR_df['decimal_time'], RISR_df['los_ion_vel'], 'r.', markersize=5)
    ax[1].errorbar(RISR_df['decimal_time'], RISR_df['los_ion_vel'], yerr=RISR_df['los_ion_vel_err'],
                   fmt='none', color='black', linewidth=0.75)

    plt.show()

    if SAVE_PLOT:
        cur_path = os.path.dirname(__file__)  # where we are
        fig.savefig(
            cur_path + '/out/vel_vs_time/' + SD_station + "_" + RISR_station + "_vel_vs_time_"
            + year + month + day + " "
            + start_time[0:5].replace(":", ".") + "-" + end_time[0:5].replace(":", ".")
            + "UT.jpg",
            format='jpg', dpi=300)
