import calendar
import math
import time
import os
import pathlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy import stats

import pandas as pd

if __name__ == '__main__':
    """
    Plot SuperDARN and RISR velocity vs time
    
    RISR Velocities are flipped and divided by sin(49.9 deg)
    
    """

    SAVE_PLOTS = True
    SHOW_PLOTS = False

    year = "2014"   # yyyy
    month = "03"    # mm
    day = "03"      # dd

    SD_station = "rkn"
    SD_beam_range = [5, 5]
    SD_gate_range = [30, 40]

    RISR_station = "ran"
    # RISR experiments can run for several days, the whole experiment is in one file labeled with the experiment start
    # date.  This experiment start date is not necessarily the date we are plotting.
    RISR_exp_start_month = "3"
    RISR_exp_start_day = "2"
    RISR_wd_beam_range = [2, 2]

    # Compute start and end epochs
    pattern = '%Y.%m.%d %H:%M:%S'  # This is the human readable time pattern we use
    start_date_time = year + "." + month + "." + day + " " + "00:00:00"
    end_date_time = year + "." + month + "." + day + " " + "23:59:59"
    start_epoch = calendar.timegm(time.strptime(start_date_time, pattern))
    end_epoch = calendar.timegm(time.strptime(end_date_time, pattern))

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
    SD_in_dir = loc_root + "/DataReading/SD/data/" + SD_station + "/" + SD_station + year + month + day
    SD_in_file = SD_in_dir + "/" + SD_station + year + month + day + ".pkl"
    SD_df = pd.read_pickle(SD_in_file)

    # Read in RISR data
    RISR_in_dir = loc_root + "/DataReading/RISR/data/" + RISR_station + year + RISR_exp_start_month + RISR_exp_start_day
    RISR_in_file = RISR_in_dir + "/" + RISR_station + year + RISR_exp_start_month + RISR_exp_start_day \
                   + ".LongPulse.pkl"
    RISR_df = pd.read_pickle(RISR_in_file)

    # Filter SuperDARN data
    SD_df = SD_df.loc[(SD_df['time'] >= start_epoch) & (SD_df['time'] <= end_epoch)  # filer times
                      & (SD_df['gate'] >= SD_gate_range[0]) & (SD_df['gate'] <= SD_gate_range[1])  # filter gates
                      & (SD_df['bmnum'] >= SD_beam_range[0]) & (SD_df['bmnum'] <= SD_beam_range[1])]  # filter beams
    SD_df.reset_index(drop=True, inplace=True)

    # Filter RISR data
    RISR_df = RISR_df.loc[(RISR_df['time'] >= start_epoch) & (RISR_df['time'] <= end_epoch)  # filer times
                          & (RISR_df['wd_bmnum'] >= RISR_wd_beam_range[0])
                          & (RISR_df['wd_bmnum'] <= RISR_wd_beam_range[1])]  # filter beams
    RISR_df.reset_index(drop=True, inplace=True)

    # RISR-N and RKN velocities are going opposite directions, flip RISR so toward is +ve
    RISR_df['los_ion_vel'] = RISR_df['los_ion_vel'].apply(lambda x: x * -1)

    # The Beams are not quite aligned
    # TODO: Confirm the angle between beams - might change based on magnetic field direction
    if SD_station == "rkn":
        RISR_df['los_ion_vel'] = RISR_df['los_ion_vel'].apply(lambda x: x / math.sin(math.radians(49.9)))
    else:
        print("Waring: No adjustments have been made to compensate for the angle between beams.")

    # Recover decimal times from epoch times
    SD_decimal_time = []
    RISR_decimal_time = []
    for t in range(len(SD_df['time'])):
        hour = str.split(time.strftime(pattern, time.gmtime(SD_df['time'][t])))[1][0:2]
        min = str.split(time.strftime(pattern, time.gmtime(SD_df['time'][t])))[1][3:5]
        sec = str.split(time.strftime(pattern, time.gmtime(SD_df['time'][t])))[1][6:8]
        SD_decimal_time.append(float(hour) + float(min) / 60.0 + float(sec) / 3600.0)
    for t in range(len(RISR_df['time'])):
        hour = str.split(time.strftime(pattern, time.gmtime(RISR_df['time'][t])))[1][0:2]
        min = str.split(time.strftime(pattern, time.gmtime(RISR_df['time'][t])))[1][3:5]
        sec = str.split(time.strftime(pattern, time.gmtime(RISR_df['time'][t])))[1][6:8]
        RISR_decimal_time.append(float(hour) + float(min) / 60.0 + float(sec) / 3600.0)
    SD_df['decimal_time'] = SD_decimal_time
    RISR_df['decimal_time'] = RISR_decimal_time

    # Build title strings
    if SD_beam_range[0] == SD_beam_range[1]:
        SD_beam_string = "Beam " + str(SD_beam_range[0])
    else:
        SD_beam_string = "Beams " + str(SD_beam_range[0]) + "-" + str(SD_beam_range[1])
    if SD_gate_range[0] == SD_gate_range[1]:
        SD_gate_string = "Gate " + str(SD_gate_range[0])
    else:
        SD_gate_string = "Gates " + str(SD_gate_range[0]) + "-" + str(SD_gate_range[1])

    if RISR_wd_beam_range[0] == RISR_wd_beam_range[1]:
        RISR_beam_string = "WD Beam " + str(RISR_wd_beam_range[0])
    else:
        RISR_beam_string = "WD Beams " + str(RISR_wd_beam_range[0]) + "-" + str(RISR_wd_beam_range[1])

    # Ensure out directory
    out_dir = loc_root + "/Plotting/vel_vs_time/out/" + year + month + day
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Loop thorough and plot 2 hour chunks of data
    length_of_chunks_h = 2
    num_of_chunks = int(24 / length_of_chunks_h)
    for chunk_num in range(num_of_chunks):

        # Computer start and end epochs and build restricted data frames
        start_hour_here = 0 + 2 * chunk_num
        end_hour_here = 2 + 2 * chunk_num
        start_date_time_here = year + "." + month + "." + day + " " + str(start_hour_here) + ":00:00"
        if end_hour_here == 24:
            end_date_time_here = year + "." + month + "." + day + " 23:59:59"
        else:
            end_date_time_here = year + "." + month + "." + day + " " + str(end_hour_here) + ":00:00"
        start_epoch_here = calendar.timegm(time.strptime(start_date_time_here, pattern))
        end_epoch_here = calendar.timegm(time.strptime(end_date_time_here, pattern))

        # Build a restricted data frame based just on the times here
        restricted_SD_df = SD_df.loc[(SD_df['time'] >= start_epoch_here) & (SD_df['time'] <= end_epoch_here)]
        restricted_RISR_df = RISR_df.loc[(RISR_df['time'] >= start_epoch_here) & (RISR_df['time'] <= end_epoch_here)]

        # We have to remove all nan values because median can't handle them
        restricted_RISR_df = restricted_RISR_df[restricted_RISR_df['los_ion_vel'].notna()]

        # Compute binned medians and standard deviations
        n_bins = int(length_of_chunks_h / (1 / 12))  # 5 minute (1/12 hour) bins
        try:
            SD_bin_medians, SD_bin_edges, SD_binnumber = stats.binned_statistic(
                restricted_SD_df['decimal_time'], restricted_SD_df['vel'], 'median', bins=n_bins,
                range=(start_hour_here, end_hour_here))
            SD_bin_stds, x, xx = stats.binned_statistic(
                restricted_SD_df['decimal_time'], restricted_SD_df['vel'], 'std', bins=n_bins,
                range=(start_hour_here, end_hour_here))
            SD_bin_width = (SD_bin_edges[1] - SD_bin_edges[0])
            SD_bin_centers = SD_bin_edges[1:] - SD_bin_width / 2
        except:
            print("Warning for " + str(start_hour_here) + "-" + str(end_hour_here) + " UT: " + SD_station
                  + " has no data here")
            SD_bin_medians, SD_bin_centers, SD_bin_stds = [], [], []

        try:
            RISR_bin_medians, RISR_bin_edges, RISR_binnumber = stats.binned_statistic(
                restricted_RISR_df['decimal_time'], restricted_RISR_df['los_ion_vel'], 'median', bins=n_bins,
                range=(start_hour_here, end_hour_here))
            RISR_bin_stds, x, xx = stats.binned_statistic(
                restricted_RISR_df['decimal_time'], restricted_RISR_df['los_ion_vel'], 'std', bins=n_bins,
                range=(start_hour_here, end_hour_here))
            RISR_bin_width = (RISR_bin_edges[1] - RISR_bin_edges[0])
            RISR_bin_centers = RISR_bin_edges[1:] - RISR_bin_width / 2
        except:
            print("Warning for " + str(start_hour_here) + "-" + str(end_hour_here) + " UT: " + RISR_station
                  + " has no data here")
            RISR_bin_medians, RISR_bin_centers, RISR_bin_stds = [], [], []

        # Set up the plot
        fig, ax = plt.subplots(figsize=(8, 9), dpi=300, nrows=3, ncols=1)
        plt.subplots_adjust(hspace=0.4)
        fig.suptitle(SD_station + " and " + RISR_station + " LOS Velocity Evolution; "
                     + str(start_hour_here) + "-" + str(end_hour_here) + " UT"
                     + "\nNote: Positive Velocity Means Towards the Radars"
                     + "\nProduced by " + str(os.path.basename(__file__)),
                     fontsize=13)

        # Apply common subplot formatting
        for i in range(ax.size):
            ax[i].set_ylim([-600, 600])
            ax[i].set_xlim([start_hour_here, end_hour_here])
            ax[i].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax[i].grid(which='major', axis='y', linestyle='--', linewidth=0.5)
            ax[i].plot(ax[i].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
            ax[i].set_xlabel('Time [UT]')
            ax[i].set_ylabel('LOS Velocity [m/s]')

        # Plot SuperDARN data on the first set of axis
        ax[0].title.set_text(SD_station + "; " + SD_beam_string + "; " + SD_gate_string)
        ax[0].scatter(SD_bin_centers, SD_bin_medians, marker='o', s=40, facecolor='none',
                      edgecolors='r', linewidths=2, label='Binned Medians')
        ax[0].errorbar(SD_bin_centers, SD_bin_medians, yerr=SD_bin_stds,
                       fmt='none', color='red', linewidth=1)
        ax[0].scatter(restricted_SD_df['decimal_time'], restricted_SD_df['vel'], s=4, facecolor='none',
                      edgecolors='k', linewidths=0.75, label='Raw Scatter')
        ax[0].errorbar(restricted_SD_df['decimal_time'], restricted_SD_df['vel'], yerr=restricted_SD_df['vel_err'],
                       fmt='none', color='black', linewidth=0.75)

        # Plot RISR data on the second plot
        ax[1].title.set_text(RISR_station + "; " + RISR_beam_string)
        ax[1].scatter(RISR_bin_centers, RISR_bin_medians, marker='D', s=40, facecolor='none',
                      edgecolors='b', linewidths=2, label='Binned Medians')
        ax[1].errorbar(RISR_bin_centers, RISR_bin_medians, yerr=RISR_bin_stds,
                       fmt='none', color='blue', linewidth=1)
        ax[1].scatter(restricted_RISR_df['decimal_time'], restricted_RISR_df['los_ion_vel'], s=4, facecolor='none',
                      edgecolors='k', linewidths=0.75, label='Raw Scatter')
        ax[1].errorbar(restricted_RISR_df['decimal_time'], restricted_RISR_df['los_ion_vel'],
                       yerr=restricted_RISR_df['los_ion_vel_err'], fmt='none', color='black', linewidth=0.75)

        # Plot both sets of medians on the third subplot
        ax[2].title.set_text(SD_station + " " + RISR_station + " Velocity Comparison")
        try:
            ax[2].scatter(SD_bin_centers + 0.005, SD_bin_medians, marker='o', s=40, facecolor='none',
                          edgecolors='r', linewidths=2, label=SD_station + ' Binned Medians')
            ax[2].errorbar(SD_bin_centers + 0.005, SD_bin_medians, yerr=SD_bin_stds,
                           fmt='none', color='red', linewidth=1)
        except:
            pass
        try:
            ax[2].scatter(RISR_bin_centers - 0.005, RISR_bin_medians, marker='D', s=40, facecolor='none',
                          edgecolors='b', linewidths=2, label=RISR_station + ' Binned Medians')
            ax[2].errorbar(RISR_bin_centers - 0.005, RISR_bin_medians, yerr=RISR_bin_stds,
                           fmt='none', color='blue', linewidth=1)
        except:
            pass

        # Add in legends
        for i in range(ax.size):
            # Shrink the current axis to make space for the legend
            box = ax[i].get_position()
            ax[i].set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])
            try:
                ax[i].legend(bbox_to_anchor=(0.5, -0.50), fancybox=True, ncol=2, loc='lower center')
            except:
                pass

        if SHOW_PLOTS:
            plt.show()

        if SAVE_PLOTS:
            cur_path = os.path.dirname(__file__)  # where we are
            fig.savefig(out_dir + "/" + SD_station + "_" + RISR_station + "_vel_vs_time_"
                        + year + month + day + " " + str(start_hour_here) + "-" + str(end_hour_here) + "UT"
                        + ".pdf", format='pdf', dpi=300)
