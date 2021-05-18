import calendar
import glob
import numpy as np
import time
import os
import pathlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy import stats
from PyPDF2 import PdfFileMerger
from matplotlib.ticker import MultipleLocator

import pandas as pd

from DataAnalysis.DataReading.SD.basic_SD_df_filter import basic_SD_df_filter

if __name__ == '__main__':
    """
    Plot SuperDARN and RISR_HDF5 velocity vs time
    
    RISR_HDF5 Velocities are flipped and divided by the sin of the elevation angle
    
    """

    SAVE_PLOTS = True
    SHOW_PLOTS = False

    year = "2014"  # yyyy
    month = "03"  # mm
    day = "04"  # dd

    SD_station = "rkn"
    SD_beam_range = [5, 5]
    SD_gate_range = [30, 40]

    RISR_station = "ran"
    # RISR_HDF5 experiments can run for several days, the whole experiment is in one file labeled with the experiment start
    # date.  This experiment start date is not necessarily the date we are plotting.
    RISR_exp_start_month = "3"
    RISR_exp_start_day = "2"
    RISR_wd_beam_range = [2, 2]

    SD_numonic = SD_station.upper()
    if RISR_station == "ran":
        RISR_numonic = "RISR_HDF5-N"
    elif RISR_station == "ras":
        RISR_numonic = "RISR_HDF5-C"
    else:
        raise Exception("Error, " + RISR_station + " not recognized.")

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

    # Read in RISR_HDF5 data
    RISR_in_dir = loc_root + "/DataReading/RISR_HDF5/data/" + RISR_station + "/" + RISR_station + year + RISR_exp_start_month + RISR_exp_start_day
    RISR_in_file = RISR_in_dir + "/" + RISR_station + year + RISR_exp_start_month + RISR_exp_start_day \
                   + ".5min.LongPulse.pkl"
    RISR_df = pd.read_pickle(RISR_in_file)

    # Filter SuperDARN data
    SD_df = basic_SD_df_filter(SD_df)
    SD_df = SD_df.loc[(SD_df['epoch'] >= start_epoch) & (SD_df['epoch'] <= end_epoch)  # filer times
                      & (SD_df['gate'] >= SD_gate_range[0]) & (SD_df['gate'] <= SD_gate_range[1])  # filter gates
                      & (SD_df['bmnum'] >= SD_beam_range[0]) & (SD_df['bmnum'] <= SD_beam_range[1])]  # filter beams
    SD_df.reset_index(drop=True, inplace=True)

    # Filter RISR_HDF5 data
    RISR_df = RISR_df.loc[(RISR_df['epoch'] >= start_epoch) & (RISR_df['epoch'] <= end_epoch)  # filer times
                          & (RISR_df['wdBmnum'] >= RISR_wd_beam_range[0])
                          & (RISR_df['wdBmnum'] <= RISR_wd_beam_range[1])]  # filter beams

    # We have to remove all nan values from RISR_HDF5 because several of the upcoming functions can't handle them
    RISR_df = RISR_df.loc[RISR_df['losIonVel'].notna()]
    RISR_df.reset_index(drop=True, inplace=True)

    # RISR_HDF5-N and RKN velocities are going opposite directions, flip RISR_HDF5 so toward is +ve
    RISR_df['losIonVel'] = RISR_df['losIonVel'].apply(lambda x: x * -1)

    # SuperDARN measures horizontal plasma flow, but RISR_HDF5 only sees a component
    # Therefore we need to divide RISR_HDF5 velocities by the sin of the elevation angle
    # Note: elevation angle is in degrees
    # TODO: This assumes the magnetic field lines are perpendicular, which might not be the case
    #   Koustov to confirm exactly what angle needs to be used
    RISR_df['losIonVel'] = np.divide(RISR_df['losIonVel'], np.sin((RISR_df['elv']) * np.pi / 180))

    # Recover decimal times from epoch times
    SD_decimal_time = []
    RISR_decimal_time = []
    for t in range(len(SD_df['epoch'])):
        hour = str.split(time.strftime(pattern, time.gmtime(SD_df['epoch'][t])))[1][0:2]
        min = str.split(time.strftime(pattern, time.gmtime(SD_df['epoch'][t])))[1][3:5]
        sec = str.split(time.strftime(pattern, time.gmtime(SD_df['epoch'][t])))[1][6:8]
        SD_decimal_time.append(float(hour) + float(min) / 60.0 + float(sec) / 3600.0)
    for t in range(len(RISR_df['epoch'])):
        hour = str.split(time.strftime(pattern, time.gmtime(RISR_df['epoch'][t])))[1][0:2]
        min = str.split(time.strftime(pattern, time.gmtime(RISR_df['epoch'][t])))[1][3:5]
        sec = str.split(time.strftime(pattern, time.gmtime(RISR_df['epoch'][t])))[1][6:8]
        RISR_decimal_time.append(float(hour) + float(min) / 60.0 + float(sec) / 3600.0)
    SD_df['decimalTime'] = SD_decimal_time
    RISR_df['decimalTime'] = RISR_decimal_time

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
    out_dir = loc_root + "/One-and-One-Half-Hop/VelocityAnalysis/out/" + year + month + day
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Loop thorough and plot 2 hour chunks of data
    length_of_chunks_h = 2
    if SHOW_PLOTS:
        # If we are looking at the plots we only need the first one
        num_of_chunks = 1
    else:
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
        restricted_SD_df = SD_df.loc[(SD_df['epoch'] >= start_epoch_here) & (SD_df['epoch'] <= end_epoch_here)]
        restricted_RISR_df = RISR_df.loc[(RISR_df['epoch'] >= start_epoch_here) & (RISR_df['epoch'] <= end_epoch_here)]

        # Compute binned medians and standard deviations
        n_bins = int(length_of_chunks_h / (1 / 12))  # 5 minute (1/12 hour) bins
        try:
            SD_bin_medians, SD_bin_edges, SD_binnumber = stats.binned_statistic(
                restricted_SD_df['decimalTime'], restricted_SD_df['vel'], 'median', bins=n_bins,
                range=(start_hour_here, end_hour_here))
            SD_bin_stds, x, xx = stats.binned_statistic(
                restricted_SD_df['decimalTime'], restricted_SD_df['vel'], 'std', bins=n_bins,
                range=(start_hour_here, end_hour_here))
            SD_counts, xxx, xxxx = stats.binned_statistic(
                restricted_SD_df['decimalTime'], restricted_SD_df['vel'], 'count', bins=n_bins,
                range=(start_hour_here, end_hour_here))
            SD_bin_width = (SD_bin_edges[1] - SD_bin_edges[0])
            SD_bin_centers = SD_bin_edges[1:] - SD_bin_width / 2
        except:
            print("Warning for " + str(start_hour_here) + "-" + str(end_hour_here) + " UT: " + SD_numonic
                  + " has no data here")
            SD_bin_medians, SD_bin_centers, SD_bin_stds, SD_counts = [], [], [], []

        try:
            RISR_bin_medians, RISR_bin_edges, RISR_binnumber = stats.binned_statistic(
                restricted_RISR_df['decimalTime'], restricted_RISR_df['losIonVel'], 'median', bins=n_bins,
                range=(start_hour_here, end_hour_here))
            RISR_bin_stds, x, xx = stats.binned_statistic(
                restricted_RISR_df['decimalTime'], restricted_RISR_df['losIonVel'], 'std', bins=n_bins,
                range=(start_hour_here, end_hour_here))
            RISR_counts, xxx, xxxx = stats.binned_statistic(
                restricted_RISR_df['decimalTime'], restricted_RISR_df['losIonVel'], 'count', bins=n_bins,
                range=(start_hour_here, end_hour_here))
            RISR_bin_width = (RISR_bin_edges[1] - RISR_bin_edges[0])
            RISR_bin_centers = RISR_bin_edges[1:] - RISR_bin_width / 2
        except:
            print("Warning for " + str(start_hour_here) + "-" + str(end_hour_here) + " UT: " + RISR_numonic
                  + " has no data here")
            RISR_bin_medians, RISR_bin_centers, RISR_bin_stds, RISR_counts = [], [], [], []

        # Put the statistics into a df
        SD_stats_df = pd.DataFrame({'binCenters': SD_bin_centers,
                                    'medians': SD_bin_medians,
                                    'stdDev': SD_bin_stds,
                                    'count': SD_counts
                                    })
        RISR_stats_df = pd.DataFrame({'binCenters': RISR_bin_centers,
                                      'medians': RISR_bin_medians,
                                      'stdDev': RISR_bin_stds,
                                      'count': RISR_counts
                                      })
        # Remove all points where there are less than three raw data points
        SD_stats_df = SD_stats_df.loc[(SD_stats_df['count'] >= 3)]
        RISR_stats_df = RISR_stats_df.loc[RISR_stats_df['count'] >= 3]
        SD_stats_df.reset_index(drop=True, inplace=True)
        RISR_stats_df.reset_index(drop=True, inplace=True)

        # Set up the plot
        fig, ax = plt.subplots(figsize=(8, 9), dpi=300, nrows=3, ncols=1)
        plt.subplots_adjust(hspace=0.4)
        fig.suptitle(SD_numonic + " and " + RISR_numonic + " LOS Velocity Evolution; " + year + "." + month + "." + day
                     + "; " + str(start_hour_here) + "-" + str(end_hour_here) + " UT"
                     + "\nNote: Positive Velocity Means Towards the Radars"
                     + "\nProduced by " + str(os.path.basename(__file__)),
                     fontsize=13)

        # Apply common subplot formatting
        for i in range(ax.size):
            ax[i].yaxis.set_minor_locator(MultipleLocator(100))
            ax[i].set_ylim([-1000, 1000])
            ax[i].set_xlim([start_hour_here, end_hour_here])
            ax[i].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax[i].grid(b=True, which='major', axis='y', linestyle='--', linewidth=0.5)
            ax[i].grid(b=True, which='minor', axis='y', linestyle='--', linewidth=0.2)
            ax[i].plot(ax[i].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
            ax[i].set_xlabel('Time [UT]')
            ax[i].set_ylabel('LOS Velocity [m/s]')

        # Plot SuperDARN data on the first set of axis
        ax[0].title.set_text(SD_numonic + "; " + SD_beam_string + "; " + SD_gate_string)
        ax[0].scatter(SD_stats_df['binCenters'], SD_stats_df['medians'], marker='o', s=40, facecolor='none',
                      edgecolors='r', linewidths=2, label='Binned Medians')
        ax[0].errorbar(SD_stats_df['binCenters'], SD_stats_df['medians'], yerr=SD_stats_df['stdDev'],
                       fmt='none', color='red', linewidth=1)
        ax[0].scatter(restricted_SD_df['decimalTime'], restricted_SD_df['vel'], s=4, facecolor='none',
                      edgecolors='k', linewidths=0.75, label='Raw Scatter')
        ax[0].errorbar(restricted_SD_df['decimalTime'], restricted_SD_df['vel'], yerr=restricted_SD_df['velErr'],
                       fmt='none', color='black', linewidth=0.75)

        # Plot RISR_HDF5 data on the second plot
        ax[1].title.set_text(RISR_numonic + "; " + RISR_beam_string)
        ax[1].scatter(RISR_stats_df['binCenters'], RISR_stats_df['medians'], marker='D', s=40, facecolor='none',
                      edgecolors='b', linewidths=2, label='Binned Medians')
        ax[1].errorbar(RISR_stats_df['binCenters'], RISR_stats_df['medians'], yerr=RISR_stats_df['stdDev'],
                       fmt='none', color='blue', linewidth=1)
        ax[1].scatter(restricted_RISR_df['decimalTime'], restricted_RISR_df['losIonVel'], s=4, facecolor='none',
                      edgecolors='k', linewidths=0.75, label='Raw Scatter')
        ax[1].errorbar(restricted_RISR_df['decimalTime'], restricted_RISR_df['losIonVel'],
                       yerr=restricted_RISR_df['losIonVelErr'], fmt='none', color='black', linewidth=0.75)

        # Plot both sets of medians on the third subplot
        ax[2].title.set_text(SD_numonic + " " + RISR_numonic + " Velocity Comparison")
        try:
            ax[2].scatter(SD_stats_df['binCenters'] + 0.005, SD_stats_df['medians'], marker='o', s=40, facecolor='none',
                          edgecolors='r', linewidths=2, label=SD_numonic + ' Binned Medians')
            ax[2].errorbar(SD_stats_df['binCenters'] + 0.005, SD_stats_df['medians'], yerr=SD_stats_df['stdDev'],
                           fmt='none', color='red', linewidth=1)
        except:
            pass
        try:
            ax[2].scatter(RISR_stats_df['binCenters'] - 0.005, RISR_stats_df['medians'], marker='D', s=40, facecolor='none',
                          edgecolors='b', linewidths=2, label=RISR_numonic + ' Binned Medians')
            ax[2].errorbar(RISR_stats_df['binCenters'] - 0.005, RISR_stats_df['medians'], yerr=RISR_stats_df['stdDev'],
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
            fig.savefig(out_dir + "/" + SD_numonic + "_" + RISR_numonic + "_vel_vs_time_" + year + month + day
                        + "_" + chr(ord('a') + chunk_num) + " " + str(start_hour_here) + "-" + str(end_hour_here) + "UT"
                        + "_temp.pdf", format='pdf', dpi=300)

    if SAVE_PLOTS:
        # Merge all the temp pdf files
        merger = PdfFileMerger()
        for pdf in glob.iglob("out/" + year + month + day + "/*_temp.pdf"):
            merger.append(pdf)
        with open(out_dir + "/" + SD_numonic + "_" + RISR_numonic + "_vel_vs_time_" + year + month + day + ".pdf",
                  "wb") as fout:
            merger.write(fout)
        merger.close()

        for file in glob.iglob("out/" + year + month + day + "/*_temp.pdf"):
            os.remove(file)
