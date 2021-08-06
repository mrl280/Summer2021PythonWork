import calendar
import glob
import bz2
import time
import os
import pathlib

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.ticker import MultipleLocator
from scipy import stats
from PyPDF2 import PdfFileMerger

from lib.basic_SD_df_filter import basic_SD_df_filter

if __name__ == '__main__':
    """
    Plot SuperDARN and RISR velocity vs velocity
    Produces plots for a day
    One scatter plot for each 2 hour time period
    """

    SAVE_PLOTS = True
    SHOW_PLOTS = False

    year = "2014"  # yyyy
    month = "03"  # mm
    day = "03"  # dd

    SD_station = "rkn"
    SD_beam_range = [5, 5]
    # SD_gate_range = [30, 40]
    SD_gate_range = [25, 32]
    # SD_gate_range = [31, 41]

    RISR_station = "ran"
    RISR_wd_beam_range = [5, 5]
    resolution = 5  # minutes

    SD_numonic = SD_station.upper()
    if RISR_station == "ran":
        RISR_numonic = "RISR-N"
    elif RISR_station == "ras":
        RISR_numonic = "RISR-C"
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
    SD_in_file = SD_in_dir + "/" + SD_station + year + month + day + ".pbz2"

    data_stream = bz2.BZ2File(SD_in_file, "rb")
    SD_df = pd.read_pickle(data_stream)

    # Read in RISR data
    RISR_in_dir = loc_root + "/DataReading/RISR/data/" + RISR_station + "/" + RISR_station + year + month + day
    RISR_in_file = RISR_in_dir + "/" + RISR_station + year + month + day + "." + str(resolution) + "min.pbz2"

    data_stream = bz2.BZ2File(RISR_in_file, "rb")
    RISR_df = pd.read_pickle(data_stream)

    # Filter SuperDARN data
    SD_df = basic_SD_df_filter(SD_df)
    SD_df = SD_df.loc[(SD_df['epoch'] >= start_epoch) & (SD_df['epoch'] <= end_epoch)  # filer times
                      & (SD_df['gate'] >= SD_gate_range[0]) & (SD_df['gate'] <= SD_gate_range[1])  # filter gates
                      & (SD_df['bmnum'] >= SD_beam_range[0]) & (SD_df['bmnum'] <= SD_beam_range[1])]  # filter beams
    SD_df.reset_index(drop=True, inplace=True)

    # Filter RISR data
    RISR_df = RISR_df.loc[(RISR_df['epoch'] >= start_epoch) & (RISR_df['epoch'] <= end_epoch)  # filer times
                          & (RISR_df['wdBmnum'] >= RISR_wd_beam_range[0])
                          & (RISR_df['wdBmnum'] <= RISR_wd_beam_range[1])]  # filter beams

    # We have to remove all nan values from RISR because several of the upcoming functions can't handle them
    RISR_df = RISR_df.loc[RISR_df['losIonVel'].notna()]
    RISR_df.reset_index(drop=True, inplace=True)

    # RISR-N and RKN velocities are going opposite directions, flip RISR so toward is +ve
    RISR_df['losIonVel'] = RISR_df['losIonVel'].apply(lambda x: x * -1)

    # Original method: divide by sin of the elevation angle
    # RISR_df['losIonVel'] = np.divide(RISR_df['losIonVel'], np.sin((RISR_df['elv']) * np.pi / 180))

    # SuperDARN measures the velocity perpendicular to the magnetic field
    # Therefore we need to find the RISR component in this direction
    RISR_df['delta'] = RISR_df['aspect'] - RISR_df['elv'] - 90
    # TODO: This might need to change based on spherical geometry of the Earth, idk
    RISR_df['losIonVel'] = np.divide(RISR_df['losIonVel'], np.cos((RISR_df['elv'] + RISR_df['delta']) * np.pi / 180))

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
    out_dir = loc_root + "/OneAndOneHalfHop/VelocityAnalysis/out/" + year + month + day
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Loop thorough and plot 2 hour chunks of data
    # That means there will be 6 chunks per page
    if SHOW_PLOTS:
        # If we are looking at the plots we only need the first one
        num_of_chunks_per_page = 1
    else:
        num_of_chunks_per_page = 6
    # There are going to be two pages of plots, each page with 6 subplots
    for page in range(2):

        # Set up the plot
        n_rows = 3
        n_cols = 2
        fig, ax = plt.subplots(figsize=(8, 9), dpi=300, nrows=n_rows, ncols=n_cols)
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
        fig.suptitle(SD_numonic + " and " + RISR_numonic + " LOS Velocity Comparison; " + year + "." + month + "." + day
                     + "\nNote: Positive Velocity Means Towards the Radars; " + SD_gate_string
                     + "\nProduced by " + str(os.path.basename(__file__)),
                     fontsize=13)

        # Apply common subplot formatting
        for row in range(n_rows):
            for col in range(n_cols):
                ax[row][col].set_ylim([-1000, 1000])
                ax[row][col].set_xlim([-1000, 1000])
                ax[row][col].yaxis.set_minor_locator(MultipleLocator(100))
                ax[row][col].xaxis.set_minor_locator(MultipleLocator(100))
                ax[row][col].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
                ax[row][col].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
                ax[row][col].plot(ax[row][col].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
                ax[row][col].plot([0, 0], ax[row][col].get_xlim(), linestyle='-', linewidth=0.5, color='black')
                ax[row][col].plot([ax[row][col].get_ylim()[0], ax[row][col].get_ylim()[1]],
                                  [ax[row][col].get_xlim()[0], ax[row][col].get_xlim()[1]],
                                  linestyle='-', linewidth=1, color='red')
                ax[row][col].set_xlabel(RISR_numonic + ' LOS Velocity [m/s]')
                ax[row][col].set_ylabel(SD_numonic + ' LOS Velocity [m/s]')

        row = -1
        col = 1
        for chunk_num in range(num_of_chunks_per_page):
            if (chunk_num - 1) % 2:
                row = row + 1
            if col == 1:
                col = 0
            else:
                col = 1

            # Computer start and end epochs and build restricted data frames
            start_hour_here = (page * 12) + 0 + 2 * chunk_num
            end_hour_here = (page * 12) + 2 + 2 * chunk_num
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
            n_bins = int(2 / (1 / 12))  # 5 minute (1/12 hour) bins
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
            stats_df = pd.DataFrame({'SDbinCenters': SD_bin_centers,
                                     'SDmedians': SD_bin_medians,
                                     'SDstdDev': SD_bin_stds,
                                     'SDcount': SD_counts,
                                     'RISRbinCenters': RISR_bin_centers,
                                     'RISRmedians': RISR_bin_medians,
                                     'RISRstdDev': RISR_bin_stds,
                                     'RISRcount': RISR_counts
                                     })

            # Remove all points where there are less than three raw data points
            stats_df = stats_df.loc[(stats_df['SDcount'] >= 3) & (stats_df['RISRcount'] >= 3)]

            # Plot SD vs RISR Medians
            ax[row][col].title.set_text(str(start_hour_here) + "-" + str(end_hour_here) + " UT")
            ax[row][col].scatter(stats_df['RISRmedians'], stats_df['SDmedians'], marker='^', s=25, facecolor='none',
                          edgecolors='g', linewidths=1, label='Binned Medians')

        if SHOW_PLOTS:
            plt.show()

        if SAVE_PLOTS:
            fig.savefig(out_dir + "/" + SD_numonic + "_" + RISR_numonic + "_vel_vs_vel_" + year + month + day
                         + "_" + str(page) + "_temp.pdf", format='pdf', dpi=300)

    if SAVE_PLOTS:
        # Merge all the temp pdf files
        merger = PdfFileMerger()
        for pdf in glob.iglob("out/" + year + month + day + "/*_temp.pdf"):
            merger.append(pdf)
        with open(out_dir + "/" + SD_numonic + "_" + RISR_numonic + "_vel_vs_vel_" + year + month + day + ".pdf",
                  "wb") as fout:
            merger.write(fout)
        merger.close()

        for file in glob.iglob("out/" + year + month + day + "/*_temp.pdf"):
            os.remove(file)

