import glob
import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.ticker import MultipleLocator
from PyPDF2 import PdfFileMerger

if __name__ == '__main__':
    """
    Produces two types of plots:
        Velocity vs Velocity scatter
        Velocity ratio vs time
    
    Data is plotted for a specific gate region
    """

    SAVE_PLOTS = True
    SHOW_PLOTS = False

    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd

    station = "rkn"
    count_min = 4
    gates = [10, 30]
    data_match_type = "Median"
    second_resolution = 60

    start_hour = 0
    end_hour = 4

    mnemonic = station.upper()
    gate_label = "gg: " + str(gates[0]) + "-" + str(gates[1])

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
    in_dir = loc_root + "/MultiFreqExperiment/VelocityAnalysis/data/" + station + "/" + station + year + month + day
    in_file = in_dir + "/" + station + year + month + day + "." + \
              data_match_type + "MatchedData.1gg" + str(second_resolution) + "s.pkl"
    df = pd.read_pickle(in_file)

    # Filter the data based on the expected gate range of the region of interest
    df = df.loc[(df['gate'] >= gates[0]) & (df['gate'] <= gates[1])]
    df.reset_index(drop=True, inplace=True)

    out_dir = loc_root + "/MultiFreqExperiment/VelocityAnalysis/out/"

    normal_ratio_c = "k"
    flagged_ratio_c = "b"
    pts_size = 3

    # For each time period, we need 2 pages
    # For 4 hours, we need 8 time period
    if SHOW_PLOTS:
        time_chunks = 1
    else:
        time_chunks = 8

    for time_chunk in range(time_chunks):
        start_time = start_hour + time_chunk * 0.5
        end_time = (start_hour + 0.5) + time_chunk * 0.5

        print("Plotting from: " + str(start_time) + " till " + str(end_time) + " UT")
        time_restricted_df = df[(df['decimalTime'] >= start_time)
                                & (df['decimalTime'] < end_time)]

        # Build frequency dependant data frames
        df_10_12 = time_restricted_df[(time_restricted_df['count10'] >= count_min)
                                      & (time_restricted_df['count12'] >= count_min)
                                      & time_restricted_df['10over12'].notna()]
        df_13_12 = time_restricted_df[(time_restricted_df['count13'] >= count_min)
                                      & (time_restricted_df['count12'] >= count_min)
                                      & time_restricted_df['13over12'].notna()]
        df_14_12 = time_restricted_df[(time_restricted_df['count14'] >= count_min)
                                      & (time_restricted_df['count12'] >= count_min)
                                      & time_restricted_df['14over12'].notna()]
        df_14_13 = time_restricted_df[(time_restricted_df['count14'] >= count_min)
                                      & (time_restricted_df['count13'] >= count_min)
                                      & time_restricted_df['14over13'].notna()]
        df_13_10 = time_restricted_df[(time_restricted_df['count13'] >= count_min)
                                      & (time_restricted_df['count10'] >= count_min)
                                      & time_restricted_df['13over10'].notna()]
        df_14_10 = time_restricted_df[(time_restricted_df['count14'] >= count_min)
                                      & (time_restricted_df['count10'] >= count_min)
                                      & time_restricted_df['14over10'].notna()]

        # We will have 2 plots, one for vel v vel ratios and one for ratios v time
        # First lets just deal with the first plot

        # Set up the first plot
        n_rows = 3
        n_cols = 2
        fig1, ax1 = plt.subplots(figsize=(8, 9), dpi=300, nrows=n_rows, ncols=n_cols)
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
        fig1.suptitle("Velocity Frequency Dependence: Vel vs Vel Scatter"
                      + "\n" + mnemonic + " " + year + "." + month + "." + day
                      + "; " + gate_label
                      + "; " + str(start_time) + "-" + str(end_time) + " UT"
                      + "\nProduced by " + str(os.path.basename(__file__)),
                      fontsize=13)

        # Format subplots
        for row in range(n_rows):
            ax1[row][0].set_xlabel("12 MHz Velocities [m/s]")
            ax1[row][1].set_xlabel("10 MHz Velocities [m/s]")

            for col in range(n_cols):
                ax1[row][col].set_ylim([-600, 600])
                ax1[row][col].set_xlim([-600, 600])
                ax1[row][col].yaxis.set_minor_locator(MultipleLocator(100))
                ax1[row][col].xaxis.set_minor_locator(MultipleLocator(100))
                ax1[row][col].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
                ax1[row][col].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
                ax1[row][col].plot(ax1[row][col].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
                ax1[row][col].plot([0, 0], ax1[row][col].get_xlim(), linestyle='-', linewidth=0.5, color='black')
                ax1[row][col].plot([ax1[row][col].get_ylim()[0], ax1[row][col].get_ylim()[1]],
                                   [ax1[row][col].get_xlim()[0], ax1[row][col].get_xlim()[1]],
                                   linestyle='-', linewidth=1, color='red')

        ax1[0][1].set_xlabel("13 MHz Velocities [m/s]")

        ax1[0][0].set_ylabel("10 MHz Velocities [m/s]")
        ax1[1][0].set_ylabel("13 MHz Velocities [m/s]")
        ax1[2][0].set_ylabel("14 MHz Velocities [m/s]")

        ax1[0][1].set_ylabel("14 MHz Velocities [m/s]")
        ax1[1][1].set_ylabel("13 MHz Velocities [m/s]")
        ax1[2][1].set_ylabel("14 MHz Velocities [m/s]")

        # Plot 10 to 12 Comparison data in ROW: 0, COL: 0
        ax1[0][0].scatter(df_10_12.loc[df_10_12['diffHeightFlag_10over12'], 'vel12'],
                          df_10_12.loc[df_10_12['diffHeightFlag_10over12'], 'vel10'],
                          s=pts_size, color=flagged_ratio_c, marker='.', zorder=4, label="Frequencies Measuring Different Heights")
        ax1[0][0].scatter(df_10_12.loc[np.logical_not(df_10_12['diffHeightFlag_10over12']), 'vel12'],
                          df_10_12.loc[np.logical_not(df_10_12['diffHeightFlag_10over12']), 'vel10'],
                          s=pts_size, color=normal_ratio_c, marker='.', zorder=3, label="Frequencies Measuring Similar Heights")
        ax1[0][0].text(-490, 425, 'n=' + str(df_10_12.shape[0]), fontsize=12)

        # Plot 13 to 12 Comparison data in ROW: 1, COL: 0
        ax1[1][0].scatter(df_13_12.loc[df_13_12['diffHeightFlag_13over12'], 'vel12'],
                          df_13_12.loc[df_13_12['diffHeightFlag_13over12'], 'vel13'],
                          s=pts_size, color=flagged_ratio_c, marker='.', zorder=4)
        ax1[1][0].scatter(df_13_12.loc[np.logical_not(df_13_12['diffHeightFlag_13over12']), 'vel12'],
                          df_13_12.loc[np.logical_not(df_13_12['diffHeightFlag_13over12']), 'vel13'],
                          s=pts_size, color=normal_ratio_c, marker='.', zorder=3)
        ax1[1][0].text(-490, 425, 'n=' + str(df_13_12.shape[0]), fontsize=12)

        # Plot 14 to 12 Comparison data in ROW: 2, COL: 0
        ax1[2][0].scatter(df_14_12.loc[df_14_12['diffHeightFlag_14over12'], 'vel12'],
                          df_14_12.loc[df_14_12['diffHeightFlag_14over12'], 'vel14'],
                          s=pts_size, color=flagged_ratio_c, marker='.', zorder=4)
        ax1[2][0].scatter(df_14_12.loc[np.logical_not(df_14_12['diffHeightFlag_14over12']), 'vel12'],
                          df_14_12.loc[np.logical_not(df_14_12['diffHeightFlag_14over12']), 'vel14'],
                          s=pts_size, color=normal_ratio_c, marker='.', zorder=3)
        ax1[2][0].text(-490, 425, 'n=' + str(df_14_12.shape[0]), fontsize=12)

        # Plot 14 to 13 Comparison data in ROW: 0, COL: 1
        ax1[0][1].scatter(df_14_13.loc[df_14_13['diffHeightFlag_14over13'], 'vel13'],
                          df_14_13.loc[df_14_13['diffHeightFlag_14over13'], 'vel14'],
                          s=pts_size, color=flagged_ratio_c, marker='.', zorder=4)
        ax1[0][1].scatter(df_14_13.loc[np.logical_not(df_14_13['diffHeightFlag_14over13']), 'vel13'],
                          df_14_13.loc[np.logical_not(df_14_13['diffHeightFlag_14over13']), 'vel14'],
                          s=pts_size, color=normal_ratio_c, marker='.', zorder=3)
        ax1[0][1].text(-490, 425, 'n=' + str(df_14_13.shape[0]), fontsize=12)

        # Plot 13 to 10 Comparison data in ROW: 1, COL: 1
        ax1[1][1].scatter(df_13_10.loc[df_13_10['diffHeightFlag_13over10'], 'vel10'],
                          df_13_10.loc[df_13_10['diffHeightFlag_13over10'], 'vel13'],
                          s=pts_size, color=flagged_ratio_c, marker='.', zorder=3)
        ax1[1][1].scatter(df_13_10.loc[np.logical_not(df_13_10['diffHeightFlag_13over10']), 'vel10'],
                          df_13_10.loc[np.logical_not(df_13_10['diffHeightFlag_13over10']), 'vel13'],
                          s=pts_size, color=normal_ratio_c, marker='.', zorder=3)
        ax1[1][1].text(-490, 425, 'n=' + str(df_13_10.shape[0]), fontsize=12)

        # Plot 14 to 10 Comparison data in ROW: 2, COL: 1
        ax1[2][1].scatter(df_14_10.loc[df_14_10['diffHeightFlag_14over10'], 'vel10'],
                          df_14_10.loc[df_14_10['diffHeightFlag_14over10'], 'vel14'],
                          s=pts_size, color=flagged_ratio_c, marker='.', zorder=4)
        ax1[2][1].scatter(df_14_10.loc[np.logical_not(df_14_10['diffHeightFlag_14over10']), 'vel10'],
                          df_14_10.loc[np.logical_not(df_14_10['diffHeightFlag_14over10']), 'vel14'],
                          s=pts_size, color=normal_ratio_c, marker='.', zorder=3)
        ax1[2][1].text(-490, 425, 'n=' + str(df_14_10.shape[0]), fontsize=12)

        # Set up the second plot
        fig2, ax2 = plt.subplots(figsize=(8, 9), dpi=300, nrows=n_rows, ncols=n_cols)
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
        fig2.suptitle("Velocity Frequency Dependence: Ratios vs Time"
                      + "\n" + mnemonic + " " + year + "." + month + "." + day
                      + " " + str(start_time) + "-" + str(end_time) + " UT"
                      + "\nProduced by " + str(os.path.basename(__file__)),
                      fontsize=13)

        # Format subplots in the second column
        for row in range(n_rows):
            for col in range(n_cols):
                ax2[row][col].set_xlabel("Time, UT")
                ax2[row][col].set_ylim([0, 2])
                ax2[row][col].set_xlim([start_time, end_time])
                ax2[row][col].yaxis.set_minor_locator(MultipleLocator(1))
                ax2[row][col].grid(b=True, which='major', axis='y', linestyle='--', linewidth=0.5)
                ax2[row][col].grid(b=True, which='minor', axis='y', linestyle='--', linewidth=0.2)
                ax2[row][col].plot(ax2[row][col].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
                ax2[row][col].plot(ax2[row][col].get_ylim(), [1, 1], linestyle='-', linewidth=1, color='red')

        ax2[0][0].set_ylabel("Velocity Ratio: 10/12")
        ax2[1][0].set_ylabel("Velocity Ratio: 13/12")
        ax2[2][0].set_ylabel("Velocity Ratio: 14/12")

        ax2[0][1].set_ylabel("Velocity Ratio: 14/13")
        ax2[1][1].set_ylabel("Velocity Ratio: 13/10")
        ax2[2][1].set_ylabel("Velocity Ratio: 14/10")

        # Plot 10 to 12 Comparison data in ROW: 0, COL: 0
        ax2[0][0].scatter(df_10_12['decimalTime'], df_10_12['10over12'], s=4, color='b', marker='.')

        # Plot 13 to 12 Comparison data in ROW: 1, COL: 0
        ax2[1][0].scatter(df_13_12['decimalTime'], df_13_12['13over12'], s=4, color='b', marker='.')

        # Plot 14 to 12 Comparison data in ROW: 2, COL: 0
        ax2[2][0].scatter(df_14_12['decimalTime'], df_14_12['14over12'], s=4, color='b', marker='.')

        # Plot 12 to 10 Comparison data in ROW: 0, COL: 1
        ax2[0][1].scatter(df_14_13['decimalTime'], df_14_13['14over13'], s=4, color='b', marker='.')

        # Plot 13 to 10 Comparison data in ROW: 1, COL: 1
        ax2[1][1].scatter(df_13_10['decimalTime'], df_13_10['13over10'], s=4, color='b', marker='.')

        # Plot 14 to 10 Comparison data in ROW: 2, COL: 1
        ax2[2][1].scatter(df_14_10['decimalTime'], df_14_10['14over10'], s=4, color='b', marker='.')

        # Add legends to the plots
        ax1[0][0].legend(loc=(0, 1.05), ncol=2)

        if SHOW_PLOTS:
            plt.show()

        if SAVE_PLOTS:
            # Save the files as a temp file
            fig1.savefig(out_dir + "/" + mnemonic + "_vel_comp_" + year + month + day
                         + "_" + chr(ord('a') + time_chunk) + " " + str(start_time) + "-" + str(end_time) + "UT"
                         + "_temp1.pdf", format='pdf', dpi=300)
            fig2.savefig(out_dir + "/" + mnemonic + "_vel_comp_" + year + month + day
                         + "_" + chr(ord('a') + time_chunk) + " " + str(start_time) + "-" + str(end_time) + "UT"
                         + "_temp2.pdf", format='pdf', dpi=300)

    if SAVE_PLOTS:
        # Merge all the temp pdf files
        merger = PdfFileMerger()
        for pdf in glob.iglob("out/*_temp1.pdf"):
            merger.append(pdf)
        for pdf in glob.iglob("out/*_temp2.pdf"):
            merger.append(pdf)
        with open(out_dir + "/" + mnemonic + "_vel_comp_" + year + month + day
                  + "_gg" + str(gates[0]) + "-" + str(gates[1]) + ".pdf",
                  "wb") as fout:
            merger.write(fout)
        merger.close()

        for file in glob.iglob("out/*_temp*.pdf"):
            os.remove(file)
