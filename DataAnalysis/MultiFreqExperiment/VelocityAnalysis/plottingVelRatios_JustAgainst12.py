import glob
import os
import pathlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from PyPDF2 import PdfFileMerger
import pandas as pd


if __name__ == '__main__':
    """
    This is an older verison of the program that just has 10/12, 13/12, and 14/12
    
    Plot velocity at one frequency against velocity at another frequency
    Plot velocity ratios against time (ratio evolution)
    """

    SAVE_PLOTS = False
    SHOW_PLOTS = True

    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd

    station = "rkn"
    count_min = 4

    numonic = station.upper()

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
    in_dir = loc_root + "/MultiFreqExperiment/VelocityAnalysis/data/" + station
    in_file = in_dir + "/" + station + year + month + day + ".MatchedData.1gg60s.pkl"
    df = pd.read_pickle(in_file)

    out_dir = loc_root + "/MultiFreqExperiment/VelocityAnalysis/out/"

    # We need 6 pages
    if SHOW_PLOTS:
        pages = 1
    else:
        pages = 8
    for time_chunk in range(pages):
        start_time = 0 + time_chunk * 0.5
        end_time = 0.5 + time_chunk * 0.5

        print("Time: " + str(start_time) + " till " + str(end_time))
        time_restricted_df = df[(df['decimalTimes'] >= start_time)
                                & (df['decimalTimes'] < end_time)]

        # Set up the plot
        n_rows = 3
        n_cols = 2
        fig, ax = plt.subplots(figsize=(8, 9), dpi=300, nrows=n_rows, ncols=n_cols)
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
        fig.suptitle("Velocity Frequency Dependence: Velocity Ratios"
                     + "\n" + numonic + " " + year + "." + month + "." + day
                     + " " + str(start_time) + "-" + str(end_time) + " UT"
                     + "\nProduced by " + str(os.path.basename(__file__)),
                     fontsize=13)

        # Format subplots in the first column
        for row in range(n_rows):
            ax[row][0].set_xlabel("12 MHz Velocities [m/s]")
            ax[row][0].set_ylim([-1000, 1000])
            ax[row][0].set_xlim([-1000, 1000])
            ax[row][0].yaxis.set_minor_locator(MultipleLocator(100))
            ax[row][0].xaxis.set_minor_locator(MultipleLocator(100))
            ax[row][0].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
            ax[row][0].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
            ax[row][0].plot(ax[row][0].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
            ax[row][0].plot([0, 0], ax[row][0].get_xlim(), linestyle='-', linewidth=0.5, color='black')
            ax[row][0].plot([ax[row][0].get_ylim()[0], ax[row][0].get_ylim()[1]],
                            [ax[row][0].get_xlim()[0], ax[row][0].get_xlim()[1]],
                            linestyle='-', linewidth=1, color='red')

        # Format subplots in the second column
        for row in range(n_rows):
            ax[row][1].set_xlabel("Time, UT")
            ax[row][1].set_ylim([0, 2])
            ax[row][1].set_xlim([start_time, end_time])
            ax[row][1].yaxis.set_minor_locator(MultipleLocator(1))
            ax[row][1].grid(b=True, which='major', axis='y', linestyle='--', linewidth=0.5)
            ax[row][1].grid(b=True, which='minor', axis='y', linestyle='--', linewidth=0.2)
            ax[row][1].plot(ax[row][1].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
            ax[row][1].plot(ax[row][1].get_ylim(), [1, 1], linestyle='-', linewidth=1, color='red')

        ax[0][0].set_ylabel("10 MHz Velocities [m/s]")
        ax[1][0].set_ylabel("13 MHz Velocities [m/s]")
        ax[2][0].set_ylabel("14 MHz Velocities [m/s]")

        ax[0][1].set_ylabel("Velocity Ratio: 10/12")
        ax[1][1].set_ylabel("Velocity Ratio: 13/12")
        ax[2][1].set_ylabel("Velocity Ratio: 14/12")

        # Plot 10 to 12 Comparison data on the first set of plots
        df_10_12 = time_restricted_df[(time_restricted_df['count10'] >= count_min)
                                      & (time_restricted_df['count12'] >= count_min)
                                      & time_restricted_df['10over12'].notna()]

        ax[0][0].scatter(df_10_12['vel12'], df_10_12['vel10'],
                         s=4, color='k', marker='.', label='gg 10-74')
        ax[0][0].text(-900, 710, 'n=' + str(df_10_12.shape[0]), fontsize=12)
        ax[0][0].legend(loc='lower right')

        ax[0][1].scatter(df_10_12['decimalTimes'], df_10_12['10over12'],
                         s=4, color='k', marker='.')

        # Plot 13 to 12 Comparison data on the second set of plots
        df_13_12 = time_restricted_df[(time_restricted_df['count13'] >= count_min)
                                      & (time_restricted_df['count12'] >= count_min)
                                      & time_restricted_df['13over12'].notna()]

        ax[1][0].scatter(df_13_12['vel12'], df_13_12['vel13'],
                        s=4, color='k', marker='.', label='gg 10-74')
        ax[1][0].text(-900, 710, 'n=' + str(df_13_12.shape[0]), fontsize=12)
        ax[1][0].legend(loc='lower right')

        ax[1][1].scatter(df_13_12['decimalTimes'], df_13_12['13over12'],
                         s=4, color='k', marker='.')

        # Plot 14 to 12 Comparison data on the third set of plots
        df_14_12 = time_restricted_df[(time_restricted_df['count14'] >= count_min)
                                      & (time_restricted_df['count12'] >= count_min)
                                      & time_restricted_df['14over12'].notna()]

        ax[2][0].scatter(df_14_12['vel12'], df_14_12['vel14'],
                         s=4, color='k', marker='.', label='gg 10-74')
        ax[2][0].text(-900, 710, 'n=' + str(df_14_12.shape[0]), fontsize=12)
        ax[2][0].legend(loc='lower right')

        ax[2][1].scatter(df_14_12['decimalTimes'], df_14_12['14over12'],
                         s=4, color='k', marker='.')

        if SHOW_PLOTS:
            plt.show()

        if SAVE_PLOTS:
            fig.savefig(out_dir + "/" + numonic + "_vel_comp_" + year + month + day
                        + "_" + chr(ord('a') + time_chunk) + " " + str(start_time) + "-" + str(end_time) + "UT"
                        + "_temp.pdf", format='pdf', dpi=300)

    if SAVE_PLOTS:
        # Merge all the temp pdf files
        merger = PdfFileMerger()
        for pdf in glob.iglob("out/*_temp.pdf"):
            merger.append(pdf)
        with open(out_dir + "/" + numonic + "_vel_comp_" + year + month + day + ".pdf",
                  "wb") as fout:
            merger.write(fout)
        merger.close()

        for file in glob.iglob("out/*_temp.pdf"):
            os.remove(file)