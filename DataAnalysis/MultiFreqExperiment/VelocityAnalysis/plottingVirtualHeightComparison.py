import glob
import os
import pathlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from PyPDF2 import PdfFileMerger
import pandas as pd

if __name__ == '__main__':
    """
    Create scatter plots comparing virtual height at one frequency to virtual height at another frequency
    """

    # TODO: Need to plot the points flagged with the height ratio in a different colour This needs to be done on both
    #  this plot and the vel ratios plots This is needed to flag points where data at one frequency is coming from a
    #  noticeably different height than data at another frequency

    SAVE_PLOTS = True
    SHOW_PLOTS = False

    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd

    station = "rkn"
    count_min = 4
    gates = [19, 21]
    data_match_type = "Median"
    second_resolution = 60

    start_hour = 0
    end_hour = 4

    gate_label = "gg: " + str(gates[0]) + "-" + str(gates[1])
    mnemonic = station.upper()

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
    in_dir = loc_root + "/MultiFreqExperiment/VelocityAnalysis/data/" + station + "/" + station + year + month + day
    in_file = in_dir + "/" + station + year + month + day + "." + \
              data_match_type + "MatchedData.1gg" + str(second_resolution) + "s.pkl"
    df = pd.read_pickle(in_file)

    # Filter the data based on the expected gate range of the region of interest
    df = df.loc[(df['gate'] >= gates[0]) & (df['gate'] <= gates[1])]
    df.reset_index(drop=True, inplace=True)
    print(df)

    out_dir = loc_root + "/MultiFreqExperiment/VelocityAnalysis/out/"

    normal_ratio_c = "k"
    flagged_ratio_c = "b"

    # For 4 hours, we need 8 time periods (half an hour periods)
    if SHOW_PLOTS:
        time_chunks = 1
    else:
        time_chunks = 8

    for time_chunk in range(time_chunks):
        start_time = start_hour + time_chunk * 0.5
        end_time = (start_hour + 0.5) + time_chunk * 0.5

        print("Plotting from: " + str(start_time) + " till " + str(end_time) + " UT")
        time_restricted_df = df[(df['decimalTime'] >= start_time) & (df['decimalTime'] < end_time)]

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

        # Set up the plot
        n_rows = 3
        n_cols = 2
        fig1, ax1 = plt.subplots(figsize=(8, 9), dpi=300, nrows=n_rows, ncols=n_cols)
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
        fig1.suptitle("Echo Frequency Dependence: Virtual Height Comparison"
                      + "\n" + mnemonic + " " + year + "." + month + "." + day
                      + "; " + gate_label
                      + "; " + str(start_time) + "-" + str(end_time) + " UT"
                      + "\nProduced by " + str(os.path.basename(__file__)),
                      fontsize=13)

        # Format subplots
        for row in range(n_rows):
            ax1[row][0].set_xlabel("12 MHz Virtual Heights [km]")
            ax1[row][1].set_xlabel("10 MHz Virtual Heights [km]")

            for col in range(n_cols):
                ax1[row][col].set_ylim([40, 200])
                ax1[row][col].set_xlim([40, 200])
                ax1[row][col].yaxis.set_minor_locator(MultipleLocator(10))
                ax1[row][col].xaxis.set_minor_locator(MultipleLocator(10))
                ax1[row][col].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
                ax1[row][col].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
                ax1[row][col].plot(ax1[row][col].get_ylim(), [120, 120], linestyle='-', linewidth=0.5, color='black')
                ax1[row][col].plot([120, 120], ax1[row][col].get_xlim(), linestyle='-', linewidth=0.5, color='black')
                ax1[row][col].plot([ax1[row][col].get_ylim()[0], ax1[row][col].get_ylim()[1]],
                                   [ax1[row][col].get_xlim()[0], ax1[row][col].get_xlim()[1]],
                                   linestyle='-', linewidth=1, color='red')

        ax1[0][1].set_xlabel("13 MHz Virtual Heights [km]")

        ax1[0][0].set_ylabel("10 MHz Virtual Heights [km]")
        ax1[1][0].set_ylabel("13 MHz Virtual Heights [km]")
        ax1[2][0].set_ylabel("14 MHz Virtual Heights [km]")

        ax1[0][1].set_ylabel("14 MHz Virtual Heights [km]")
        ax1[1][1].set_ylabel("13 MHz Virtual Heights [km]")
        ax1[2][1].set_ylabel("14 MHz Virtual Heights [km]")

        # Plot 10 to 12 Comparison data in ROW: 0, COL: 0
        ax1[0][0].scatter(df_10_12.loc[df_10_12['diffHeightFlag_10over12'], 'height12'],
                          df_10_12.loc[df_10_12['diffHeightFlag_10over12'], 'height10'],
                          s=4, color=flagged_ratio_c, marker='.', zorder=3, label="Ratios Diverging from Unity")
        ax1[0][0].scatter(df_10_12.loc[np.logical_not(df_10_12['diffHeightFlag_10over12']), 'height12'],
                          df_10_12.loc[np.logical_not(df_10_12['diffHeightFlag_10over12']), 'height10'],
                          s=4, color=normal_ratio_c, marker='.', zorder=3, label="Ratios near Unity")
        ax1[0][0].text(52, 172, 'n=' + str(df_10_12.shape[0]), fontsize=12)

        # Plot 13 to 12 Comparison data in ROW: 1, COL: 0
        ax1[1][0].scatter(df_13_12.loc[df_13_12['diffHeightFlag_13over12'], 'height12'],
                          df_13_12.loc[df_13_12['diffHeightFlag_13over12'], 'height13'],
                          s=4, color=flagged_ratio_c, marker='.', zorder=3)
        ax1[1][0].scatter(df_13_12.loc[np.logical_not(df_13_12['diffHeightFlag_13over12']), 'height12'],
                          df_13_12.loc[np.logical_not(df_13_12['diffHeightFlag_13over12']), 'height13'],
                          s=4, color=normal_ratio_c, marker='.', zorder=3)
        ax1[1][0].text(52, 172, 'n=' + str(df_13_12.shape[0]), fontsize=12)

        # Plot 14 to 12 Comparison data in ROW: 2, COL: 0
        ax1[2][0].scatter(df_14_12.loc[df_14_12['diffHeightFlag_14over12'], 'height12'],
                          df_14_12.loc[df_14_12['diffHeightFlag_14over12'], 'height14'],
                          s=4, color=flagged_ratio_c, marker='.', zorder=3)
        ax1[2][0].scatter(df_14_12.loc[np.logical_not(df_14_12['diffHeightFlag_14over12']), 'height12'],
                          df_14_12.loc[np.logical_not(df_14_12['diffHeightFlag_14over12']), 'height14'],
                          s=4, color=normal_ratio_c, marker='.', zorder=3)
        ax1[2][0].text(52, 172, 'n=' + str(df_14_12.shape[0]), fontsize=12)

        # Plot 14 to 13 Comparison data in ROW: 0, COL: 1
        ax1[0][1].scatter(df_14_13.loc[df_14_13['diffHeightFlag_14over13'], 'height13'],
                          df_14_13.loc[df_14_13['diffHeightFlag_14over13'], 'height14'],
                          s=4, color=flagged_ratio_c, marker='.', zorder=3)
        ax1[0][1].scatter(df_14_13.loc[np.logical_not(df_14_13['diffHeightFlag_14over13']), 'height13'],
                          df_14_13.loc[np.logical_not(df_14_13['diffHeightFlag_14over13']), 'height14'],
                          s=4, color=normal_ratio_c, marker='.', zorder=3)
        ax1[0][1].text(52, 172, 'n=' + str(df_14_13.shape[0]), fontsize=12)

        # Plot 13 to 10 Comparison data in ROW: 1, COL: 1
        ax1[1][1].scatter(df_13_10.loc[df_13_10['diffHeightFlag_13over10'], 'height10'],
                          df_13_10.loc[df_13_10['diffHeightFlag_13over10'], 'height13'],
                          s=4, color=flagged_ratio_c, marker='.', zorder=3)
        ax1[1][1].scatter(df_13_10.loc[np.logical_not(df_13_10['diffHeightFlag_13over10']), 'height10'],
                          df_13_10.loc[np.logical_not(df_13_10['diffHeightFlag_13over10']), 'height13'],
                          s=4, color=normal_ratio_c, marker='.', zorder=3)
        ax1[1][1].text(52, 172, 'n=' + str(df_13_10.shape[0]), fontsize=12)

        # Plot 14 to 10 Comparison data in ROW: 2, COL: 1
        ax1[2][1].scatter(df_14_10.loc[df_14_10['diffHeightFlag_14over10'], 'height10'],
                          df_14_10.loc[df_14_10['diffHeightFlag_14over10'], 'height14'],
                          s=4, color=flagged_ratio_c, marker='.', zorder=3)
        ax1[2][1].scatter(df_14_10.loc[np.logical_not(df_14_10['diffHeightFlag_14over10']), 'height10'],
                          df_14_10.loc[np.logical_not(df_14_10['diffHeightFlag_14over10']), 'height14'],
                          s=4, color=normal_ratio_c, marker='.', zorder=3)
        ax1[2][1].text(52, 172, 'n=' + str(df_14_10.shape[0]), fontsize=12)

        # Add legends to the plots
        ax1[0][0].legend(loc=(0, 1.03))

        if SHOW_PLOTS:
            plt.show()

        if SAVE_PLOTS:
            # Save the files as a temp file
            fig1.savefig(out_dir + "/" + mnemonic + "_height_comp_" + year + month + day
                         + "_" + chr(ord('a') + time_chunk) + " " + str(start_time) + "-" + str(end_time) + "UT"
                         + "_temp1.pdf", format='pdf', dpi=300)

    if SAVE_PLOTS:
        # Merge all the temp pdf files
        merger = PdfFileMerger()
        for pdf in glob.iglob("out/*_temp1.pdf"):
            merger.append(pdf)
        with open(out_dir + "/" + mnemonic + "_height_comp_" + year + month + day +
                  "_gg" + str(gates[0]) + "-" + str(gates[1]) + ".pdf",
                  "wb") as fout:
            merger.write(fout)
        merger.close()

        for file in glob.iglob("out/*_temp*.pdf"):
            os.remove(file)
