import glob
import os
import pathlib

import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.ticker import MultipleLocator
from PyPDF2 import PdfFileMerger

if __name__ == '__main__':
    """
    Create scatter plots comparing elevation angles at one frequency to elevation angles at another frequency
    """

    SAVE_PLOTS = True
    SHOW_PLOTS = False

    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd

    station = "rkn"
    gates = [10, 30]
    data_match_type = "Median"     # "Median" or "Raw"
    count_min = 4   # Only used for median matched data

    start_hour = 0
    end_hour = 4

    gate_label = "gg: " + str(gates[0]) + "-" + str(gates[1])
    mnemonic = station.upper()

    if data_match_type == "Median":
        second_resolution = 60
    elif data_match_type == "Raw":
        second_resolution = 14
    else:
        raise Exception("Error: data match type " + data_match_type + " not recognized.")

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
    in_dir = loc_root + "/MultiFreqExperiment/RatioAnalysis/data/" + station + "/" + station + year + month + day
    in_file = in_dir + "/" + station + year + month + day + "." + \
              data_match_type + "MatchedData.1gg" + str(second_resolution) + "s.pkl"
    df = pd.read_pickle(in_file)

    # Filter the data based on the expected gate range of the region of interest
    df = df.loc[(df['gate'] >= gates[0]) & (df['gate'] <= gates[1])]
    df.reset_index(drop=True, inplace=True)

    out_dir = loc_root + "/MultiFreqExperiment/RatioAnalysis/out/"

    normal_ratio_c = "k"
    ylarger_ratio_c = "b"
    xlarger_ratio_c = "r"
    pts_size = 2

    # For 4 hours, we need 8 time periods (half an hour periods)
    if SHOW_PLOTS:
        time_chunks = 1
    else:
        time_chunks = (end_hour - start_hour) * 2

    for time_chunk in range(time_chunks):
        start_time = start_hour + time_chunk * 0.5
        end_time = (start_hour + 0.5) + time_chunk * 0.5

        print("Plotting from: " + str(start_time) + " till " + str(end_time) + " UT")
        time_restricted_df = df[(df['decimalTime'] >= start_time) & (df['decimalTime'] < end_time)]

        # Build frequency dependant data frames
        if data_match_type == "Median":
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
        elif data_match_type == "Raw":
            df_10_12 = time_restricted_df[time_restricted_df['10over12'].notna()]
            df_13_12 = time_restricted_df[time_restricted_df['13over12'].notna()]
            df_14_12 = time_restricted_df[time_restricted_df['14over12'].notna()]
            df_14_13 = time_restricted_df[time_restricted_df['14over13'].notna()]
            df_13_10 = time_restricted_df[time_restricted_df['13over10'].notna()]
            df_14_10 = time_restricted_df[time_restricted_df['14over10'].notna()]
        else:
            raise Exception("Error: data match type " + data_match_type + " not recognized.")

        # Set up the plot
        n_rows = 3
        n_cols = 2
        fig1, ax1 = plt.subplots(figsize=(8, 9), dpi=300, nrows=n_rows, ncols=n_cols)
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
        fig1.suptitle("Echo Frequency Dependence: Elevation Angle Comparison"
                      + "\n" + mnemonic + " " + year + "." + month + "." + day
                      + "; " + gate_label + "; " + data_match_type + " Matched Data"
                      + "; " + str(start_time) + "-" + str(end_time) + " UT"
                      + "\nProduced by " + str(os.path.basename(__file__)),
                      fontsize=13)

        # Format subplots
        for row in range(n_rows):
            ax1[row][0].set_xlabel("12 MHz Elevation Angles [deg]")
            ax1[row][1].set_xlabel("10 MHz Elevation Angles [deg]")

            for col in range(n_cols):
                ax1[row][col].set_ylim([0, 40])
                ax1[row][col].set_xlim([0, 40])
                ax1[row][col].yaxis.set_minor_locator(MultipleLocator(10))
                ax1[row][col].xaxis.set_minor_locator(MultipleLocator(10))
                ax1[row][col].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
                ax1[row][col].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
                # ax1[row][col].plot(ax1[row][col].get_ylim(), [, 120], linestyle='-', linewidth=0.5, color='black')
                # ax1[row][col].plot([120, 120], ax1[row][col].get_xlim(), linestyle='-', linewidth=0.5, color='black')
                ax1[row][col].plot([ax1[row][col].get_ylim()[0], ax1[row][col].get_ylim()[1]],
                                   [ax1[row][col].get_xlim()[0], ax1[row][col].get_xlim()[1]],
                                   linestyle='-', linewidth=1, color='m')
                ax1[row][col].plot([ax1[row][col].get_ylim()[0] / 0.88, ax1[row][col].get_ylim()[1] / 0.88],
                                   [ax1[row][col].get_xlim()[0], ax1[row][col].get_xlim()[1]],
                                   linestyle='--', linewidth=0.5, color='m')
                ax1[row][col].plot([ax1[row][col].get_ylim()[0] / 1.12, ax1[row][col].get_ylim()[1] / 1.12],
                                   [ax1[row][col].get_xlim()[0], ax1[row][col].get_xlim()[1]],
                                   linestyle='--', linewidth=0.5, color='m')

        ax1[0][1].set_xlabel("13 MHz Elevation Angles [deg]")

        ax1[0][0].set_ylabel("10 MHz Elevation Angles [deg]")
        ax1[1][0].set_ylabel("13 MHz Elevation Angles [deg]")
        ax1[2][0].set_ylabel("14 MHz Elevation Angles [deg]")

        ax1[0][1].set_ylabel("14 MHz Elevation Angles [deg]")
        ax1[1][1].set_ylabel("13 MHz Elevation Angles [deg]")
        ax1[2][1].set_ylabel("14 MHz Elevation Angles [deg]")

        # Plot 10 to 12 Comparison data in ROW: 0, COL: 0
        ax1[0][0].scatter(df_10_12.loc[df_10_12['diffElvFlag_10largerthan12'], 'elv12'],
                          df_10_12.loc[df_10_12['diffElvFlag_10largerthan12'], 'elv10'],
                          s=pts_size, color=ylarger_ratio_c, marker='.', zorder=4, label="y larger than x")
        ax1[0][0].scatter(df_10_12.loc[df_10_12['diffElvFlag_10lessthan12'], 'elv12'],
                          df_10_12.loc[df_10_12['diffElvFlag_10lessthan12'], 'elv10'],
                          s=pts_size, color=xlarger_ratio_c, marker='.', zorder=4, label="x larger than y")
        ax1[0][0].scatter(df_10_12.loc[df_10_12['diffElvFlag_10about12'], 'elv12'],
                          df_10_12.loc[df_10_12['diffElvFlag_10about12'], 'elv10'],
                          s=pts_size, color=normal_ratio_c, marker='.', zorder=3,
                          label="Similar Elevation Angles")
        ax1[0][0].text(2, 34, 'n=' + str(df_10_12.shape[0]), fontsize=12, c='purple')

        # Plot 13 to 12 Comparison data in ROW: 1, COL: 0
        ax1[1][0].scatter(df_13_12.loc[df_13_12['diffElvFlag_13largerthan12'], 'elv12'],
                          df_13_12.loc[df_13_12['diffElvFlag_13largerthan12'], 'elv13'],
                          s=4, color=ylarger_ratio_c, marker='.', zorder=3)
        ax1[1][0].scatter(df_13_12.loc[df_13_12['diffElvFlag_13lessthan12'], 'elv12'],
                          df_13_12.loc[df_13_12['diffElvFlag_13lessthan12'], 'elv13'],
                          s=4, color=xlarger_ratio_c, marker='.', zorder=3)
        ax1[1][0].scatter(df_13_12.loc[df_13_12['diffElvFlag_13about12'], 'elv12'],
                          df_13_12.loc[df_13_12['diffElvFlag_13about12'], 'elv13'],
                          s=4, color=normal_ratio_c, marker='.', zorder=3)
        ax1[1][0].text(2, 34, 'n=' + str(df_13_12.shape[0]), fontsize=12, c='purple')

        # Plot 14 to 12 Comparison data in ROW: 2, COL: 0
        ax1[2][0].scatter(df_14_12.loc[df_14_12['diffElvFlag_14largerthan12'], 'elv12'],
                          df_14_12.loc[df_14_12['diffElvFlag_14largerthan12'], 'elv14'],
                          s=4, color=ylarger_ratio_c, marker='.', zorder=3)
        ax1[2][0].scatter(df_14_12.loc[df_14_12['diffElvFlag_14lessthan12'], 'elv12'],
                          df_14_12.loc[df_14_12['diffElvFlag_14lessthan12'], 'elv14'],
                          s=4, color=xlarger_ratio_c, marker='.', zorder=3)
        ax1[2][0].scatter(df_14_12.loc[df_14_12['diffElvFlag_14about12'], 'elv12'],
                          df_14_12.loc[df_14_12['diffElvFlag_14about12'], 'elv14'],
                          s=4, color=normal_ratio_c, marker='.', zorder=3)
        ax1[2][0].text(2, 34, 'n=' + str(df_14_12.shape[0]), fontsize=12, c='purple')

        # Plot 14 to 13 Comparison data in ROW: 0, COL: 1
        ax1[0][1].scatter(df_14_13.loc[df_14_13['diffElvFlag_14largerthan13'], 'elv13'],
                          df_14_13.loc[df_14_13['diffElvFlag_14largerthan13'], 'elv14'],
                          s=4, color=ylarger_ratio_c, marker='.', zorder=3)
        ax1[0][1].scatter(df_14_13.loc[df_14_13['diffElvFlag_14lessthan13'], 'elv13'],
                          df_14_13.loc[df_14_13['diffElvFlag_14lessthan13'], 'elv14'],
                          s=4, color=xlarger_ratio_c, marker='.', zorder=3)
        ax1[0][1].scatter(df_14_13.loc[df_14_13['diffElvFlag_14about13'], 'elv13'],
                          df_14_13.loc[df_14_13['diffElvFlag_14about13'], 'elv14'],
                          s=4, color=normal_ratio_c, marker='.', zorder=3)
        ax1[0][1].text(2, 34, 'n=' + str(df_14_13.shape[0]), fontsize=12, c='purple')

        # Plot 13 to 10 Comparison data in ROW: 1, COL: 1
        ax1[1][1].scatter(df_13_10.loc[df_13_10['diffElvFlag_13largerthan10'], 'elv10'],
                          df_13_10.loc[df_13_10['diffElvFlag_13largerthan10'], 'elv13'],
                          s=4, color=ylarger_ratio_c, marker='.', zorder=3)
        ax1[1][1].scatter(df_13_10.loc[df_13_10['diffElvFlag_13lessthan10'], 'elv10'],
                          df_13_10.loc[df_13_10['diffElvFlag_13lessthan10'], 'elv13'],
                          s=4, color=xlarger_ratio_c, marker='.', zorder=3)
        ax1[1][1].scatter(df_13_10.loc[df_13_10['diffElvFlag_13about10'], 'elv10'],
                          df_13_10.loc[df_13_10['diffElvFlag_13about10'], 'elv13'],
                          s=4, color=normal_ratio_c, marker='.', zorder=3)
        ax1[1][1].text(2, 34, 'n=' + str(df_13_10.shape[0]), fontsize=12, c='purple')

        # Plot 14 to 10 Comparison data in ROW: 2, COL: 1
        ax1[2][1].scatter(df_14_10.loc[df_14_10['diffElvFlag_14largerthan10'], 'elv10'],
                          df_14_10.loc[df_14_10['diffElvFlag_14largerthan10'], 'elv14'],
                          s=4, color=ylarger_ratio_c, marker='.', zorder=3)
        ax1[2][1].scatter(df_14_10.loc[df_14_10['diffElvFlag_14lessthan10'], 'elv10'],
                          df_14_10.loc[df_14_10['diffElvFlag_14lessthan10'], 'elv14'],
                          s=4, color=xlarger_ratio_c, marker='.', zorder=3)
        ax1[2][1].scatter(df_14_10.loc[df_14_10['diffElvFlag_14about10'], 'elv10'],
                          df_14_10.loc[df_14_10['diffElvFlag_14about10'], 'elv14'],
                          s=4, color=normal_ratio_c, marker='.', zorder=3)
        ax1[2][1].text(2, 34, 'n=' + str(df_14_10.shape[0]), fontsize=12, c='purple')

        # Add legends to the plots
        ax1[0][0].legend(loc=(0.2, 1.03), ncol=3)

        if SHOW_PLOTS:
            plt.show()

        if SAVE_PLOTS:
            # Save the files as a temp file
            fig1.savefig(out_dir + "/" + mnemonic + "_elv_comp_" + year + month + day
                         + "_" + chr(ord('a') + time_chunk) + " " + str(start_time) + "-" + str(end_time) + "UT"
                         + "_temp1.pdf", format='pdf', dpi=300)

    if SAVE_PLOTS:
        # Merge all the temp pdf files
        merger = PdfFileMerger()
        for pdf in glob.iglob("out/*_temp1.pdf"):
            merger.append(pdf)
        with open(out_dir + "/" + mnemonic + "_elv_comp_" + year + month + day +
                  "_gg" + str(gates[0]) + "-" + str(gates[1]) + "_" + data_match_type + ".pdf",
                  "wb") as fout:
            merger.write(fout)
        merger.close()

        for file in glob.iglob("out/*_temp*.pdf"):
            os.remove(file)
