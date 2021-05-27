import glob
import os
import pathlib

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.cm as cm

from matplotlib.colors import ListedColormap
from matplotlib.ticker import MultipleLocator
from PyPDF2 import PdfFileMerger
from scipy import stats

if __name__ == '__main__':
    """
    Plot velocity comparison (as contour) beside elevation comparison (also as contour)
    """

    SAVE_PLOTS = True
    SHOW_PLOTS = False

    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd

    station = "rkn"
    gates = [10, 30]
    data_match_type = "Raw"  # "Matched" or "Raw"
    count_min = 4  # Only used for median matched data

    start_hour = 0
    end_hour = 4

    show_scatter = False

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
    in_dir = loc_root + "/MultiFreqExperiment/VelocityAnalysis/data/" + station + "/" + station + year + month + day
    in_file = in_dir + "/" + station + year + month + day + "." + \
              data_match_type + "MatchedData.1gg" + str(second_resolution) + "s.pkl"
    df = pd.read_pickle(in_file)

    # Filter the data based on the expected gate range of the region of interest
    df = df.loc[(df['gate'] >= gates[0]) & (df['gate'] <= gates[1])]
    df.reset_index(drop=True, inplace=True)

    out_dir = loc_root + "/MultiFreqExperiment/VelocityAnalysis/out/"

    height_scat_c = 'k'
    vel_scat_c = 'k'
    pts_size = 3

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
        fig1.suptitle("Echo Frequency Dependence: Velocity Comparison w Elevation Angle"
                      + "\n" + mnemonic + " " + year + "." + month + "." + day
                      + "; " + gate_label + "; " + data_match_type + " Matched Data"
                      + "; " + str(start_time) + "-" + str(end_time) + " UT"
                      + "\nProduced by " + str(os.path.basename(__file__)),
                      fontsize=13)

        # Format subplots
        for row in range(n_rows):
            # First row will be velocity contour plots
            ax1[row][0].set_xlabel("12 MHz Velocities [m/s]")
            ax1[row][0].set_ylim([-600, 600])
            ax1[row][0].set_xlim([-600, 600])
            ax1[row][0].yaxis.set_minor_locator(MultipleLocator(100))
            ax1[row][0].xaxis.set_minor_locator(MultipleLocator(100))
            ax1[row][0].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
            ax1[row][0].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
            ax1[row][0].plot(ax1[row][0].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
            ax1[row][0].plot([0, 0], ax1[row][0].get_xlim(), linestyle='-', linewidth=0.5, color='black')
            ax1[row][0].plot([ax1[row][0].get_ylim()[0], ax1[row][0].get_ylim()[1]],
                             [ax1[row][0].get_xlim()[0], ax1[row][0].get_xlim()[1]],
                             linestyle='--', linewidth=0.75, color='red')
            plt.sca(ax1[row][0])
            plt.xticks([-600, -400, -200, 0, 200, 400, 600])
            ax1[row][0].tick_params(axis='x', labelrotation=25)

            # Second row will be height comparison plots
            ax1[row][1].set_xlabel("12 MHz Elevation Angles [deg]")
            ax1[row][1].set_ylim([0, 40])
            ax1[row][1].set_xlim([0, 40])
            ax1[row][1].yaxis.set_minor_locator(MultipleLocator(10))
            ax1[row][1].xaxis.set_minor_locator(MultipleLocator(10))
            ax1[row][1].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
            ax1[row][1].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
            ax1[row][1].plot([ax1[row][1].get_ylim()[0], ax1[row][1].get_ylim()[1]],
                             [ax1[row][1].get_xlim()[0], ax1[row][1].get_xlim()[1]],
                             linestyle='--', linewidth=0.75, color='red')


        ax1[0][0].set_ylabel("10 MHz Velocities [m/s]")
        ax1[1][0].set_ylabel("13 MHz Velocities [m/s]")
        ax1[2][0].set_ylabel("14 MHz Velocities [m/s]")

        ax1[0][1].set_ylabel("10 MHz Elevation Angles [deg]")
        ax1[1][1].set_ylabel("13 MHz Elevation Angles [deg]")
        ax1[2][1].set_ylabel("14 MHz Elevation Angles [deg]")

        # Compute binned velocity counts
        n_bins_vel = 48  # 25 m/s bins
        contour_range_vel = [ax1[0][0].get_ylim(), ax1[0][0].get_xlim()]
        try:
            binned_counts_10_12_vel, bin_xedges_vel, bin_yedges_vel, bin_numbers_vel = stats.binned_statistic_2d(
                df_10_12['vel12'], df_10_12['vel10'], values=None,
                statistic='count', bins=[n_bins_vel, n_bins_vel], range=contour_range_vel)
        except:
            pass

        try:
            binned_counts_13_12_vel, i, ii, iii = stats.binned_statistic_2d(
                df_13_12['vel12'], df_13_12['vel13'], values=None,
                statistic='count', bins=[n_bins_vel, n_bins_vel], range=contour_range_vel)
        except:
            pass

        try:
            binned_counts_14_12_vel, j, jj, jjj = stats.binned_statistic_2d(
                df_14_12['vel12'], df_14_12['vel14'], values=None,
                statistic='count', bins=[n_bins_vel, n_bins_vel], range=contour_range_vel)
        except:
            pass

        # Compute bin centers, these will be the same for all frequency comparisons
        try:
            bin_xwidth_vel = (bin_xedges_vel[1] - bin_xedges_vel[0])
            bin_ywidth_vel = (bin_yedges_vel[1] - bin_yedges_vel[0])
            bin_xcenters_vel = bin_xedges_vel[1:] - bin_xwidth_vel / 2
            bin_ycenters_vel = bin_yedges_vel[1:] - bin_ywidth_vel / 2
        except:
            continue

        # Compute binned height counts
        n_bins_h = 20  # 2 deg bins
        contour_range_h = [ax1[0][1].get_ylim(), ax1[0][1].get_xlim()]
        try:
            binned_counts_10_12_elv, bin_xedges_elv, bin_yedges_elv, bin_numbers_elv = stats.binned_statistic_2d(
                df_10_12['elv12'], df_10_12['elv10'], values=None,
                statistic='count', bins=[n_bins_h, n_bins_h], range=contour_range_h)
        except:
            pass

        try:
            binned_counts_13_12_elv, k, kk, kkk = stats.binned_statistic_2d(
                df_13_12['elv12'], df_13_12['elv13'], values=None,
                statistic='count', bins=[n_bins_h, n_bins_h], range=contour_range_h)
        except:
            pass

        try:
            binned_counts_14_12_elv, l, ll, lll = stats.binned_statistic_2d(
                df_14_12['elv12'], df_14_12['elv14'], values=None,
                statistic='count', bins=[n_bins_h, n_bins_h], range=contour_range_h)
        except:
            pass

        # Compute bin centers, these will be the same for all frequency comparisons
        bin_xwidth_elv = (bin_xedges_elv[1] - bin_xedges_elv[0])
        bin_ywidth_elv = (bin_yedges_elv[1] - bin_yedges_elv[0])
        bin_xcenters_elv = bin_xedges_elv[1:] - bin_xwidth_elv / 2
        bin_ycenters_elv = bin_yedges_elv[1:] - bin_ywidth_elv / 2

        # Modify the 'jet' colour map
        jet = cm.get_cmap('jet', 256)
        newcolours = jet(np.linspace(0, 1, 256))
        purple = np.array([155/256, 53/256, 161/256, 1])  # RGBA colours
        white = np.array([0, 0, 0, 0])
        newcolours[:35, :] = white  # Make the first few colours white
        newcolours[220:255, :] = purple  # Make the last few colours purple
        newcmp = ListedColormap(newcolours)

        # Plot 10 to 12 Comparison data in ROW: 0
        cont = ax1[0][0].contourf(bin_xcenters_vel, bin_ycenters_vel, binned_counts_10_12_vel.transpose(),
                                  5, cmap=newcmp)
        fig1.colorbar(cont, ax=ax1[0][0])
        cont = ax1[0][1].contourf(bin_xcenters_elv, bin_ycenters_elv, binned_counts_10_12_elv.transpose(),
                                  5, cmap=newcmp)
        fig1.colorbar(cont, ax=ax1[0][1])
        if show_scatter:
            ax1[0][0].scatter(df_10_12['vel12'], df_10_12['vel10'],
                              s=pts_size, color=height_scat_c, marker='.', zorder=3)
            ax1[0][1].scatter(df_10_12['elv12'], df_10_12['elv10'],
                              s=pts_size, color=height_scat_c, marker='.', zorder=3)
        ax1[0][0].text(-490, 425, 'n=' + str(df_10_12.shape[0]), fontsize=12, c='k')
        ax1[0][1].text(2, 34, 'n=' + str(df_10_12.shape[0]), fontsize=12, c='k')

        # Plot 13 to 12 Comparison data in ROW: 1
        cont = ax1[1][0].contourf(bin_xcenters_vel, bin_ycenters_vel, binned_counts_13_12_vel.transpose(),
                                  5, cmap=newcmp)
        fig1.colorbar(cont, ax=ax1[1][0])
        cont = ax1[1][1].contourf(bin_xcenters_elv, bin_ycenters_elv, binned_counts_13_12_elv.transpose(),
                                  5, cmap=newcmp)
        fig1.colorbar(cont, ax=ax1[1][1])
        if show_scatter:
            ax1[1][0].scatter(df_13_12['vel12'], df_13_12['vel13'],
                              s=pts_size, color=height_scat_c, marker='.', zorder=3)
            ax1[1][1].scatter(df_13_12['elv12'], df_13_12['elv13'],
                              s=pts_size, color=height_scat_c, marker='.', zorder=3)
        ax1[1][0].text(-490, 425, 'n=' + str(df_13_12.shape[0]), fontsize=12, c='k')
        ax1[1][1].text(2, 34, 'n=' + str(df_13_12.shape[0]), fontsize=12, c='k')

        # Plot 14 to 12 Comparison data in ROW: 2
        cont = ax1[2][0].contourf(bin_xcenters_vel, bin_ycenters_vel, binned_counts_14_12_vel.transpose(),
                                  5, cmap=newcmp)
        fig1.colorbar(cont, ax=ax1[2][0])
        cont = ax1[2][1].contourf(bin_xcenters_elv, bin_ycenters_elv, binned_counts_14_12_elv.transpose(),
                                  5, cmap=newcmp)
        fig1.colorbar(cont, ax=ax1[2][1])
        if show_scatter:
            ax1[2][0].scatter(df_14_12['vel12'], df_14_12['vel14'],
                              s=pts_size, color=height_scat_c, marker='.', zorder=3)
            ax1[2][1].scatter(df_14_12['elv12'], df_14_12['elv14'],
                              s=pts_size, color=height_scat_c, marker='.', zorder=3)
        ax1[2][0].text(-490, 425, 'n=' + str(df_14_12.shape[0]), fontsize=12, c='k')
        ax1[2][1].text(2, 34, 'n=' + str(df_14_12.shape[0]), fontsize=12, c='k')

        if SHOW_PLOTS:
            plt.show()

        if SAVE_PLOTS:
            # Save the files as a temp file
            fig1.savefig(out_dir + "/" + mnemonic + "_velContourWithElv_" + year + month + day
                         + "_" + chr(ord('a') + time_chunk) + " " + str(start_time) + "-" + str(end_time) + "UT"
                         + "_temp.pdf", format='pdf', dpi=300)

    if SAVE_PLOTS:
        # Merge all the temp pdf files
        merger = PdfFileMerger()
        for pdf in glob.iglob("out/*_temp.pdf"):
            merger.append(pdf)
        with open(out_dir + "/" + mnemonic + "_velContourWithElv_" + year + month + day
                  + "_gg" + str(gates[0]) + "-" + str(gates[1]) + "_" + data_match_type + ".pdf",
                  "wb") as fout:
            merger.write(fout)
        merger.close()

        for file in glob.iglob("out/*_temp*.pdf"):
            os.remove(file)
