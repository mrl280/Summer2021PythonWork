import calendar
import math
import os
import time
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.colors import ListedColormap

from scipy import stats
from DataAnalysis.DataReading.SD.basic_SD_df_filter import basic_SD_df_filter
from DataAnalysis.DataReading.SD.elevation_v2 import elevation_v2

if __name__ == '__main__':
    """
    Produces a file with 2 types of plots:
        In the first column, plot occurrence rate against velocity
        In the second column, plot range vs velocity with occurrence rates as contours
    """

    PLOT_OUT = "show"  # "save" or "show"

    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd

    station = "rkn"
    start_hour = 0
    end_hour = 4

    max_height = 180  # The maximum allowed height, points coming from above this height will be assumed to be F
    #  region scatter and will not be included

    t_diff = 0.003  # Elevation angle correction in microseconds
    Re = 6370  # Radius of the Earth, [km]
    gate_range = [10, 30]  # Inclusive
    frequencies = [10, 12, 13, 14]

    mnemonic = station.upper()
    gate_label = "gg: " + str(gate_range[0]) + "-" + str(gate_range[1])
    vel_lim = [-600, 600]
    occ_lim = [-100, 2600]
    n_bins = int((vel_lim[1] - vel_lim[0]) / 50)  # 50 m/s bins

    # Compute start and end epoch
    pattern = '%Y-%m-%d %H:%M:%S'
    start_date_time = year + "-" + month + "-" + day + " " + str(start_hour) + ":00:00"
    end_date_time = year + "-" + month + "-" + day + " " + str(end_hour) + ":00:00"
    start_epoch = calendar.timegm(time.strptime(start_date_time, pattern))
    end_epoch = calendar.timegm(time.strptime(end_date_time, pattern))

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
    in_dir = loc_root + "/DataReading/SD/data/" + station + "/" + station + year + month + day
    in_file = in_dir + "/" + station + year + month + day + ".pkl"
    df = pd.read_pickle(in_file)

    # We are only interested in 15 km resolution data (from the multi-freq analysis)
    # This double filter should be redundant but better to be safe
    df = df.loc[(df['firstRang'] == 90) & (df['rangeSep'] == 15)]

    df = basic_SD_df_filter(df)

    # Filter for a gate range and time
    df = df.loc[(df['gate'] >= gate_range[0]) & (df['gate'] <= gate_range[1])]
    df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]

    df = df.loc[df['vel'].notna()]  # There should not be any nan velocities, but just to be safe

    # Compute adjusted elevation angle and use it to compute virtual height
    df.reset_index(drop=True, inplace=True)
    elevation_v2(df, t_diff)
    slant_range = 90 + 15 * df['gate'] + 15 / 2  # Slant range [km]
    df['height'] = np.sqrt(Re * Re + slant_range * slant_range
                           + 2 * Re * slant_range * np.sin(np.radians(np.asarray(df['adjElv'])))) - Re

    # Drop everything that is above the max height (this will be assumed to be F region scatter)
    df = df.loc[(df['height'] <= max_height)]
    df.reset_index(drop=True, inplace=True)

    # Put frequencies in MHz and round, this makes them easier to compare, drop any rows with unrecognized frequencies
    df['transFreq'] = round(df['transFreq'] * 1e-3, 0)
    df = df.loc[np.isin(df['transFreq'], frequencies)]
    df.reset_index(drop=True, inplace=True)

    # Set up the plot
    n_rows = len(frequencies)
    n_cols = 2
    fig, ax = plt.subplots(figsize=(8, 9), sharex='col', dpi=300, nrows=n_rows, ncols=n_cols)
    plt.subplots_adjust(hspace=0.05, wspace=0.4)
    fig.suptitle("Velocity Frequency Dependence: Occurrence Distributions"
                 + "\n" + mnemonic + " " + year + "." + month + "." + day
                 + ";  " + gate_label
                 + ";  " + str(start_hour) + "-" + str(end_hour) + " UT"
                 + "\nProduced by " + str(os.path.basename(__file__)),
                 fontsize=13)

    # Plot occurrence data in the first column
    for row in range(n_rows):
        freq_df = df[df['transFreq'] == frequencies[row]]
        ax[row][0].hist(freq_df['vel'], bins=n_bins, range=vel_lim, align='mid', histtype='step', zorder=3)
        # ax[row][0].text(vel_lim[0] + 0.05 * (vel_lim[1] - vel_lim[0]),
        #                 occ_lim[0] + 0.5 * (occ_lim[1] - occ_lim[0]), "f=" + str(frequencies[row]))
        ax[row][0].text(vel_lim[0] + 0.05 * (vel_lim[1] - vel_lim[0]),
                        occ_lim[0] + 0.85 * (occ_lim[1] - occ_lim[0]), "n=" + str(freq_df.shape[0]))
        print(ax[row][0].get_ylim()[0])

    # Format the first column
    ax[n_rows - 1][0].set_xlabel("Velocity [m/s]")
    for row in range(n_rows):
        ax[row][0].set_xlim(vel_lim)
        ax[row][0].set_ylim(occ_lim)
        ax[row][0].set_ylabel(str(frequencies[row]) + " MHz Occurrence")
        ax[row][0].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
        ax[row][0].plot([0, 0], ax[row][0].get_ylim(), linestyle='-', linewidth=0.5, color='black')
        # ax[row][0].plot([-300, -300], ax[row][0].get_ylim(), linestyle='--', linewidth=0.5, color='black')
        # ax[row][0].plot([300, 300], ax[row][0].get_ylim(), linestyle='--', linewidth=0.5, color='black')
        ax[row][0].tick_params(direction="in")

    # Modify the 'jet' colour map
    jet = cm.get_cmap('jet', 256)
    newcolours = jet(np.linspace(0, 1, 256))
    purple = np.array([155 / 256, 53 / 256, 161 / 256, 1])  # RGBA colours
    white = np.array([0, 0, 0, 0])
    newcolours[:35, :] = white  # Make the first few colours white
    newcolours[220:255, :] = purple  # Make the last few colours purple
    newcmp = ListedColormap(newcolours)

    # Plot something else in the second column
    for row in range(n_rows):
        freq_df = df[df['transFreq'] == frequencies[row]]
        freq_df = freq_df.loc[(freq_df['vel'] <= 600) & (freq_df['vel'] >= -600)]
        contour_range = [[start_hour, end_hour], gate_range]
        binned_counts, bin_xedges, bin_yedges, bin_numbers = stats.binned_statistic_2d(
            freq_df['decimalTime'], freq_df['gate'], values=freq_df['vel'],
            statistic='median', bins=[80, 20], range=contour_range)

        # TODO: Finish fixing up this contour plot

        pixel = ax[row][1].imshow(np.flip(binned_counts.transpose(), axis=0), aspect='auto', cmap="seismic_r",
                                  extent=(start_hour, end_hour, gate_range[0], gate_range[1]))
        fig.colorbar(pixel, ax=ax[row][1])

    # Format the second column
    ax[n_rows - 1][1].set_xlabel("Time, UT")
    ax[n_rows - 1][1].set_ylabel("Range Gate")
    for row in range(n_rows):
        # ax[row][1].set_xlim([start_hour, end_hour])
        ax[row][1].set_ylim(gate_range)

        ax[row][1].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
        ax[row][1].plot([0, 0], ax[row][0].get_ylim(), linestyle='-', linewidth=0.5, color='black')
        # ax[row][0].plot([-300, -300], ax[row][0].get_ylim(), linestyle='--', linewidth=0.5, color='black')
        # ax[row][0].plot([300, 300], ax[row][0].get_ylim(), linestyle='--', linewidth=0.5, color='black')
        # ax[row][0].tick_params(direction="in")


    if PLOT_OUT == "show":
        plt.show()
