import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.ticker import MultipleLocator
from lib.basic_SD_df_filter import basic_SD_df_filter
from DataAnalysis.MultiFreqExperiment.get_df_multi_event import get_df_multi_event

if __name__ == '__main__':
    """
    Produce a velocity histogram that includes data for several events
    """

    PLOT_OUT = "show"  # "save" or "show"

    list_of_events = "event_summary"
    list_of_events_flag = "for_vel_hist"

    max_height = 160  # The maximum allowed height, points coming from above this height will be assumed to be F
    #  region scatter and will not be included

    Re = 6370  # Radius of the Earth, [km]
    gate_range = [10, 30]  # Inclusive
    frequencies = [10, 12, 13, 14]  # MHz

    gate_label = "gg: " + str(gate_range[0]) + "-" + str(gate_range[1])
    vel_lim = [-600, 600]
    occ_lim = [-100, 8600]
    n_bins = int((vel_lim[1] - vel_lim[0]) / 50)  # 50 m/s bins

    # Get multi-event data
    print("Getting multi-event data...")
    df = get_df_multi_event(file_name=list_of_events, flag=list_of_events_flag, include_adj_elv=True)

    print("Filtering data...")

    # We are only interested in 15 km resolution data (from the multi-freq analysis)
    # This double filter should be redundant but better to be safe
    df = basic_SD_df_filter(df)
    df = df.loc[(df['firstRang'] == 90) & (df['rangeSep'] == 15)]
    df = df.loc[(df['gate'] >= gate_range[0]) & (df['gate'] <= gate_range[1])]  # Filter for a gate range
    df = df.loc[df['vel'].notna()]  # There should not be any nan velocities, but just to be safe
    df.reset_index(drop=True, inplace=True)

    # Compute virtual height
    slant_range = 90 + 15 * df['gate'] + 15 / 2  # Slant range [km]
    df['height'] = np.sqrt(Re * Re + slant_range * slant_range
                           + 2 * Re * slant_range * np.sin(np.radians(np.asarray(df['adjElv'])))) - Re
    # Drop everything that is above the max height (this will be assumed to be F region scatter)
    df = df.loc[(df['height'] <= max_height)]
    df.reset_index(drop=True, inplace=True)

    # Put frequencies in MHz and round, this makes them easier to compare
    df['transFreq'] = round(df['transFreq'] * 1e-3, 0)
    df = df.loc[np.isin(df['transFreq'], frequencies)]  # drop any rows with unrecognized frequencies
    df.reset_index(drop=True, inplace=True)

    # Build an event summary that can be printed onto the plot
    df['eventString'] = df['station'] + " " + \
                        df['year'].astype(str) + "." + df['month'].astype(str) + "." + df['day'].astype(str)
    list_of_included_events = [i for i in df['eventString'].unique()]

    # Set up the plot
    n_rows = len(frequencies)
    n_cols = 1
    fig, ax = plt.subplots(figsize=(8, 9), sharex='col', dpi=300, nrows=n_rows, ncols=n_cols)
    plt.subplots_adjust(hspace=0.05)
    fig.suptitle("Velocity Frequency Dependence: Multi-Event Occurrence Distributions"
                 + "\n" + gate_label + "; Max Virtual Height: " + str(max_height) + " km"
                 + "\n Events: " + str(list_of_included_events)
                 + "\nProduced by " + str(os.path.basename(__file__)),
                 fontsize=12)

    # start by just plotting the first event
    first_event = list_of_included_events[0]
    event_df = df[df['eventString'] == first_event]

    # Plot the histograms
    for row in range(n_rows):
        freq_df = event_df[event_df['transFreq'] == frequencies[row]]
        ax[row].hist(freq_df['vel'], bins=n_bins, range=vel_lim, align='mid', histtype='step', zorder=3, label=first_event)

    # Loop thorough all the other events and add them on one at a time
    for event in list_of_included_events:
        if event == first_event:
            continue  # we have already done this one

        # add in the next event
        event_df = pd.concat([event_df, df[df['eventString'] == event]])

        # Plot the histogram
        for row in range(n_rows):
            freq_df = event_df[event_df['transFreq'] == frequencies[row]]
            ax[row].hist(freq_df['vel'], bins=n_bins, range=vel_lim, align='mid', histtype='step', zorder=3, label=event)

    # Print out frequency and count info
    for row in range(n_rows):
        freq_df = df[df['transFreq'] == frequencies[row]]
        ax[row].text(vel_lim[1] - 0.12 * (vel_lim[1] - vel_lim[0]),
                        occ_lim[0] + 0.80 * (occ_lim[1] - occ_lim[0]), "n=" + str(freq_df.shape[0]))
        ax[row].text(vel_lim[1] - 0.12 * (vel_lim[1] - vel_lim[0]),
                        occ_lim[0] + 0.88 * (occ_lim[1] - occ_lim[0]), str(frequencies[row]) + " MHz")

    # Format the sub-plots
    ax[n_rows - 1].set_xlabel("Velocity [m/s]")
    ax[n_rows - 1].set_ylabel("Occurrence")
    for row in range(n_rows):
        ax[row].set_xlim(vel_lim)
        ax[row].set_ylim(occ_lim)
        ax[row].xaxis.set_minor_locator(MultipleLocator(50))
        ax[row].grid(b=True, which='major', axis='x', linestyle='--', linewidth=0.5)
        ax[row].grid(b=True, which='minor', axis='x', linestyle='--', linewidth=0.1)
        ax[row].plot([0, 0], ax[row].get_ylim(), linestyle='-', linewidth=0.5, color='black')
        # ax[row][0].plot([-300, -300], ax[row][0].get_ylim(), linestyle='--', linewidth=0.5, color='black')
        # ax[row][0].plot([300, 300], ax[row][0].get_ylim(), linestyle='--', linewidth=0.5, color='black')
        ax[row].tick_params(direction="in")

    # Add a legend to the first plot
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(handles[::-1], labels[::-1], loc='upper left')

    if PLOT_OUT == "show":
        plt.show()

    if PLOT_OUT == "save":
        # Save to file
        loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
        out_dir = loc_root + "/MultiFreqExperiment/OccRates/out/"
        fig.savefig(out_dir + "/multi_event_occ_dist_of_vel_"
                    + "_gg" + str(gate_range[0]) + "-" + str(gate_range[1])
                    + ".pdf", format='pdf', dpi=300)
