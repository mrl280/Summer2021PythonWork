import os
import pathlib

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator
from lib.basic_SD_df_filter import basic_SD_df_filter
from DataAnalysis.OneAndOneHalfHop.get_df_multi_event import get_df_multi_event


def SD_velocity_histogram(file_name, flag, plot_out):
    """
    Create a velocity histogram
    :param plot_out: str: show or save, what to do with the build plot
    :param file_name: str: The name of the pickled data frame you would like to use.  e.g. "event_summary"
    :param flag: str: All events in the event summary that have this flag will be considered
    :return:
    """

    if plot_out != 'show' and plot_out != 'save':
        raise Exception("Invalid plot_out, options are 'show' and 'save'.")

    vel_lim = [-600, 600]
    occ_lim = [-100, 8600]
    n_bins = int((vel_lim[1] - vel_lim[0]) / 40)  # 50 m/s bins

    # Check to make sure inputs are okay
    if flag is None or file_name is None:
        raise Exception("Error: you must pass SDvelocityHistogram() a file_name and flag "
                        "so it knows which events to consider.")
    if not isinstance(flag, str) or not isinstance(file_name, str):
        raise Exception("Error: both file_name and flag must be strings.")

    df = get_df_multi_event(file_name, flag)
    df = basic_SD_df_filter(df)

    # Restrict to the overlapping data
    # At RKN the overlap is beam 5 and gates 25-32
    # At INV the overlap is beam 12 and gates 30-38
    starting_gate = np.where(df['station'] == 'rkn', 25, 30)
    ending_gate = np.where(df['station'] == 'rkn', 32, 38)
    overlapping_beam = np.where(df['station'] == 'rkn', 5, 12)
    df = df.loc[(df['gate'] >= starting_gate) & (df['gate'] <= ending_gate) & (df['bmnum'] == overlapping_beam)]
    df = df.loc[df['vel'].notna()]  # There should not be any nan velocities, but just to be safe
    df.reset_index(drop=True, inplace=True)

    # Build an event summary that can be printed onto the plot
    df['eventString'] = df['station'] + " " + \
                        df['year'].astype(str) + "." + df['month'].astype(str) + "." + df['day'].astype(str)
    list_of_included_events = [i for i in df['eventString'].unique()]

    fig, ax = plt.subplots(figsize=(8, 9), sharex='col', dpi=300, nrows=1, ncols=1)
    plt.subplots_adjust(hspace=0.05)
    fig.suptitle("One-And-One-Half-Hop: Multi-Event Occurrence Distributions"
                 + "\n Events: " + str(list_of_included_events)
                 + "\nProduced by " + str(os.path.basename(__file__)),
                 fontsize=12)

    # start by just plotting the first event
    first_event = list_of_included_events[0]
    event_df = df[df['eventString'] == first_event]
    ax.hist(event_df['vel'], bins=n_bins, range=vel_lim, align='mid', histtype='step', zorder=3, label=first_event)

    # Loop thorough all the other events and add them on one at a time
    for event in list_of_included_events:
        if event == first_event:
            continue  # we have already done this one

        # add in the next event
        event_df = pd.concat([event_df, df[df['eventString'] == event]])

        # Plot the histogram
        ax.hist(event_df['vel'], bins=n_bins, range=vel_lim, align='mid', histtype='step', zorder=3, label=event)

    # Format the plot
    ax.set_xlim(vel_lim)
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.grid(b=True, which='major', axis='x', linestyle='--', linewidth=0.5)
    ax.grid(b=True, which='minor', axis='x', linestyle='--', linewidth=0.1)
    # ax.plot([0, 0], ax.get_ylim(), linestyle='-', linewidth=0.5, color='black')

    # Add a legend to the plot
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper left')

    if plot_out == "show":
        plt.show()

    if plot_out == "save":
        # Save to file
        loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).absolute())
        out_dir = loc_root + "/VelocityAnalysis/out/"
        fig.savefig(out_dir + "/multi_event_vel_hist_" + ".pdf", format='pdf', dpi=300)


if __name__ == '__main__':
    """
    """

    # Create a velocity histogram for a set of flagged events
    SD_velocity_histogram(file_name="event_summary", flag='for_vel_hist', plot_out='save')
