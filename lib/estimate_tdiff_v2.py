import bz2
import os
import pathlib
import warnings

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from lib.elevation_v2 import elevation_v2
from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch
from lib.basic_SD_df_filter import basic_SD_df_filter


def estimate_tdiff_v2(df):
    """

    Use the percent difference method to estimate a suggested tdiff.

    # TODO: This program is the same as estimate_tdiff().  The only difference is in the helper function.  Add the
       parameter to use (phase or elevation) as a function parameter and combine these two programs.
    # TODO: It might be worth the speed increment to create a new program that just computes adjusted phase.
       elevation_v2() computes both adjusted elevation and adjusted phase.

    The program considers the percent difference between phase, instead of elevation (elevation is used in
     estimate_tdiff()).  The goal is that this will avoid problems caused by the nonlinear conversion from
     phase -> elevation.

    This is effectively just an automation of the Ponomarenko et al. (2015) technique.  For more information see:
        Application of ground scatter returns for calibration of HF interferometry data. Ponomarenko et al.
        Earth, Planets and Space (2015). 67:138 DOI 10.1186/s40623-015-0310-3.

    To perform this calculation, elevation angles are calculated for a range of tdiffs.  These columns are left in the
     dataframe and keyed as "'elv_for_' + str(tdiff)".  Therefore, the calibrated elevation angles are already
     in the dataframe, access with "'elv_for_' + str(preferred_tdiff)"

    # TODO: Right now this dataframe uses frame.insert in a loop.  This results in a highly fragmented dataframe.

    :param df: pandas.DataFrame:
            A SuperDARN fit dataframe.
            # TODO: The df might need to be restricted to certain echo bands.  This requires more investigation.

    :return suggested_tdiff, fig: float, matplotlib.pyplot.figure:
            An estimate for tdiff (in microseconds) and a supporting figure.  The tdiff estimate found can be used to
             calibrate the dataframes elevation angles (see elevation_v2() for more on elevation angle calibration).
             The  can then be modified, added to, printed out, or saved in whatever format is desired.
    """

    tdiff_range_ns = (0, 50)  # the range of t_diffs to run through, in nanoseconds
    label_fontsize = 14
    title_font_size = 18

    if len(df) <= 0:
        warnings.warn("There is no data in the provided dataframe, estimate_tdiff() is returning None, None",
                      category=Warning)
        return None, None

    # Make sure all of the data is of the same spatial resolution
    spatial_resolution = df['rangeSep'].iat[0]
    warnings.warn("estimate_tdiff() is only considering " + str(spatial_resolution) + " km data.", category=Warning)
    df = df.loc[(df['rangeSep'] == spatial_resolution)]

    # Filter out ground scatter, low quality echoes, and low power echoes
    df = basic_SD_df_filter(df)

    # Prepare a figure, we will return some supporting plots to aid in validation of the provided tdiff
    fig, axes = plt.subplots(figsize=[10, 8], dpi=300, nrows=1, ncols=2, constrained_layout=True)

    # Format the subplots
    axes[0].set_ylabel("Percent Change [%]", fontsize=label_fontsize)
    axes[1].set_ylabel("First Derivative of Percent Change", fontsize=label_fontsize)
    for ax in axes:
        ax.grid(b=True, which='both', axis='both', zorder=5)
        ax.set_xlabel("Artificial Time Delay, tdiff [ns]", fontsize=label_fontsize)  # We will plot in nanoseconds
        ax.set_xlim(np.asarray(tdiff_range_ns))

    suggested_tdiff = estimate_tdiff_helper(axes=axes, df=df, tdiff_range_ns=tdiff_range_ns)

    # title the plot
    starting_datetime = df['datetime'].iat[0]
    ending_datetime = df['datetime'].iat[-1]
    station = df['station'].iat[0]
    fig.suptitle("Event: " + str(starting_datetime) + " to " + str(ending_datetime) + " at " + station.upper() + "." +
                 "\nSuggested tdiff=" + str(round(suggested_tdiff, 1)) + " ns." +
                 "\nProduced by " + str(os.path.basename(__file__)),
                 fontsize=title_font_size)

    # Add plot legends
    for ax in axes:
        ax.legend(loc='upper right', prop={'size': 15})

    return suggested_tdiff, fig


def estimate_tdiff_helper(axes, df, tdiff_range_ns):
    """

    Helper function to run through all tdiffs in the provided range, compute percent differences, fill in the plots,
     and actually find a suggested tdiff.

    :param axes: matplotlib.axes:
            The axes to draw on.
    :param df: pandas.DataFrame:
            A SuperDARN fit dataframe.
    :param tdiff_range_ns: (int, int):
            The range of tdiffs to consider in nanoseconds.  This will be used as the x-axis limits

    :return tdiff: float:
            The suggested tdiff, in microseconds.
    """

    # When we smooth the data, the edges are going to become unusable, extend the range to compensate for this
    extended_tdiff_range_ns = (tdiff_range_ns[0] - 5, tdiff_range_ns[1] + 5)

    # Compute elevation angles for each tdiff in the extended range
    num = int((extended_tdiff_range_ns[1] - extended_tdiff_range_ns[0] + 1) * 3)
    tdiffs_ns = np.linspace(start=extended_tdiff_range_ns[0], stop=extended_tdiff_range_ns[1], num=num)

    for tdiff in tdiffs_ns:
        print("Computing elevation angles for tdiff=" + str(round(tdiff, 2)) + " ns.")
        elevation_v2(df=df, t_diff=tdiff / 1000)
        df["phase_for_" + str(tdiff)] = df['adjPhase']

    # Compute the mean percent change between elevation angles resulting from adjacent tdiffs
    mean_percent_change = []
    for t in range(len(tdiffs_ns) - 1):
        tdiff_init = tdiffs_ns[t]
        tdiff_final = tdiffs_ns[t + 1]
        elv_init_key = "phase_for_" + str(tdiff_init)
        elv_final_key = "phase_for_" + str(tdiff_final)

        mean_percent_change.append(np.mean((df[elv_final_key] - df[elv_init_key]) / df[elv_init_key]) * 100)

    # There is one less point in mean_percent_change than there is in tdiffs.
    # Make a tdiff array that has the same number of points so we can plot the data
    delta = tdiffs_ns[1] - tdiffs_ns[0]
    tdiffs_for_plotting = tdiffs_ns[1:] - delta / 2

    # Smooth the signal
    kernel_size = 5
    kernel = np.ones(kernel_size) / kernel_size
    smoothed_mean_percent_change = np.convolve(mean_percent_change, kernel, mode='same')

    # Lightly plot the raw data
    axes[0].plot(tdiffs_for_plotting, mean_percent_change, color="red", linewidth=0.2, label="Raw")
    # Plot the smoothed signal
    axes[0].plot(tdiffs_for_plotting, smoothed_mean_percent_change, color="blue", linewidth=1, label="Smooth")

    # Compute the derivative of the smoothed signal
    mean_percent_change_diff = np.diff(smoothed_mean_percent_change)

    # The derivative is one element smaller, make a tdiff array that has the same number of points so we can plot
    tdiffs_for_plotting_diff = tdiffs_for_plotting[1:] - delta / 2

    # Smooth the derivative
    smoothed_mean_percent_change_diff = np.convolve(mean_percent_change_diff, kernel, mode='same')

    # Find the minimum of this derivative, it will be our suggested tdiff
    # Make sure we are not way off base be restricting the candidates
    in_range_candidates = np.where((tdiffs_for_plotting_diff > tdiff_range_ns[0]) &
                                   (tdiffs_for_plotting_diff < tdiff_range_ns[1]))
    min_candidates_y = smoothed_mean_percent_change_diff[in_range_candidates]
    min_candidates_x = tdiffs_for_plotting_diff[in_range_candidates]

    min_y = np.min(min_candidates_y)
    x_index = np.where(min_candidates_y == min_y)
    min_x = min_candidates_x[x_index]

    # Lightly plot the raw data
    axes[1].plot(tdiffs_for_plotting_diff, mean_percent_change_diff, color="red", linewidth=0.2, label="Raw")
    # Plot the smoothed signal
    axes[1].plot(tdiffs_for_plotting_diff, smoothed_mean_percent_change_diff, color="blue", linewidth=1,
                 label="Smooth")
    # Plot the candidate value
    axes[1].plot(min_x, min_y, 'go', label="Suggested tdiff")

    return min_x[0]


if __name__ == "__main__":
    """ Testing """

    station = "hok"
    year = "2014"
    month = "02"
    day = "23"

    time_units = 'ut'
    start_hour = 6
    end_hour = 12

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()))
    in_dir = loc_root + "/DataAnalysis/DataReading/SD/data/" + station + "/" + station + year + month + day
    in_file = in_dir + "/" + station + year + month + day + ".pbz2"
    data_stream = bz2.BZ2File(in_file, "rb")
    df = pd.read_pickle(data_stream)

    # Restrict the data to within the desired hour range / beam
    _, start_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=start_hour)
    _, end_epoch = build_datetime_epoch(year=int(year), month=int(month), day=int(day), hour=end_hour)
    df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]
    df = df.loc[df['bmnum'] == 7]
    df.reset_index(drop=True, inplace=True)

    tdiff, fig = estimate_tdiff_v2(df=df)

    print("The suggested tdiff is: " + str(tdiff) + " \u03BCs")
    plt.show()

    # Save the figure to file
    loc_root = str(pathlib.Path().parent.absolute())
    out_dir = loc_root + "/out"
    out_fig = out_dir + "/estimate_tidff_v2_" + station + year + month + day
    print("Saving plot as " + out_fig)
    fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)

