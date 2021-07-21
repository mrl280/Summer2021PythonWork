import os
import pathlib
import bz2
import warnings

import _pickle as cPickle
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import MultipleLocator
from scipy import stats

from DataAnalysis.EchoOccurrence.lib.cm.modified_jet import modified_jet


def vel_and_height_contours(single_day_df, beam_range, gate_range, hour_range, count_min=4, area=None):
    """

    Plot velocity and height comparisons beside each other, both as contours (only against 12)
    This program runs for all the data in an event, for half hour intervals use plottingVelAsContwHeightComparison.py

    Note: If point count is greater than points you see it is because the rest are out of frame.

    :param single_day_df: pandas.Data_Frame:
            A single days worth of data that we want to range profile
    :param gate_range: (int, int):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (int, int):
             Just for the figure title.
    :param hour_range: (int, int) (optional):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not quite inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param count_min: int (optional):
            The minimum number of points required to count a median match.
            This parameter has no effect on raw matched data
    :param area: int:
            The numbered area of interest.

    :return: matplotlib.pyplot.figure:
            The figure, it can then be viewed, modified, or saved to file
    """

    n_levels = 5

    if len(single_day_df) <= 0:
        warnings.warn("There is no data in the provided dataframe", category=Warning)

    df = single_day_df.copy()  # Convenience

    # Restrict data to within the desired hour range
    df = df.loc[(df['decimalTime'] >= hour_range[0]) & (df['decimalTime'] <= hour_range[1])]

    # Filter the data based on the expected gate range of the region of interest
    df = df.loc[(df['gate'] >= gate_range[0]) & (df['gate'] <= gate_range[1])]

    df.reset_index(drop=True, inplace=True)

    print("Setting up the figure...")

    resolution_string = "15 km data"
    if beam_range[0] == beam_range[1]:
        beam_string = "Beam " + str(beam_range[0])
    else:
        beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    gate_string = "Gates " + str(gate_range[0]) + "-" + str(gate_range[1])
    hour_string = "Hours " + str(hour_range[0]) + "-" + str(hour_range[1]) + " UT"

    date_string = "September 26, 2016"
    title_font_size = 14

    fig, axes = plt.subplots(figsize=[8, 9], nrows=3, ncols=2, constrained_layout=True, dpi=300)
    format_subplots(axes)

    if area is None:
        fig.suptitle(date_string + " at " + station.upper() + "; " + beam_string + "; " + gate_string
                     + "; " + data_match_type + " Matched Data"
                     + "\n" + resolution_string + "; " + hour_string + "; "
                     + "Produced by " + str(os.path.basename(__file__)),
                     fontsize=title_font_size)
    else:
        fig.suptitle(date_string + " at " + station.upper() + "; Area " + str(area)
                     + "; " + data_match_type + " Matched Data"
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)),
                     fontsize=title_font_size)

    print("Building frequency dependent dataframes...")

    # Build frequency dependant data frames
    if data_match_type == "Median":
        df_10_12 = df[(df['count10'] >= count_min) & (df['count12'] >= count_min) & (df['10over12'].notna())]

        df_13_12 = df[(df['count13'] >= count_min) & (df['count12'] >= count_min) & (df['13over12'].notna())]

        df_14_12 = df[(df['count14'] >= count_min) & (df['count12'] >= count_min) & (df['14over12'].notna())]

    elif data_match_type == "Raw":
        df_10_12 = df[df['10over12'].notna()]

        df_13_12 = df[df['13over12'].notna()]

        df_14_12 = df[df['14over12'].notna()]

    else:
        raise Exception("Error: data match type " + data_match_type + " not recognized.")

    print("And now we are plotting the data...")

    # Compute binned velocity counts
    n_bins_vel = 48  # 25 m/s bins
    contour_range_vel = [axes[0][0].get_ylim(), axes[0][0].get_xlim()]
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
        pass

    # Compute binned height counts
    n_bins_h = 32  # 5 km bins
    contour_range_h = [axes[0][1].get_ylim(), axes[0][1].get_xlim()]
    try:
        binned_counts_10_12_h, bin_xedges_h, bin_yedges_h, bin_numbers_h = stats.binned_statistic_2d(
            df_10_12['height12'], df_10_12['height10'], values=None,
            statistic='count', bins=[n_bins_h, n_bins_h], range=contour_range_h)
    except:
        pass

    try:
        binned_counts_13_12_h, k, kk, kkk = stats.binned_statistic_2d(
            df_13_12['height12'], df_13_12['height13'], values=None,
            statistic='count', bins=[n_bins_h, n_bins_h], range=contour_range_h)
    except:
        pass

    try:
        binned_counts_14_12_h, l, ll, lll = stats.binned_statistic_2d(
            df_14_12['height12'], df_14_12['height14'], values=None,
            statistic='count', bins=[n_bins_h, n_bins_h], range=contour_range_h)
    except:
        pass

    # Compute bin centers, these will be the same for all frequency comparisons
    bin_xwidth_h = (bin_xedges_h[1] - bin_xedges_h[0])
    bin_ywidth_h = (bin_yedges_h[1] - bin_yedges_h[0])
    bin_xcenters_h = bin_xedges_h[1:] - bin_xwidth_h / 2
    bin_ycenters_h = bin_yedges_h[1:] - bin_ywidth_h / 2

    cmap = modified_jet(levels=n_levels - 1)

    # Plot 10 to 12 Comparison data in ROW: 0
    cont = axes[0][0].contourf(bin_xcenters_vel, bin_ycenters_vel, binned_counts_10_12_vel.transpose(),
                               levels=n_levels, cmap=cmap)
    fig.colorbar(cont, ax=axes[0][0])
    cont = axes[0][1].contourf(bin_xcenters_h, bin_ycenters_h, binned_counts_10_12_h.transpose(),
                               levels=n_levels, cmap=cmap)
    fig.colorbar(cont, ax=axes[0][1])

    axes[0][0].text(-490, 425, 'n=' + str(df_10_12.shape[0]), fontsize=12, c='k')
    axes[0][1].text(52, 177, 'n=' + str(df_10_12.shape[0]), fontsize=12, c='k')

    # Plot 13 to 12 Comparison data in ROW: 1
    cont = axes[1][0].contourf(bin_xcenters_vel, bin_ycenters_vel, binned_counts_13_12_vel.transpose(),
                               levels=n_levels, cmap=cmap)
    fig.colorbar(cont, ax=axes[1][0])
    cont = axes[1][1].contourf(bin_xcenters_h, bin_ycenters_h, binned_counts_13_12_h.transpose(),
                               levels=n_levels, cmap=cmap)
    fig.colorbar(cont, ax=axes[1][1])

    axes[1][0].text(-490, 425, 'n=' + str(df_13_12.shape[0]), fontsize=12, c='k')
    axes[1][1].text(52, 177, 'n=' + str(df_13_12.shape[0]), fontsize=12, c='k')

    # Plot 14 to 12 Comparison data in ROW: 2
    cont = axes[2][0].contourf(bin_xcenters_vel, bin_ycenters_vel, binned_counts_14_12_vel.transpose(),
                               levels=n_levels, cmap=cmap)
    fig.colorbar(cont, ax=axes[2][0])
    cont = axes[2][1].contourf(bin_xcenters_h, bin_ycenters_h, binned_counts_14_12_h.transpose(),
                               levels=n_levels, cmap=cmap)
    fig.colorbar(cont, ax=axes[2][1])

    axes[2][0].text(-490, 425, 'n=' + str(df_14_12.shape[0]), fontsize=12, c='k')
    axes[2][1].text(52, 177, 'n=' + str(df_14_12.shape[0]), fontsize=12, c='k')

    return fig


def format_subplots(axes):
    """
    :param axes: array of matplotlib.axis:
        The axes to format
    """

    bisector_colour = "purple"

    n_rows = axes.shape[0]
    n_cols = axes.shape[1]

    for row in range(n_rows):
        # First row will be velocity contour plots
        axes[row][0].set_xlabel("12 MHz Velocities [m/s]")
        axes[row][0].set_ylim([-600, 600])
        axes[row][0].set_xlim([-600, 600])
        axes[row][0].yaxis.set_minor_locator(MultipleLocator(100))
        axes[row][0].xaxis.set_minor_locator(MultipleLocator(100))
        axes[row][0].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
        axes[row][0].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
        axes[row][0].plot(axes[row][0].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
        axes[row][0].plot([0, 0], axes[row][0].get_xlim(), linestyle='-', linewidth=0.5, color='black')
        axes[row][0].plot([axes[row][0].get_ylim()[0], axes[row][0].get_ylim()[1]],
                          [axes[row][0].get_xlim()[0], axes[row][0].get_xlim()[1]],
                          linestyle='--', linewidth=0.75, color=bisector_colour)
        plt.sca(axes[row][0])
        plt.xticks([-600, -400, -200, 0, 200, 400, 600])
        axes[row][0].tick_params(axis='x', labelrotation=25)

        # Second row will be height comparison plots
        axes[row][1].set_xlabel("12 MHz Virtual Heights [km]")
        axes[row][1].set_ylim([40, 200])
        axes[row][1].set_xlim([40, 200])
        axes[row][1].yaxis.set_minor_locator(MultipleLocator(10))
        axes[row][1].xaxis.set_minor_locator(MultipleLocator(10))
        axes[row][1].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
        axes[row][1].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
        axes[row][1].plot(axes[row][1].get_ylim(), [120, 120], linestyle='-', linewidth=0.5, color='black')
        axes[row][1].plot([120, 120], axes[row][1].get_xlim(), linestyle='-', linewidth=0.5, color='black')
        axes[row][1].plot([axes[row][1].get_ylim()[0], axes[row][1].get_ylim()[1]],
                          [axes[row][1].get_xlim()[0], axes[row][1].get_xlim()[1]],
                          linestyle='--', linewidth=0.75, color=bisector_colour)

        axes[0][0].set_ylabel("10 MHz Velocities [m/s]")
        axes[1][0].set_ylabel("13 MHz Velocities [m/s]")
        axes[2][0].set_ylabel("14 MHz Velocities [m/s]")

        axes[0][1].set_ylabel("10 MHz Virtual Heights [km]")
        axes[1][1].set_ylabel("13 MHz Virtual Heights [km]")
        axes[2][1].set_ylabel("14 MHz Virtual Heights [km]")


if __name__ == '__main__':
    """

    """

    testing = True

    # station = "rkn"
    # area = None
    # year = "2016"  # yyyy
    # month = "09"  # mm
    # day = "26"  # dd
    # gate_range = (10, 40)
    # beam_range = (7, 7)
    # data_match_type = "Raw"  # "Median" or "Raw"
    # count_min = 4  # Only used for median matched data
    # start_hour = 0
    # end_hour = 4

    station = "rkn"
    area = 5
    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd
    gate_range = (0, 74)
    beam_range = (7, 7)
    data_match_type = "Median"  # "Median" or "Raw"
    count_min = 4  # Only used for median matched data
    start_hour = 0
    end_hour = 4

    if data_match_type == "Median":
        second_resolution = 60
    elif data_match_type == "Raw":
        second_resolution = 14
    else:
        raise Exception("Error: data match type " + data_match_type + " not recognized.")

    # Read in SuperDARN data
    loc_root = str((pathlib.Path().parent.absolute()).parent.absolute())
    in_dir = loc_root + "/RatioAnalysis/data/" + station + "/" + station + year + month + day

    if area is None:
        in_file = in_dir + "/" + station + year + month + day + "." + data_match_type + "MatchedData.1gg" \
                  + str(second_resolution) + "s.pbz2"
    else:
        in_file = in_dir + "/" + station + year + month + day + "." + data_match_type + "MatchedData.1gg" \
                  + str(second_resolution) + "s_area" + str(area) + ".pbz2"

    print("Reading in file: " + in_file)
    data_stream = bz2.BZ2File(in_file, "rb")
    df = cPickle.load(data_stream)

    fig = vel_and_height_contours(single_day_df=df, gate_range=gate_range, beam_range=beam_range,
                                  hour_range=(start_hour, end_hour), count_min=count_min, area=area)

    if testing:
        plt.show()
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        if area is None:
            out_fig = out_dir + "/" + "vel_and_height_contours-" + station + year + month + day \
                      + "-" + data_match_type + "MatchedData"
        else:
            out_fig = out_dir + "/" + "vel_and_height_contours-" + station + year + month + day \
                      + "-" + data_match_type + "MatchedData" + "_area" + str(area)

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
