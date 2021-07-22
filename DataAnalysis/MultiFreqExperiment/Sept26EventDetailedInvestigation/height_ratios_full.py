import os
import pathlib
import bz2

import _pickle as cPickle
import warnings

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import MultipleLocator


def height_ratios_scatter(single_day_df, beam_range, gate_range, hour_range, count_min=4, area=None):
    """

    Create scatter plots comparing virtual height at one frequency to virtual height at another frequency
    This program runs for all the data in an event, for half hour intervals use plottingHeightRatios.py

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

    normal_ratio_c = "black"
    ylarger_ratio_c = "blue"
    xlarger_ratio_c = "red"
    pts_size = 2

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

    date_string = "September 26, 2016"

    fig, axes = plt.subplots(figsize=[8, 9], nrows=3, ncols=2, constrained_layout=True, dpi=300)
    format_subplots(axes)

    if area is None:
        fig.suptitle(date_string + " at " + station.upper() + "; " + beam_string + "; " + gate_string
                     + "; " + data_match_type + " Matched Data"
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)
    else:
        fig.suptitle(date_string + " at " + station.upper() + "; Area " + str(area)
                     + "; " + data_match_type + " Matched Data"
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)), fontsize=18)

    print("Building frequency dependent dataframes...")

    # Build frequency dependant data frames
    if data_match_type == "Median":
        df_10_12 = df[(df['count10'] >= count_min) & (df['count12'] >= count_min) & (df['10over12'].notna())]

        df_13_12 = df[(df['count13'] >= count_min) & (df['count12'] >= count_min) & (df['13over12'].notna())]

        df_14_12 = df[(df['count14'] >= count_min) & (df['count12'] >= count_min) & (df['14over12'].notna())]

        df_14_13 = df[(df['count14'] >= count_min) & (df['count13'] >= count_min) & (df['14over13'].notna())]

        df_13_10 = df[(df['count13'] >= count_min) & (df['count10'] >= count_min) & (df['13over10'].notna())]

        df_14_10 = df[(df['count14'] >= count_min) & (df['count10'] >= count_min) & (df['14over10'].notna())]

    elif data_match_type == "Raw":
        df_10_12 = df[df['10over12'].notna()]

        df_13_12 = df[df['13over12'].notna()]

        df_14_12 = df[df['14over12'].notna()]

        df_14_13 = df[df['14over13'].notna()]

        df_13_10 = df[df['13over10'].notna()]

        df_14_10 = df[df['14over10'].notna()]

    else:
        raise Exception("Error: data match type " + data_match_type + " not recognized.")

    print("And now we are plotting the data...")

    # Plot 10 to 12 Comparison data in ROW: 0, COL: 0
    axes[0][0].scatter(df_10_12.loc[df_10_12['diffHeightFlag_10largerthan12'], 'height12'],
                       df_10_12.loc[df_10_12['diffHeightFlag_10largerthan12'], 'height10'],
                       s=pts_size, color=ylarger_ratio_c, marker='.', zorder=3, label="y > x")
    axes[0][0].scatter(df_10_12.loc[df_10_12['diffHeightFlag_10lessthan12'], 'height12'],
                       df_10_12.loc[df_10_12['diffHeightFlag_10lessthan12'], 'height10'],
                       s=pts_size, color=xlarger_ratio_c, marker='.', zorder=3, label="x > y")
    axes[0][0].scatter(df_10_12.loc[df_10_12['diffHeightFlag_10about12'], 'height12'],
                       df_10_12.loc[df_10_12['diffHeightFlag_10about12'], 'height10'],
                       s=pts_size, color=normal_ratio_c, marker='.', zorder=3,
                       label="x ~ y")
    axes[0][0].text(52, 172, 'n=' + str(df_10_12.shape[0]), fontsize=12, c='purple')

    # Plot 13 to 12 Comparison data in ROW: 1, COL: 0
    axes[1][0].scatter(df_13_12.loc[df_13_12['diffHeightFlag_13largerthan12'], 'height12'],
                       df_13_12.loc[df_13_12['diffHeightFlag_13largerthan12'], 'height13'],
                       s=4, color=ylarger_ratio_c, marker='.', zorder=3)
    axes[1][0].scatter(df_13_12.loc[df_13_12['diffHeightFlag_13lessthan12'], 'height12'],
                       df_13_12.loc[df_13_12['diffHeightFlag_13lessthan12'], 'height13'],
                       s=4, color=xlarger_ratio_c, marker='.', zorder=3)
    axes[1][0].scatter(df_13_12.loc[df_13_12['diffHeightFlag_13about12'], 'height12'],
                       df_13_12.loc[df_13_12['diffHeightFlag_13about12'], 'height13'],
                       s=4, color=normal_ratio_c, marker='.', zorder=3)
    axes[1][0].text(52, 172, 'n=' + str(df_13_12.shape[0]), fontsize=12, c='purple')

    # Plot 14 to 12 Comparison data in ROW: 2, COL: 0
    axes[2][0].scatter(df_14_12.loc[df_14_12['diffHeightFlag_14largerthan12'], 'height12'],
                       df_14_12.loc[df_14_12['diffHeightFlag_14largerthan12'], 'height14'],
                       s=4, color=ylarger_ratio_c, marker='.', zorder=3)
    axes[2][0].scatter(df_14_12.loc[df_14_12['diffHeightFlag_14lessthan12'], 'height12'],
                       df_14_12.loc[df_14_12['diffHeightFlag_14lessthan12'], 'height14'],
                       s=4, color=xlarger_ratio_c, marker='.', zorder=3)
    axes[2][0].scatter(df_14_12.loc[df_14_12['diffHeightFlag_14about12'], 'height12'],
                       df_14_12.loc[df_14_12['diffHeightFlag_14about12'], 'height14'],
                       s=4, color=normal_ratio_c, marker='.', zorder=3)
    axes[2][0].text(52, 172, 'n=' + str(df_14_12.shape[0]), fontsize=12, c='purple')

    # Plot 14 to 13 Comparison data in ROW: 0, COL: 1
    axes[0][1].scatter(df_14_13.loc[df_14_13['diffHeightFlag_14largerthan13'], 'height13'],
                       df_14_13.loc[df_14_13['diffHeightFlag_14largerthan13'], 'height14'],
                       s=4, color=ylarger_ratio_c, marker='.', zorder=3)
    axes[0][1].scatter(df_14_13.loc[df_14_13['diffHeightFlag_14lessthan13'], 'height13'],
                       df_14_13.loc[df_14_13['diffHeightFlag_14lessthan13'], 'height14'],
                       s=4, color=xlarger_ratio_c, marker='.', zorder=3)
    axes[0][1].scatter(df_14_13.loc[df_14_13['diffHeightFlag_14about13'], 'height13'],
                       df_14_13.loc[df_14_13['diffHeightFlag_14about13'], 'height14'],
                       s=4, color=normal_ratio_c, marker='.', zorder=3)
    axes[0][1].text(52, 172, 'n=' + str(df_14_13.shape[0]), fontsize=12, c='purple')

    # Plot 13 to 10 Comparison data in ROW: 1, COL: 1
    axes[1][1].scatter(df_13_10.loc[df_13_10['diffHeightFlag_13largerthan10'], 'height10'],
                       df_13_10.loc[df_13_10['diffHeightFlag_13largerthan10'], 'height13'],
                       s=4, color=ylarger_ratio_c, marker='.', zorder=3)
    axes[1][1].scatter(df_13_10.loc[df_13_10['diffHeightFlag_13lessthan10'], 'height10'],
                       df_13_10.loc[df_13_10['diffHeightFlag_13lessthan10'], 'height13'],
                       s=4, color=xlarger_ratio_c, marker='.', zorder=3)
    axes[1][1].scatter(df_13_10.loc[df_13_10['diffHeightFlag_13about10'], 'height10'],
                       df_13_10.loc[df_13_10['diffHeightFlag_13about10'], 'height13'],
                       s=4, color=normal_ratio_c, marker='.', zorder=3)
    axes[1][1].text(52, 172, 'n=' + str(df_13_10.shape[0]), fontsize=12, c='purple')

    # Plot 14 to 10 Comparison data in ROW: 2, COL: 1
    axes[2][1].scatter(df_14_10.loc[df_14_10['diffHeightFlag_14largerthan10'], 'height10'],
                       df_14_10.loc[df_14_10['diffHeightFlag_14largerthan10'], 'height14'],
                       s=4, color=ylarger_ratio_c, marker='.', zorder=3)
    axes[2][1].scatter(df_14_10.loc[df_14_10['diffHeightFlag_14lessthan10'], 'height10'],
                       df_14_10.loc[df_14_10['diffHeightFlag_14lessthan10'], 'height14'],
                       s=4, color=xlarger_ratio_c, marker='.', zorder=3)
    axes[2][1].scatter(df_14_10.loc[df_14_10['diffHeightFlag_14about10'], 'height10'],
                       df_14_10.loc[df_14_10['diffHeightFlag_14about10'], 'height14'],
                       s=4, color=normal_ratio_c, marker='.', zorder=3)
    axes[2][1].text(52, 172, 'n=' + str(df_14_10.shape[0]), fontsize=12, c='purple')

    # Add legends to the plots
    # axes[0][0].legend(loc=(0.2, 1.03), ncol=3)
    axes[0][0].legend(loc='lower left')

    # df = df.loc[df['height10'].notna()]
    # print(np.min(df['height10'].unique()))

    return fig


def format_subplots(axes):
    """
    :param axes: array of matplotlib.axis:
        The axes to format
    """

    n_rows = axes.shape[0]
    n_cols = axes.shape[1]

    for row in range(n_rows):
        axes[row][0].set_xlabel("12 MHz Virtual Heights [km]")
        axes[row][1].set_xlabel("10 MHz Virtual Heights [km]")

        for col in range(n_cols):
            axes[row][col].set_ylim([40, 200])
            axes[row][col].set_xlim([40, 200])
            axes[row][col].yaxis.set_minor_locator(MultipleLocator(10))
            axes[row][col].xaxis.set_minor_locator(MultipleLocator(10))
            axes[row][col].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
            axes[row][col].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
            axes[row][col].plot(axes[row][col].get_ylim(), [120, 120], linestyle='-', linewidth=0.5, color='black')
            axes[row][col].plot([120, 120], axes[row][col].get_xlim(), linestyle='-', linewidth=0.5, color='black')
            axes[row][col].plot([axes[row][col].get_ylim()[0], axes[row][col].get_ylim()[1]],
                                [axes[row][col].get_xlim()[0], axes[row][col].get_xlim()[1]],
                                linestyle='-', linewidth=1, color='m', zorder=5)
            axes[row][col].plot([axes[row][col].get_ylim()[0] / 0.92, axes[row][col].get_ylim()[1] / 0.92],
                                [axes[row][col].get_xlim()[0], axes[row][col].get_xlim()[1]],
                                linestyle='--', linewidth=0.5, color='m', zorder=5)
            axes[row][col].plot([axes[row][col].get_ylim()[0] / 1.08, axes[row][col].get_ylim()[1] / 1.08],
                                [axes[row][col].get_xlim()[0], axes[row][col].get_xlim()[1]],
                                linestyle='--', linewidth=0.5, color='m', zorder=5)

    axes[0][1].set_xlabel("13 MHz Virtual Heights [km]")

    axes[0][0].set_ylabel("10 MHz Virtual Heights [km]")
    axes[1][0].set_ylabel("13 MHz Virtual Heights [km]")
    axes[2][0].set_ylabel("14 MHz Virtual Heights [km]")

    axes[0][1].set_ylabel("14 MHz Virtual Heights [km]")
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
    area = "3c"  # options: [1, 2, 3, 4, 5, "3a", "3b", "3c"]
    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd
    gate_range = (0, 74)
    beam_range = (7, 7)
    data_match_type = "Raw"  # "Median" or "Raw"
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

    fig = height_ratios_scatter(single_day_df=df, gate_range=gate_range, beam_range=beam_range,
                                hour_range=(start_hour, end_hour), count_min=count_min, area=area)

    if testing:
        plt.show()
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        if area is None:
            out_fig = out_dir + "/" + "height_ratios_full-" + station + year + month + day \
                      + "-" + data_match_type + "MatchedData"
        else:
            out_fig = out_dir + "/" + "height_ratios_full-" + station + year + month + day \
                      + "-" + data_match_type + "MatchedData" + "_area" + str(area)

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
