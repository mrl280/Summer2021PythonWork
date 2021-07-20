import glob
import os
import pathlib
import warnings
import bz2

import _pickle as cPickle
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator
from PyPDF2 import PdfFileMerger


def add_data_to_plot(ax, time_restricted_df, freq_x, freq_y, data_match_type,
                     plot_similar_heights, plot_y_larger_than_x, plot_x_larger_than_y):
    """

    Add the scatter data to vel vs vel plot.

    :param ax: pyplot.axis: The axis to draw on
    :param time_restricted_df: pandas.DataFrame: The time restricted dataframe
    :param freq_x: str: The frequency to plot along x
    :param freq_y: str: The frequency to plot along y
    :param data_match_type: str: "Median" or "Raw"
    :param plot_x_larger_than_y: bool:
        If True, velocity scatter where x heights are larger than y heights is added to the plot
    :param plot_y_larger_than_x:
        If True, velocity scatter where y heights are larger than x heights is added to the plot
    :param plot_similar_heights:
        If True, velocity scatter where both frequencies see about the same heights is added to the plot
    """

    count_min = 4  # Only used for median matched data
    height_max = 180  # km  All points from heights above this are not considered

    similar_ratio_c = "black"
    ylarger_ratio_c = "blue"
    xlarger_ratio_c = "red"
    pts_size = 2

    txt_x_loc = -490  # For printing out point counts
    txt_y_loc = 425
    txt_fontsize = 12

    about_flag = "diffHeightFlag_" + freq_y + "about" + freq_x
    y_larger_flag = "diffHeightFlag_" + freq_y + "largerthan" + freq_x
    x_larger_flag = "diffHeightFlag_" + freq_y + "lessthan" + freq_x

    x_vel_str = "vel" + freq_x
    y_vel_str = "vel" + freq_y

    printed_out_count = False

    # Filter out bad ratio points
    if data_match_type == "Median":
        df = time_restricted_df[(time_restricted_df['count' + freq_y] >= count_min)
                                & (time_restricted_df['count' + freq_x] >= count_min)
                                & (time_restricted_df[freq_y + 'over' + freq_x].notna())]
    elif data_match_type == "Raw":
        df = time_restricted_df[time_restricted_df[freq_y + 'over' + freq_x].notna()]
    else:
        raise Exception("Error: data match type " + data_match_type + " not recognized.")

    # Remove points that are above the height maximum
    df = df.loc[(df['height' + freq_y] <= height_max) & (df['height' + freq_x] <= height_max)]
    df.reset_index(drop=True, inplace=True)

    # Go ahead and add the data to the plot
    if plot_similar_heights:
        ax.scatter(df.loc[df[about_flag], x_vel_str], df.loc[df[about_flag], y_vel_str],
                   s=pts_size, color=similar_ratio_c, marker='.', zorder=3, label="Similar Heights")
        if not plot_y_larger_than_x and not plot_x_larger_than_y:
            # We are just plotting similar heights data
            ax.text(txt_x_loc, txt_y_loc, 'n=' + str((df.loc[df[about_flag]]).shape[0]),
                    fontsize=txt_fontsize, color=similar_ratio_c)
            printed_out_count = True

    if plot_y_larger_than_x:
        ax.scatter(df.loc[df[y_larger_flag], x_vel_str], df.loc[df[y_larger_flag], y_vel_str],
                   s=pts_size, color=ylarger_ratio_c, marker='.', zorder=4, label="y height larger than x height")
        if not plot_similar_heights and not plot_x_larger_than_y:
            # We are just plotting y larger than x
            ax.text(txt_x_loc, txt_y_loc, 'n=' + str((df.loc[df[y_larger_flag]]).shape[0]),
                    fontsize=txt_fontsize, color=ylarger_ratio_c)
            printed_out_count = True

    if plot_x_larger_than_y:
        ax.scatter(df.loc[df[x_larger_flag], x_vel_str], df.loc[df[x_larger_flag], y_vel_str],
                   s=pts_size, color=xlarger_ratio_c, marker='.', zorder=4, label="x height larger than y height")
        if not plot_similar_heights and not plot_y_larger_than_x:
            # We are just plotting x larger than y
            ax.text(txt_x_loc, txt_y_loc, 'n=' + str((df.loc[df[x_larger_flag]]).shape[0]),
                    fontsize=txt_fontsize, color=xlarger_ratio_c)
            printed_out_count = True

    if plot_similar_heights and plot_y_larger_than_x and plot_x_larger_than_y:
        ax.text(txt_x_loc, txt_y_loc, 'n=' + str(df.shape[0]), fontsize=txt_fontsize, color='purple')
        printed_out_count = True

    if not printed_out_count:
        warnings.warn("Point count not supported for this combination of flags.", category=Warning)


def plottingVelRatios(station, year, month, day, gate_range, data_match_type, start_hour, end_hour,
                      plot_similar_heights, plot_y_larger_than_x, plot_x_larger_than_y,
                      testing=False):
    """

    Produces two types of plots:
        Velocity vs Velocity scatter
        Velocity ratio vs time

    :param station: str: The radar station as a three character string (e.g. 'rkn')
    :param year: str: The year to consider 'yyyy'
    :param month: str: The zero-padded month to consider 'mm'
    :param day: str: The zero-padded day to consider 'dd'
    :param gate_range: (int, int)
        Inclusive. The gate range to consider.  Gates start at 0.
    :param data_match_type: "Median" or "Raw" matched data
    :param start_hour: int: The starting hour in UT time
    :param end_hour: int: The ending hour in UT time
    :param plot_x_larger_than_y: bool:
        If True, velocity scatter where x heights are larger than y heights is added to the plot
    :param plot_y_larger_than_x: bool:
        If True, velocity scatter where y heights are larger than x heights is added to the plot
    :param plot_similar_heights: bool:
        If True, velocity scatter where both frequencies see about the same heights is added to the plot
    :param testing: bool: (optional - default: False)
        Set this to true if you are testing.
        If True, then instead of saving all time periods to file, the first time period will be shown
    """

    mnemonic = station.upper()
    gate_label = "gg: " + str(gate_range[0]) + "-" + str(gate_range[1])

    if data_match_type == "Median":
        second_resolution = 60
    elif data_match_type == "Raw":
        second_resolution = 14
    else:
        raise Exception("Error: data match type " + str(data_match_type) + " not recognized.")

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
    in_dir = loc_root + "/MultiFreqExperiment/RatioAnalysis/data/" + station + "/" + station + year + month + day
    in_file = in_dir + "/" + station + year + month + day + "." + \
              data_match_type + "MatchedData.1gg" + str(second_resolution) + "s.pbz2"
    data_stream = bz2.BZ2File(in_file, "rb")
    df = cPickle.load(data_stream)

    # Filter the data based on the expected gate range of the region of interest
    df = df.loc[(df['gate'] >= gate_range[0]) & (df['gate'] <= gate_range[1])]
    df.reset_index(drop=True, inplace=True)

    out_dir = loc_root + "/MultiFreqExperiment/RatioAnalysis/out/"

    if testing:
        time_chunks = 1
    else:
        time_chunks = (end_hour - start_hour) * 2  # Half hour periods

    for time_chunk in range(time_chunks):
        start_time = start_hour + time_chunk * 0.5
        end_time = (start_hour + 0.5) + time_chunk * 0.5

        print("     Plotting from: " + str(start_time) + " till " + str(end_time) + " UT")
        time_restricted_df = df[(df['decimalTime'] >= start_time) & (df['decimalTime'] < end_time)]

        # We will have 2 plots, one for vel v vel ratios and one for ratios v time
        # First lets just deal with the first plot

        """ Set up the first plot """
        n_rows = 3
        n_cols = 2
        fig1, ax1 = plt.subplots(figsize=(8, 9), dpi=300, nrows=n_rows, ncols=n_cols)
        plt.subplots_adjust(hspace=0.4, wspace=0.4)
        fig1.suptitle("Velocity Frequency Dependence: Vel vs Vel Scatter"
                      + "\n" + mnemonic + " " + year + "." + month + "." + day
                      + "; " + gate_label + "; " + data_match_type + " Matched Data"
                      + "; " + str(start_time) + "-" + str(end_time) + " UT; Max Height: 180 km"
                      + "\n Produced by " + str(os.path.basename(__file__)),
                      fontsize=13)

        # Format subplots
        for row in range(n_rows):
            ax1[row][0].set_xlabel("12 MHz Velocities [m/s]")
            ax1[row][1].set_xlabel("10 MHz Velocities [m/s]")

            for col in range(n_cols):
                ax1[row][col].set_ylim([-600, 600])
                ax1[row][col].set_xlim([-600, 600])
                ax1[row][col].yaxis.set_minor_locator(MultipleLocator(100))
                ax1[row][col].xaxis.set_minor_locator(MultipleLocator(100))
                ax1[row][col].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
                ax1[row][col].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
                ax1[row][col].plot(ax1[row][col].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
                ax1[row][col].plot([0, 0], ax1[row][col].get_xlim(), linestyle='-', linewidth=0.5, color='black')
                ax1[row][col].plot([ax1[row][col].get_ylim()[0], ax1[row][col].get_ylim()[1]],
                                   [ax1[row][col].get_xlim()[0], ax1[row][col].get_xlim()[1]],
                                   linestyle='-', linewidth=1, color='m', zorder=4)

        ax1[0][1].set_xlabel("13 MHz Velocities [m/s]")

        ax1[0][0].set_ylabel("10 MHz Velocities [m/s]")
        ax1[1][0].set_ylabel("13 MHz Velocities [m/s]")
        ax1[2][0].set_ylabel("14 MHz Velocities [m/s]")

        ax1[0][1].set_ylabel("14 MHz Velocities [m/s]")
        ax1[1][1].set_ylabel("13 MHz Velocities [m/s]")
        ax1[2][1].set_ylabel("14 MHz Velocities [m/s]")

        # Plot 10 to 12 Comparison data in ROW: 0, COL: 0
        add_data_to_plot(ax=ax1[0][0], time_restricted_df=time_restricted_df, freq_x='12', freq_y='10',
                         data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                         plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y)

        # Plot 13 to 12 Comparison data in ROW: 1, COL: 0
        add_data_to_plot(ax=ax1[1][0], time_restricted_df=time_restricted_df, freq_x='12', freq_y='13',
                         data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                         plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y)

        # Plot 14 to 12 Comparison data in ROW: 2, COL: 0
        add_data_to_plot(ax=ax1[2][0], time_restricted_df=time_restricted_df, freq_x='12', freq_y='14',
                         data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                         plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y)

        # Plot 14 to 13 Comparison data in ROW: 0, COL: 1
        add_data_to_plot(ax=ax1[0][1], time_restricted_df=time_restricted_df, freq_x='13', freq_y='14',
                         data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                         plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y)

        # Plot 13 to 10 Comparison data in ROW: 1, COL: 1
        add_data_to_plot(ax=ax1[1][1], time_restricted_df=time_restricted_df, freq_x='10', freq_y='13',
                         data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                         plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y)

        # Plot 14 to 10 Comparison data in ROW: 2, COL: 1
        add_data_to_plot(ax=ax1[2][1], time_restricted_df=time_restricted_df, freq_x='10', freq_y='14',
                         data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                         plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y)

        # Add legends to the plots
        ax1[0][0].legend(loc=(-0.2, 1.03), ncol=3)

        # """ Set up the second plot """
        # fig2, ax2 = plt.subplots(figsize=(8, 9), dpi=300, nrows=n_rows, ncols=n_cols)
        # plt.subplots_adjust(hspace=0.4, wspace=0.4)
        # fig2.suptitle("Velocity Frequency Dependence: Ratios vs Time"
        #               + "\n" + mnemonic + " " + year + "." + month + "." + day
        #               + " " + str(start_time) + "-" + str(end_time) + " UT"
        #               + "\nProduced by " + str(os.path.basename(__file__)),
        #               fontsize=13)
        #
        # # Format subplots in the second column
        # for row in range(n_rows):
        #     for col in range(n_cols):
        #         ax2[row][col].set_xlabel("Time, UT")
        #         ax2[row][col].set_ylim([0, 2])
        #         ax2[row][col].set_xlim([start_time, end_time])
        #         ax2[row][col].yaxis.set_minor_locator(MultipleLocator(1))
        #         ax2[row][col].grid(b=True, which='major', axis='y', linestyle='--', linewidth=0.5)
        #         ax2[row][col].grid(b=True, which='minor', axis='y', linestyle='--', linewidth=0.2)
        #         ax2[row][col].plot(ax2[row][col].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
        #         ax2[row][col].plot(ax2[row][col].get_ylim(), [1, 1], linestyle='-', linewidth=1, color='m')
        #
        # ax2[0][0].set_ylabel("Velocity Ratio: 10/12")
        # ax2[1][0].set_ylabel("Velocity Ratio: 13/12")
        # ax2[2][0].set_ylabel("Velocity Ratio: 14/12")
        #
        # ax2[0][1].set_ylabel("Velocity Ratio: 14/13")
        # ax2[1][1].set_ylabel("Velocity Ratio: 13/10")
        # ax2[2][1].set_ylabel("Velocity Ratio: 14/10")
        #
        # # Plot 10 to 12 Comparison data in ROW: 0, COL: 0
        # ax2[0][0].scatter(df_10_12['decimalTime'], df_10_12['10over12'], s=4, color='k', marker='.')
        #
        # # Plot 13 to 12 Comparison data in ROW: 1, COL: 0
        # ax2[1][0].scatter(df_13_12['decimalTime'], df_13_12['13over12'], s=4, color='k', marker='.')
        #
        # # Plot 14 to 12 Comparison data in ROW: 2, COL: 0
        # ax2[2][0].scatter(df_14_12['decimalTime'], df_14_12['14over12'], s=4, color='k', marker='.')
        #
        # # Plot 12 to 10 Comparison data in ROW: 0, COL: 1
        # ax2[0][1].scatter(df_14_13['decimalTime'], df_14_13['14over13'], s=4, color='k', marker='.')
        #
        # # Plot 13 to 10 Comparison data in ROW: 1, COL: 1
        # ax2[1][1].scatter(df_13_10['decimalTime'], df_13_10['13over10'], s=4, color='k', marker='.')
        #
        # # Plot 14 to 10 Comparison data in ROW: 2, COL: 1
        # ax2[2][1].scatter(df_14_10['decimalTime'], df_14_10['14over10'], s=4, color='k', marker='.')

        if testing:
            plt.show()
        else:
            # Save the files as a temp file
            fig1.savefig(out_dir + "/" + mnemonic + "_vel_comp_" + year + month + day
                         + "_" + chr(ord('a') + time_chunk) + " " + str(start_time) + "-" + str(end_time) + "UT"
                         + "_temp1.pdf", format='pdf', dpi=300)
            # fig2.savefig(out_dir + "/" + mnemonic + "_vel_comp_" + year + month + day
            #              + "_" + chr(ord('a') + time_chunk) + " " + str(start_time) + "-" + str(end_time) + "UT"
            #              + "_temp2.pdf", format='pdf', dpi=300)
            plt.close(fig1)

    if plot_similar_heights and not plot_y_larger_than_x and not plot_x_larger_than_y:
        color_str = "_similarHeightsOnly"
    elif not plot_similar_heights and plot_y_larger_than_x and not plot_x_larger_than_y:
        color_str = "_YLargerThanXOnly"
    elif not plot_similar_heights and not plot_y_larger_than_x and plot_x_larger_than_y:
        color_str = "_XLargerThanYOnly"
    else:
        warnings.warn("No colour string added to outfile", category=Warning)
        color_str = ""

    if not testing:
        # Merge all the temp pdf files
        merger = PdfFileMerger()
        for pdf in glob.iglob("out/*_temp1.pdf"):
            merger.append(pdf)
        for pdf in glob.iglob("out/*_temp2.pdf"):
            merger.append(pdf)
        with open(out_dir + "/" + mnemonic + "_vel_comp_" + year + month + day
                  + "_gg" + str(gate_range[0]) + "-" + str(gate_range[1]) + "_" + data_match_type + color_str + ".pdf",
                  "wb") as fout:
            merger.write(fout)
        merger.close()

        for file in glob.iglob("out/*_temp*.pdf"):
            os.remove(file)


if __name__ == '__main__':
    """ Testing and running """

    testing = False

    year = "2017"  # yyyy
    month = "10"  # mm
    day = "23"  # dd

    station = "rkn"
    gate_range = [10, 30]
    data_match_type = "Raw"  # "Median" or "Raw"

    start_hour = 4
    end_hour = 8

    print("Running for just similar heights - just black dots")
    plottingVelRatios(station=station, year=year, month=month, day=day, gate_range=gate_range,
                      data_match_type=data_match_type, start_hour=start_hour, end_hour=end_hour,
                      plot_similar_heights=True, plot_y_larger_than_x=False,
                      plot_x_larger_than_y=False,
                      testing=testing)

    print("Running for just y larger than x - just blue dots")
    plottingVelRatios(station=station, year=year, month=month, day=day, gate_range=gate_range,
                      data_match_type=data_match_type, start_hour=start_hour, end_hour=end_hour,
                      plot_similar_heights=False, plot_y_larger_than_x=True,
                      plot_x_larger_than_y=False,
                      testing=testing)

    print("Running for just x larger than y - just red dots")
    plottingVelRatios(station=station, year=year, month=month, day=day, gate_range=gate_range,
                      data_match_type=data_match_type, start_hour=start_hour, end_hour=end_hour,
                      plot_similar_heights=False, plot_y_larger_than_x=False,
                      plot_x_larger_than_y=True,
                      testing=testing)

    print("Running for all points - all colours of dots")
    plottingVelRatios(station=station, year=year, month=month, day=day, gate_range=gate_range,
                      data_match_type=data_match_type, start_hour=start_hour, end_hour=end_hour,
                      plot_similar_heights=True, plot_y_larger_than_x=True,
                      plot_x_larger_than_y=True,
                      testing=testing)
