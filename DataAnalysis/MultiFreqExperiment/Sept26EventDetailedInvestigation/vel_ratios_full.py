import os
import pathlib
import bz2
import warnings

import matplotlib.pyplot as plt
import _pickle as cPickle

from matplotlib.ticker import MultipleLocator

from DataAnalysis.MultiFreqExperiment.RatioAnalysis.plottingVelRatios import add_data_to_plot


def height_ratios_scatter(single_day_df, beam_range, gate_range, hour_range,
                          plot_similar_heights, plot_y_larger_than_x, plot_x_larger_than_y,
                          count_min=4, area=None):
    """

    Create scatter plots comparing velocities at one frequency to velocities at another frequency
    This program runs for all the data in an event, for half hour intervals use plottingVelRatios.py

    Note: If point count is greater than points you see it is because the rest are out of frame.

    :param single_day_df: pandas.Data_Frame:
            A single days worth of data that we want to range profile
    :param gate_range: (int, int):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (int, int):
             Just for the figure title.
    :param hour_range: (int, int):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not quite inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param plot_x_larger_than_y: bool:
            If True, velocity scatter where x heights are larger than y heights is added to the plot
    :param plot_y_larger_than_x: bool:
            If True, velocity scatter where y heights are larger than x heights is added to the plot
    :param plot_similar_heights: bool:
          If True, velocity scatter where both frequencies see about the same heights is added to the plot

    :param count_min: int (optional):
            The minimum number of points required to count a median match.
            This parameter has no effect on raw matched data
    :param area: int:
            The numbered area of interest.

    :return: matplotlib.pyplot.figure:
            The figure, it can then be viewed, modified, or saved to file
    """

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
    title_font_size = 14

    fig, axes = plt.subplots(figsize=[8, 9], nrows=3, ncols=2, constrained_layout=True, dpi=300)
    format_subplots(axes)

    if area is None:
        fig.suptitle(date_string + " at " + station.upper() + "; " + beam_string + "; " + gate_string
                     + "; " + data_match_type + " Matched Data"
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)),
                     fontsize=title_font_size)
    else:
        fig.suptitle(date_string + " at " + station.upper() + "; Area " + str(area)
                     + "; " + data_match_type + " Matched Data"
                     + "\n" + resolution_string + "; " + "Produced by " + str(os.path.basename(__file__)),
                     fontsize=title_font_size)

    print("Adding data to plots...")

    # Plot 10 to 12 Comparison data in ROW: 0, COL: 0
    add_data_to_plot(ax=axes[0][0], time_restricted_df=df, freq_x='12', freq_y='10',
                     data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                     plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y,
                     count_min=count_min)

    # Plot 13 to 12 Comparison data in ROW: 1, COL: 0
    add_data_to_plot(ax=axes[1][0], time_restricted_df=df, freq_x='12', freq_y='13',
                     data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                     plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y,
                     count_min=count_min)

    # Plot 14 to 12 Comparison data in ROW: 2, COL: 0
    add_data_to_plot(ax=axes[2][0], time_restricted_df=df, freq_x='12', freq_y='14',
                     data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                     plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y,
                     count_min=count_min)

    # Plot 14 to 13 Comparison data in ROW: 0, COL: 1
    add_data_to_plot(ax=axes[0][1], time_restricted_df=df, freq_x='13', freq_y='14',
                     data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                     plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y,
                     count_min=count_min)

    # Plot 13 to 10 Comparison data in ROW: 1, COL: 1
    add_data_to_plot(ax=axes[1][1], time_restricted_df=df, freq_x='10', freq_y='13',
                     data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                     plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y,
                     count_min=count_min)

    # Plot 14 to 10 Comparison data in ROW: 2, COL: 1
    add_data_to_plot(ax=axes[2][1], time_restricted_df=df, freq_x='10', freq_y='14',
                     data_match_type=data_match_type, plot_similar_heights=plot_similar_heights,
                     plot_y_larger_than_x=plot_y_larger_than_x, plot_x_larger_than_y=plot_x_larger_than_y,
                     count_min=count_min)

    axes[0][0].legend(loc=(0.1, 1.03))

    return fig


def format_subplots(axes):
    """
    :param axes: array of matplotlib.axis:
        The axes to format
    """

    n_rows = axes.shape[0]
    n_cols = axes.shape[1]

    # Format subplots
    for row in range(n_rows):
        axes[row][0].set_xlabel("12 MHz Velocities [m/s]")
        axes[row][1].set_xlabel("10 MHz Velocities [m/s]")

        for col in range(n_cols):
            axes[row][col].set_ylim([-600, 600])
            axes[row][col].set_xlim([-600, 600])
            axes[row][col].yaxis.set_minor_locator(MultipleLocator(100))
            axes[row][col].xaxis.set_minor_locator(MultipleLocator(100))
            axes[row][col].grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5)
            axes[row][col].grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2)
            axes[row][col].plot(axes[row][col].get_ylim(), [0, 0], linestyle='-', linewidth=0.5, color='black')
            axes[row][col].plot([0, 0], axes[row][col].get_xlim(), linestyle='-', linewidth=0.5, color='black')
            axes[row][col].plot([axes[row][col].get_ylim()[0], axes[row][col].get_ylim()[1]],
                                [axes[row][col].get_xlim()[0], axes[row][col].get_xlim()[1]],
                                linestyle='-', linewidth=1, color='m', zorder=4)

    axes[0][1].set_xlabel("13 MHz Velocities [m/s]")

    axes[0][0].set_ylabel("10 MHz Velocities [m/s]")
    axes[1][0].set_ylabel("13 MHz Velocities [m/s]")
    axes[2][0].set_ylabel("14 MHz Velocities [m/s]")

    axes[0][1].set_ylabel("14 MHz Velocities [m/s]")
    axes[1][1].set_ylabel("13 MHz Velocities [m/s]")
    axes[2][1].set_ylabel("14 MHz Velocities [m/s]")


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
    # data_match_type = "Median"  # "Median" or "Raw"
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

    loc_root = str((pathlib.Path().parent.absolute()))
    out_dir = loc_root + "/out"
    if area is None:
        out_fig_base = out_dir + "/" + "vel_ratios_full-" + station + year + month + day \
                  + "-" + data_match_type + "MatchedData"
    else:
        out_fig_base = out_dir + "/" + "vel_ratios_full-" + station + year + month + day \
                  + "-" + data_match_type + "MatchedData" + "_area" + str(area)

    """ RUN SIMILAR HEIGHTS """
    fig_similar_heights = height_ratios_scatter(single_day_df=df, gate_range=gate_range, beam_range=beam_range,
                                                hour_range=(start_hour, end_hour),
                                                plot_similar_heights=True,
                                                plot_y_larger_than_x=False,
                                                plot_x_larger_than_y=False,
                                                count_min=count_min, area=area)
    if testing:
        plt.show()
    else:
        out_fig = out_fig_base + "_similar_heights.jpg"
        print("Saving plot as " + out_fig)
        fig_similar_heights.savefig(out_fig, format='jpg', dpi=300)

    """ RUN Y LARGER """
    fig_y_larger_than_x = height_ratios_scatter(single_day_df=df, gate_range=gate_range, beam_range=beam_range,
                                                hour_range=(start_hour, end_hour),
                                                plot_similar_heights=False,
                                                plot_y_larger_than_x=True,
                                                plot_x_larger_than_y=False,
                                                count_min=count_min, area=area)
    if testing:
        plt.show()
    else:
        out_fig = out_fig_base + "_y_larger_than_x.jpg"
        print("Saving plot as " + out_fig)
        fig_y_larger_than_x.savefig(out_fig, format='jpg', dpi=300)

    """ RUN X LARGER """
    fig_x_larger_than_y = height_ratios_scatter(single_day_df=df, gate_range=gate_range, beam_range=beam_range,
                                                hour_range=(start_hour, end_hour),
                                                plot_similar_heights=False,
                                                plot_y_larger_than_x=False,
                                                plot_x_larger_than_y=True,
                                                count_min=count_min, area=area)
    if testing:
        plt.show()
    else:
        out_fig = out_fig_base + "_x_larger_than_y.jpg"
        print("Saving plot as " + out_fig)
        fig_x_larger_than_y.savefig(out_fig, format='jpg', dpi=300)

    """ RUN FOR ALL POINTS """
    fig_all_points = height_ratios_scatter(single_day_df=df, gate_range=gate_range, beam_range=beam_range,
                                           hour_range=(start_hour, end_hour),
                                           plot_similar_heights=True,
                                           plot_y_larger_than_x=True,
                                           plot_x_larger_than_y=True,
                                           count_min=count_min, area=area)
    if testing:
        plt.show()
    else:
        out_fig = out_fig_base + "_all_points.jpg"
        print("Saving plot as " + out_fig)
        fig_all_points.savefig(out_fig, format='jpg', dpi=300)
