import calendar
import math
import os
import pathlib
import pydarn

import numpy as np

from bisect import bisect
from matplotlib.ticker import MultipleLocator
from pydarn import SuperDARNRadars
from matplotlib import pyplot as plt

from lib.data_getters.get_IMF_data import get_IMF_data
from lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def occ_imf_variation(station, year=None, day_range=None, hour_range=None,
                      gate_range=None, beam_range=None, freq_range=None,
                      time_units='ut', plot_type='contour', local_testing=False, even_odd_days=None):
    """

    Plot the y-component of the interplanetary magnetic field (By) along x, and the z-component (Bz) along y.
      Occurrence is then plotted as contours or pixels.  One plot for every month.

    Idea is to show that IMF drives echo occurrence rates.

    Uses OMNI IMF data.  IMF data should be pickled with pickle_imf() before running this program.
     An error is raised if no pickled data is found.

    More on OMNI data here: https://omniweb.gsfc.nasa.gov/html/omni_min_data.html
    OMNI data files were created and downloaded from here: https://omniweb.gsfc.nasa.gov/form/omni_min.html

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - To check which fitACF program is being used, refer to the data readers in lib.data_getters

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year: int:
            The year to consider.
    :param day_range: (int, int) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param hour_range: (int, int) (optional):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not quite inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
            Note: Choosing a limited hour range slows things down because time has to be computed in time_units
             just so we can perform hour restriction.
    :param gate_range: (int, int) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (int, int) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param freq_range: (float, float) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.

    :param time_units: str (optional; default is 'ut')
            time_units is only used if hour_range is not (0, 24)
                'ut' for universal time
                'mlt' for magnetic local time
                'lt' for local time (based on longitude)
                'lst' for local standard time (based on time zones)
                'ast' for apparent solar time (based on the apparent angular motion of the sun across the sky)
    :param even_odd_days: (optional; default is None)
            'even': only even days are read in
            'odd': only odd days are read in
            None: all days are read in
    :param local_testing: bool (optional; default is False)
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :param plot_type: str (optional; default is 'contour'):
            The type of plot, either 'contour' or 'pixel'

    :return fig: matplotlib.pyplot.figure
            The figure created, it can then be modified, added to, printed out, or saved in whichever file format is
            desired.
    """

    vmaxs = {"is": 0.6,  # The top of the colour bar
             "gs": 0.2}
    title_fontsize = 18
    cbar_label_size = 10
    cbar_text_format = '%.2f'
    x_lim = (-10, 10)
    y_lim = (-10, 10)

    year = check_year(year=year)
    freq_range = check_freq_range(freq_range=freq_range)
    day_range = check_day_range(day_range=day_range)
    hour_range = check_hour_range(hour_range=hour_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid

    gate_range = check_gate_range(gate_range=gate_range, hdw_info=this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range=beam_range, hdw_info=this_radars_info.hardware_info)

    print("     Retrieving IMF data...")
    imf_df = get_IMF_data(year_range=(year, year))

    print("     Retrieving SuperDARN data...")
    df = get_data_handler(station, year_range=(year, year), month_range=None, day_range=day_range,
                          gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                          occ_data=True, local_testing=local_testing, even_odd_days=even_odd_days)
    df = only_keep_45km_res_data(df)

    # We will need ordered date-times for bisect matching, so just to be safe..
    df.sort_values(by='datetime', ascending=True, inplace=True)

    # Perform hour restriction if required
    if hour_range != (0, 24):
        # Add decimal hour to df in whatever units were requested
        # Use the middle of the year as magnetic field estimate
        date_time_est, _ = build_datetime_epoch(year=year, month=6, day=15, hour=0)
        df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id, date_time_est=date_time_est)
        df = df.loc[(df[time_units] >= hour_range[0]) & (df[time_units] <= hour_range[1])]
        df.reset_index(drop=True, inplace=True)

    print("     Preparing the figure...")
    fig = plt.figure(figsize=[10, 14], constrained_layout=True, dpi=300)
    month_axes, cbar_axes = add_axes(fig=fig)
    apply_subplot_formatting(axes=month_axes, x_lim=x_lim, y_lim=y_lim)

    title_figure(fig=fig, station=station, year=year, time_units=time_units, beam_range=beam_range,
                 gate_range=gate_range, freq_range=freq_range, hour_range=hour_range, day_range=day_range,
                 title_fontsize=title_fontsize)

    print("     Adding month to df...")
    month = []
    for i in range(len(df)):
        month.append(df['datetime'].iat[i].month)
    df['month'] = month

    # Compute By edges (x-axis edges)
    n_bins_By = 20
    By_edges = np.linspace(x_lim[0], x_lim[1], num=(n_bins_By + 1))

    # Compute Bz edges (y-axis edges)
    n_bins_Bz = 20
    Bz_edges = np.linspace(y_lim[0], y_lim[1], num=(n_bins_Bz + 1))

    for month_str, month_axis in month_axes.items():

        month = list(calendar.month_abbr).index(month_str[:3])  # Pull the month number from the month string

        # Month restrict both dataframes, this should speed up the date matching required for IMF assignment
        df_mm = df[df['month'] == month].copy()
        df_mm.reset_index(drop=True, inplace=True)
        imf_df_mm = imf_df[imf_df['month'] == month].copy()
        imf_df_mm.reset_index(drop=True, inplace=True)

        if len(df_mm) <= 0 or len(imf_df_mm) <= 0:
            # Then there is nothing we can do here
            continue

        print("     Assigning IMF values for " + month_str)
        df_mm = assign_imf_values(df=df_mm, imf_df=imf_df_mm)

        print("     Computing binned occ rates...")
        contour_data_is, contour_data_gs = bin_data(df=df_mm, By_edges=By_edges, Bz_edges=Bz_edges)

        print("     Plotting the data...")
        if plot_type == "contour":

            delta_By = By_edges[1] - By_edges[0]
            delta_Bz = Bz_edges[1] - Bz_edges[0]

            # Compute bin centers
            bin_xcenters = By_edges[1:] - delta_By / 2
            bin_ycenters = Bz_edges[1:] - delta_Bz / 2

            n_levels = 12
            is_levels = np.linspace(start=0, stop=vmaxs['is'], num=(n_levels + 1))
            gs_levels = np.linspace(start=0, stop=vmaxs['gs'], num=(n_levels + 1))
            cmap = 'jet'

            # cmap = modified_jet(levels=len(levels) - 1)

            if vmaxs['is'] < 1:
                is_plot = month_axis['is'].contourf(bin_xcenters, bin_ycenters, contour_data_is.transpose(),
                                                    cmap=cmap, levels=is_levels, zorder=0, extend='max')
            else:
                is_plot = month_axis['is'].contourf(bin_xcenters, bin_ycenters, contour_data_is.transpose(),
                                                    cmap=cmap, levels=is_levels, zorder=0)

            if vmaxs['gs'] < 1:
                gs_plot = month_axis['gs'].contourf(bin_xcenters, bin_ycenters, contour_data_gs.transpose(),
                                                    cmap=cmap, levels=gs_levels, zorder=0, extend='max')
            else:
                gs_plot = month_axis['gs'].contourf(bin_xcenters, bin_ycenters, contour_data_gs.transpose(),
                                                    cmap=cmap, levels=gs_levels, zorder=0)

        elif plot_type == "pixel":

            cmap = 'jet'
            is_plot = month_axis['is'].pcolormesh(By_edges, Bz_edges, contour_data_is.transpose(),
                                                  cmap=cmap, vmin=0, vmax=vmaxs['is'], zorder=0)
            gs_plot = month_axis['gs'].pcolormesh(By_edges, Bz_edges, contour_data_gs.transpose(),
                                                  cmap=cmap, vmin=0, vmax=vmaxs['gs'], zorder=0)

        else:
            raise Exception("plot_type not recognized")

        # Add in colour bars for each plot
        # The extend here is needed for the pixel option, it has no effect on contours,
        #  contour extend was set when plotting
        if vmaxs['is'] < 1:
            is_cbar = fig.colorbar(is_plot, cax=cbar_axes['is'], orientation='horizontal', format=cbar_text_format,
                                   extend='max')
        else:
            is_cbar = fig.colorbar(is_plot, cax=cbar_axes['is'], orientation='horizontal', format=cbar_text_format)
        is_cbar.ax.tick_params(labelsize=cbar_label_size)

        if vmaxs['gs'] < 1:
            cbar1 = fig.colorbar(gs_plot, cax=cbar_axes['gs'], orientation='horizontal', format=cbar_text_format,
                                 extend='max')
        else:
            cbar1 = fig.colorbar(gs_plot, cax=cbar_axes['gs'], orientation='horizontal', format=cbar_text_format)
        cbar1.ax.tick_params(labelsize=cbar_label_size)

    return fig


def bin_data(df, By_edges, Bz_edges):
    """
    :param df: pandas.DataFrame:
        Dataframe containing the data to be binned
    :param By_edges: numpy.array:
            The By bin edges, including the last bin
    :param Bz_edges: numpy.array:
            The Bz bin edges, including the last bin

    :return contour_data_is, contour_data_gs: 2d numpy.array, 2d numpy.array:
        Ionospheric and ground scatter binned data
    """

    n_bins_By = len(By_edges) - 1
    n_bins_Bz = len(Bz_edges) - 1

    delta_By = By_edges[1] - By_edges[0]
    delta_Bz = Bz_edges[1] - Bz_edges[0]

    contour_data_is = np.empty(shape=(n_bins_By, n_bins_Bz))
    contour_data_is[:] = math.nan

    contour_data_gs = np.empty(shape=(n_bins_By, n_bins_Bz))
    contour_data_gs[:] = math.nan

    for By_idx, By_start in enumerate(By_edges):
        if By_start == By_edges[-1]:
            continue  # The last edge is not a slice start

        By_end = By_start + delta_By
        df_By = df[(df['By_nT_GSM'] >= By_start) & (df['By_nT_GSM'] <= By_end)]

        for Bz_idx, Bz_start in enumerate(Bz_edges):
            if Bz_start == Bz_edges[-1]:
                continue  # The last edge is not a slice start

            Bz_end = Bz_start + delta_Bz
            df_By_Bz = df_By[(df_By['Bz_nT_GSM'] >= Bz_start) & (df_By['Bz_nT_GSM'] <= Bz_end)]

            try:
                contour_data_is[By_idx][Bz_idx] = sum(df_By_Bz['good_iono_echo']) / len(df_By_Bz)
                contour_data_gs[By_idx][Bz_idx] = sum(df_By_Bz['good_grndscat_echo']) / len(df_By_Bz)
            except ZeroDivisionError:
                # There are no points in this interval
                contour_data_is[By_idx][Bz_idx] = math.nan
                contour_data_gs[By_idx][Bz_idx] = math.nan
            except BaseException as e:
                print("By start: " + str(By_start))
                print("Bz start: " + str(Bz_start))
                raise e

    # # By default, Nan cells will be plotted white.  If this is not what is wanted, replace all nans with 0
    # contour_data_is = np.nan_to_num(contour_data_is, nan=0)
    # contour_data_gs = np.nan_to_num(contour_data_gs, nan=0)

    return contour_data_is, contour_data_gs


def assign_imf_values(df, imf_df):
    """

    We need to assign an IMF value to every row in SuperDARN df

    We will use the bisect method because apparently it is O(log n) on a sorted list:
    https://stackoverflow.com/questions/29700214/get-the-closest-datetime-from-a-list

    I am not really sure why this works but it seems to

    :param df: pandas.DataFrame: A SuperDARN dataframe that has a 'datetime' column
    :param imf_df: pandas.DataFrame: A dataframe with imf data for the same time period

    :return df: pandas.DataFrame: The provided SuperDARN dataframe, except now with the following additional columns:
                    'By_nT_GSM' (y-component of the IMF)
                    'Bz_nT_GSM' (z-component of the IMF)
    """

    By_nT_GSM, Bz_nT_GSM = [], []
    for i in range(len(df)):
        index = bisect(imf_df['datetime'], df['datetime'].iat[i], hi=len(imf_df['datetime']) - 1)

        # print("")
        # print("Looking for a matching date to " + str(df_mm['datetime'].iat[i]))
        # print("What we found was " + str(imf_df_mm['datetime'].iat[index]))

        By_nT_GSM.append(imf_df['By_nT_GSM'].iat[index])
        Bz_nT_GSM.append(imf_df['Bz_nT_GSM'].iat[index])

    df['By_nT_GSM'] = By_nT_GSM
    df['Bz_nT_GSM'] = Bz_nT_GSM

    # print("The dataframe here is of length: " + str(len(month_restricted_df)))

    # Remove all rows that have nan magnetic field data
    df = df.loc[df['By_nT_GSM'].notna() & df['Bz_nT_GSM'].notna()]

    # print("After removing nan IMF values, the dataframe here is of length: " + str(len(df_mm)))

    # print(df.head())

    return df


def title_figure(fig, station, year, time_units, beam_range, gate_range, freq_range, hour_range, day_range,
                 title_fontsize):
    """
    :param fig: The figure to title
    :param station: See occ_imf_variation() docstring
    :param year:  See occ_imf_variation() docstring
    :param time_units: See occ_imf_variation() docstring
    :param beam_range: See occ_imf_variation() docstring
    :param gate_range: See occ_imf_variation() docstring
    :param freq_range: See occ_imf_variation() docstring
    :param hour_range: See occ_imf_variation() docstring
    :param title_fontsize: Font title size
    """

    beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    gate_string = "Gates " + str(gate_range[0]) + "-" + str(gate_range[1])
    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"
    if day_range == (1, 31):
        day_string = "All Days"
    else:
        day_string = "Days " + str(day_range[0]) + "-" + str(day_range[1])
    if hour_range == (0, 24):
        hour_string = "All Hours"
    else:
        hour_string = "Hours " + str(hour_range[0]) + "-" + str(hour_range[1]) + " " + time_units.upper()

    fig.suptitle(str(year) + " at " + station.upper() + "; " + day_string + "; " + hour_string +
                 "\n" + gate_string + "; " + beam_string + "; " + freq_string +
                 "\nProduced by " + str(os.path.basename(__file__)), fontsize=title_fontsize)


def add_axes(fig):
    """
    :param fig: The figure to draw the axes on

    :return month_axes, cbar_axes: dictionary of axes, dictionary of axes
            Month axes will be keyed by month
    """

    # Leave a blank column in the middle to separate the is and gs subplots
    gs = fig.add_gridspec(ncols=9, nrows=7)

    month_axes, cbar_axes = dict(), dict()

    month_axes["January"] = {"is": fig.add_subplot(gs[0, 0:2]),
                             "gs": fig.add_subplot(gs[0, 5:7])}
    month_axes["February"] = {"is": fig.add_subplot(gs[0, 2:4]),
                              "gs": fig.add_subplot(gs[0, 7:9])}

    month_axes["March"] = {"is": fig.add_subplot(gs[1, 0:2]),
                           "gs": fig.add_subplot(gs[1, 5:7])}
    month_axes["April"] = {"is": fig.add_subplot(gs[1, 2:4]),
                           "gs": fig.add_subplot(gs[1, 7:9])}

    month_axes["May"] = {"is": fig.add_subplot(gs[2, 0:2]),
                         "gs": fig.add_subplot(gs[2, 5:7])}
    month_axes["June"] = {"is": fig.add_subplot(gs[2, 2:4]),
                          "gs": fig.add_subplot(gs[2, 7:9])}

    month_axes["July"] = {"is": fig.add_subplot(gs[3, 0:2]),
                          "gs": fig.add_subplot(gs[3, 5:7])}
    month_axes["August"] = {"is": fig.add_subplot(gs[3, 2:4]),
                            "gs": fig.add_subplot(gs[3, 7:9])}

    month_axes["September"] = {"is": fig.add_subplot(gs[4, 0:2]),
                               "gs": fig.add_subplot(gs[4, 5:7])}
    month_axes["October"] = {"is": fig.add_subplot(gs[4, 2:4]),
                             "gs": fig.add_subplot(gs[4, 7:9])}

    month_axes["November"] = {"is": fig.add_subplot(gs[5, 0:2]),
                              "gs": fig.add_subplot(gs[5, 5:7])}
    month_axes["December"] = {"is": fig.add_subplot(gs[5, 2:4]),
                              "gs": fig.add_subplot(gs[5, 7:9])}

    # # Add two colour bar axes, one for is and one for gs
    cbar_axes["is"] = fig.add_subplot(gs[6, 0:3])
    cbar_axes["gs"] = fig.add_subplot(gs[6, 5:8])

    return month_axes, cbar_axes


def apply_subplot_formatting(axes, x_lim, y_lim):
    """
    :param axes: matplotlib.axes: All of the smaller monthly axis requiring formatting
    :param x_lim: (int, int): x-axis limits
    :param y_lim: (int, int): y-axis limits
    """

    grid_colour = "black"

    label_font_size = 10
    title_font_size = 10

    for month, month_dict_item in axes.items():
        for subplot_type, ax in month_dict_item.items():
            ax.set_title(month[:3] + "; " + subplot_type.upper(), fontsize=title_font_size)
            ax.set_xlabel("By [nT] (GSM)", fontsize=label_font_size)
            ax.set_ylabel("Bz [nT] (GSM)", fontsize=label_font_size)

            ax.set_xlim(x_lim)
            ax.xaxis.set_major_locator(MultipleLocator(5))
            ax.xaxis.set_minor_locator(MultipleLocator(1))

            ax.set_ylim(y_lim)
            ax.yaxis.set_major_locator(MultipleLocator(5))
            ax.yaxis.set_minor_locator(MultipleLocator(1))

            ax.tick_params(axis='both', which='major', direction='in', color=grid_colour)
            ax.tick_params(axis='both', which='minor', direction='in', color=grid_colour)

            ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5, color=grid_colour)
            ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2, color=grid_colour)


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    if local_testing:
        station = "rkn"
        plot_type = 'contour'

        # Note: year, month, and day don't matter for local testing
        fig = occ_imf_variation(station=station, year=2011, day_range=(12, 12), hour_range=None,
                                gate_range=(10, 30), beam_range=(6, 8), freq_range=(11, 13),
                                plot_type=plot_type, local_testing=local_testing)

        plt.show()


    else:
        station = "dcn"
        even_odd_days = None
        years = [2019, 2020, 2021]
        plot_type = 'contour'
        freq_range = (8, 11)

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        for year in years:

            fig = occ_imf_variation(station=station, year=year, day_range=(1, 15), hour_range=None,
                                    gate_range=(10, 30), beam_range=(6, 8), freq_range=freq_range,
                                    local_testing=local_testing, plot_type=plot_type, even_odd_days=even_odd_days)

            out_fig = out_dir + "/occ_imf_variation_" + station + \
                      "_" + str(year) + "_" + str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz_" + plot_type

            print("Saving plot as " + out_fig)
            fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
