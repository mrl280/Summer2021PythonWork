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
      Occurrence is then plotted as contours or pixels.

    There one plot for every month, plus there are a couple of plots that summarize data for the whole year.

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

    vmaxs = {"is": 0.5,  # The top of the colour bar
             "gs": 0.05}
    title_fontsize = 18
    x_lim = (-20, 20)
    y_lim = (-20, 20)
    n_bins_By = 20
    n_bins_Bz = 20

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
    fig = plt.figure(figsize=[11, 26], constrained_layout=True, dpi=300)
    month_axes, year_axes = add_axes(fig=fig)
    apply_subplot_formatting(month_axes=month_axes, year_axes=year_axes, x_lim=x_lim, y_lim=y_lim, year=year)

    title_figure(fig=fig, station=station, year=year, time_units=time_units, beam_range=beam_range,
                 gate_range=gate_range, freq_range=freq_range, hour_range=hour_range, day_range=day_range,
                 title_fontsize=title_fontsize)

    print("     Adding month to df...")
    month = []
    for i in range(len(df)):
        month.append(df['datetime'].iat[i].month)
    df['month'] = month

    # Compute By edges (x-axis edges)
    By_edges = np.linspace(x_lim[0], x_lim[1], num=(n_bins_By + 1))

    # Compute Bz edges (y-axis edges)
    Bz_edges = np.linspace(y_lim[0], y_lim[1], num=(n_bins_Bz + 1))

    # We have to complete IMF assignment for the whole dataframe, even if it takes longer it is simplest to do the whole
    # year all at once
    print("     Assigning IMF values for " + str(year) + "...")
    df = assign_imf_values(df=df, imf_df=imf_df)


    """ Complete the large yearly contour plots"""
    print("     Computing binned occ rates for the whole year...")
    contour_data_is, contour_data_gs = bin_data(df=df, By_edges=By_edges, Bz_edges=Bz_edges)

    print("     Plotting data for the whole year...")
    plot_data(fig=fig, axes=year_axes['whole_year'], contour_data_is=contour_data_is, contour_data_gs=contour_data_gs,
              By_edges=By_edges, Bz_edges=Bz_edges, plot_type=plot_type, cbar=True, vmaxs=vmaxs)


    """ Complete the histogram plots"""
    complete_hist_plot(axes=year_axes['By_hist'], df=df, param='By_nT_GSM', edges=By_edges)
    complete_hist_plot(axes=year_axes['Bz_hist'], df=df, param='Bz_nT_GSM', edges=By_edges)

    # Compute magnetic field magnitude and plot it on the bigger histogram plot
    df['B_magnitude'] = np.sqrt(np.square(df['By_nT_GSM']) + np.square(df['Bz_nT_GSM']))

    # We are going to need a new set of edges here because we won't have any negative magnitudes
    magnitude_edges = np.linspace(0, x_lim[1], num=int(n_bins_By + 1))
    complete_hist_plot(axes=year_axes['big_hist'], df=df, param='B_magnitude', edges=magnitude_edges)


    """ Complete all of the small monthly axes """
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

        print("     Computing binned occ rates for " + month_str + "...")
        contour_data_is, contour_data_gs = bin_data(df=df_mm, By_edges=By_edges, Bz_edges=Bz_edges)

        print("     Plotting data for " + month_str + "...")
        plot_data(fig=fig, axes=month_axis, contour_data_is=contour_data_is, contour_data_gs=contour_data_gs,
                  By_edges=By_edges, Bz_edges=Bz_edges, plot_type=plot_type, cbar=False, vmaxs=vmaxs)

    return fig


def complete_hist_plot(axes, df, param, edges):
    """
    :param axes: dictionary of matplotlib.axes:
            The set of axes to plot on
    :param df: pandas.DataFrame:
            SuperDARN occurrence dataframe
    :param param: str:
            The df parameter to use for binning - this should be what is plotted along the x-axis
    :param edges: numpy.array:
            The histogram edges
    """

    # We need to bin data along a single direction
    n_bins = len(edges) - 1
    delta = edges[1] - edges[0]

    hist_data_is = np.empty(shape=n_bins)
    hist_data_is[:] = math.nan

    hist_data_gs = np.empty(shape=n_bins)
    hist_data_gs[:] = math.nan

    for idx, start in enumerate(edges):
        if start == edges[-1]:
            continue  # The last edge is not a slice start

        end = start + delta
        df_rr = df[(df[param] >= start) & (df[param] <= end)]

        try:
            hist_data_is[idx] = sum(df_rr['good_iono_echo']) / len(df_rr)
            hist_data_gs[idx] = sum(df_rr['good_grndscat_echo']) / len(df_rr)
        except ZeroDivisionError:
            # There are no points in this interval
            hist_data_is[idx] = 0
            hist_data_gs[idx] = 0
        except BaseException as e:
            print("start: " + str(start))
            raise e

    bin_centers = edges[1:] - delta / 2
    axes['is'].hist(bin_centers, weights=hist_data_is, histtype='step', align='mid', linewidth=1)
    axes['gs'].hist(bin_centers, weights=hist_data_gs, histtype='step', align='mid', linewidth=1)


def plot_data(fig, axes, contour_data_is, contour_data_gs, By_edges, Bz_edges, plot_type, vmaxs, cbar=True):
    """
    :param axes: dictionary of matplotlib.axes:
            The set of axes to plot on
    :param contour_data_is: 2d numpy.array:
            Ionospheric binned data
    :param contour_data_gs: 2d numpy.array:
            Ground scatter binned data
    :param By_edges: numpy.array:
            The histogram edges along the x-direction
    :param Bz_edges: numpy.array:
            The histogram edges along the y-direction
    :param plot_type: str:
            The type of plot, either 'contour' or 'pixel'
    :param vmaxs:

    :param cbar: bool (optional; default is True):
            Whether or not to add a colour bar to plot
    """

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
            is_plot = axes['is'].contourf(bin_xcenters, bin_ycenters, contour_data_is.transpose(),
                                          cmap=cmap, levels=is_levels, zorder=0, extend='max')
        else:
            is_plot = axes['is'].contourf(bin_xcenters, bin_ycenters, contour_data_is.transpose(),
                                          cmap=cmap, levels=is_levels, zorder=0)

        if vmaxs['gs'] < 1:
            gs_plot = axes['gs'].contourf(bin_xcenters, bin_ycenters, contour_data_gs.transpose(),
                                          cmap=cmap, levels=gs_levels, zorder=0, extend='max')
        else:
            gs_plot = axes['gs'].contourf(bin_xcenters, bin_ycenters, contour_data_gs.transpose(),
                                          cmap=cmap, levels=gs_levels, zorder=0)

    elif plot_type == "pixel":

        cmap = 'jet'
        is_plot = axes['is'].pcolormesh(By_edges, Bz_edges, contour_data_is.transpose(),
                                        cmap=cmap, vmin=0, vmax=vmaxs['is'], zorder=0)
        gs_plot = axes['gs'].pcolormesh(By_edges, Bz_edges, contour_data_gs.transpose(),
                                        cmap=cmap, vmin=0, vmax=vmaxs['gs'], zorder=0)

    else:
        raise Exception("plot_type not recognized")

    if cbar is True:

        cbar_label_size = 8

        # The extend here is needed for the pixel option, it has no effect on contours,
        #  contour extend was set when plotting
        if vmaxs['is'] < 1:
            is_cbar = fig.colorbar(is_plot, ax=axes['is'], format='%.2f',
                                   extend='max')
        else:
            is_cbar = fig.colorbar(is_plot, ax=axes['is'], format='%.2f')
        is_cbar.ax.tick_params(labelsize=cbar_label_size)

        if vmaxs['gs'] < 1:
            cbar1 = fig.colorbar(gs_plot, ax=axes['gs'], format='%.3f',
                                 extend='max')
        else:
            cbar1 = fig.colorbar(gs_plot, ax=axes['gs'], format='%.3f')
        cbar1.ax.tick_params(labelsize=cbar_label_size)


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
        datetime_here = df['datetime'].iat[i]

        # Some matches will be off by up to 30 seconds
        index = bisect(imf_df['datetime'], datetime_here, hi=len(imf_df['datetime'])) - 1

        # print("")
        # print("Looking for a matching date to " + str(datetime_here))
        # print("What we found was " + str(imf_df['datetime'].iat[index]))

        # There could be areas missing in the IMF data - only assign a value here if the match is acceptably close
        if math.fabs((datetime_here - imf_df['datetime'].iat[index]).total_seconds()) <= 180:
            By_nT_GSM.append(imf_df['By_nT_GSM'].iat[index])
            Bz_nT_GSM.append(imf_df['Bz_nT_GSM'].iat[index])
        else:
            By_nT_GSM.append(np.nan)
            Bz_nT_GSM.append(np.nan)

    df['By_nT_GSM'] = By_nT_GSM
    df['Bz_nT_GSM'] = Bz_nT_GSM

    # print("The dataframe here is of length: " + str(len(month_restricted_df)))

    # Remove all rows that have nan magnetic field data
    df = df.loc[df['By_nT_GSM'].notna() & df['Bz_nT_GSM'].notna()]
    df.reset_index(drop=True, inplace=True)

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

    gs = fig.add_gridspec(ncols=8, nrows=12)

    # Will will have axes for monthly reporting, and axes for yearly reporting
    month_axes, year_axes = dict(), dict()

    # There will be one set of contour plots for each month
    month_axes["January"] = {"is": fig.add_subplot(gs[0, 0:2]),
                             "gs": fig.add_subplot(gs[0, 4:6])}
    month_axes["February"] = {"is": fig.add_subplot(gs[0, 2:4]),
                              "gs": fig.add_subplot(gs[0, 6:8])}

    month_axes["March"] = {"is": fig.add_subplot(gs[1, 0:2]),
                           "gs": fig.add_subplot(gs[1, 4:6])}
    month_axes["April"] = {"is": fig.add_subplot(gs[1, 2:4]),
                           "gs": fig.add_subplot(gs[1, 6:8])}

    month_axes["May"] = {"is": fig.add_subplot(gs[2, 0:2]),
                         "gs": fig.add_subplot(gs[2, 4:6])}
    month_axes["June"] = {"is": fig.add_subplot(gs[2, 2:4]),
                          "gs": fig.add_subplot(gs[2, 6:8])}

    month_axes["July"] = {"is": fig.add_subplot(gs[3, 0:2]),
                          "gs": fig.add_subplot(gs[3, 4:6])}
    month_axes["August"] = {"is": fig.add_subplot(gs[3, 2:4]),
                            "gs": fig.add_subplot(gs[3, 6:8])}

    month_axes["September"] = {"is": fig.add_subplot(gs[4, 0:2]),
                               "gs": fig.add_subplot(gs[4, 4:6])}
    month_axes["October"] = {"is": fig.add_subplot(gs[4, 2:4]),
                             "gs": fig.add_subplot(gs[4, 6:8])}

    month_axes["November"] = {"is": fig.add_subplot(gs[5, 0:2]),
                              "gs": fig.add_subplot(gs[5, 4:6])}
    month_axes["December"] = {"is": fig.add_subplot(gs[5, 2:4]),
                              "gs": fig.add_subplot(gs[5, 6:8])}

    # There will be one set of big contour plots for the whole year
    year_axes["whole_year"] = {"is": fig.add_subplot(gs[6:9, 0:4]),
                               "gs": fig.add_subplot(gs[6:9, 4:8])}

    # There will be a couple of histogram plots for the whole year
    year_axes["By_hist"] = {"is": fig.add_subplot(gs[9, 0:2]),
                            "gs": fig.add_subplot(gs[9, 4:6])}
    year_axes["Bz_hist"] = {"is": fig.add_subplot(gs[9, 2:4]),
                            "gs": fig.add_subplot(gs[9, 6:8])}

    year_axes["big_hist"] = {"is": fig.add_subplot(gs[10:12, 0:4]),
                             "gs": fig.add_subplot(gs[10:12, 4:8])}

    return month_axes, year_axes


def apply_subplot_formatting(month_axes, year_axes, x_lim, y_lim, year):
    """
    :param month_axes: dictionary of matplotlib.axes:
            All of the smaller monthly axis requiring formatting
    :param year_axes: dictionary of matplotlib.axes:
            All of the year axes that require formatting
    :param x_lim: (int, int):
            x-axis limits
    :param y_lim: (int, int):
            y-axis limits
    :param year: int:
            The year to consider.
    """

    grid_colour = "black"

    label_font_size = 10
    title_font_size = 10

    for month_dict_key, month_dict_item in month_axes.items():
        for subplot_type, ax in month_dict_item.items():
            ax.set_title(month_dict_key[:3] + "; " + subplot_type.upper(), fontsize=title_font_size)
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
            # ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2, color=grid_colour)

    for year_dict_key, year_dict_item in year_axes.items():
        if year_dict_key == "whole_year":
            # Then we are formatting the larger contour plots that represent the whole year
            for subplot_type, ax in year_dict_item.items():
                ax.set_title("All of " + str(year) + "; " + subplot_type.upper(), fontsize=title_font_size)
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

        else:
            # We are formatting the histogram plots
            for subplot_type, ax in year_dict_item.items():
                ax.set_title("All of " + str(year) + "; " + subplot_type.upper(), fontsize=title_font_size)
                ax.set_ylabel("Occurrence", fontsize=label_font_size)

                ax.set_xlim(x_lim)
                ax.xaxis.set_major_locator(MultipleLocator(5))
                ax.xaxis.set_minor_locator(MultipleLocator(1))

                if year_dict_key == "By_hist":
                    ax.set_xlabel("By [nT] (GSM)", fontsize=label_font_size)

                elif year_dict_key == "Bz_hist":
                    ax.set_xlabel("Bz [nT] (GSM)", fontsize=label_font_size)

                elif year_dict_key == "big_hist":
                    ax.set_xlabel("sqrt(By^2 + Bz^2) [nT] (GSM)", fontsize=label_font_size)
                    ax.set_xlim([0, x_lim[1]])  # We won't have any negative values here


if __name__ == '__main__':
    """ Testing """

    local_testing = False

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
