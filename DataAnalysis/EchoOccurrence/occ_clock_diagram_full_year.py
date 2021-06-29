import math
import os
import pathlib
import pydarn

from cartopy.util import add_cyclic_point
from matplotlib import pyplot as plt
from pydarn import radar_fov, SuperDARNRadars

import cartopy.crs as ccrs
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import numpy as np

from lib.cm.modified_jet import modified_jet
from lib.add_mlt_to_df import add_mlt_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def occ_clock_diagram_full_year(station, year, day_range=None, gate_range=None, beam_range=None,
                                freq_range=None, plot_type='pixel', local_testing=False):
    """

    Produce a full years worth of full circle stereographic occurrence plots.
    All monthly and seasonal plots, all on one page

    If you would like a single clock diagram produced for an arbitrary amount of time, please see occ_clock_diagram.

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only considers 45 km data.
            (a warning will be printed if other spatial resolution data is stripped from the dataset)
        - To check which fitACF program is being used, refer to the data readers in lib.data_getters

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year: int:
            The year to consider.
    :param day_range: (<int>, <int>) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param gate_range: (<int>, <int>) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param freq_range: (<float>, <float>) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.
    :param plot_type: str (optional):
            The type of plot, either 'contour' or 'pixel', default is 'contour'
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.

    :return: pandas.DataFrame, matplotlib.pyplot.figure: The dataframe and the figure. The dataframe can then be
    analyzed, and the figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    time_units = "mlt"  # TODO: Add compatibility with other time units, if it makes sense

    year = check_year(year)
    freq_range = check_freq_range(freq_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"

    print("     Retrieving data...")
    df = get_data_handler(station, year_range=(year, year), month_range=None, day_range=day_range,
                          gate_range=gate_range, beam_range=beam_range, freq_range=freq_range, occ_data=True,
                          local_testing=local_testing)
    df = only_keep_45km_res_data(df)

    print("     Computing MLTs...")
    # Use the middle of the year as a magnetic field estimate
    date_time_est, _ = build_datetime_epoch(year=year, month=6, day=15, hour=0)
    cell_corners_aacgm_lats, cell_corners_aacgm_lons = radar_fov(stid=radar_id, coords='aacgm', date=date_time_est)

    df = add_mlt_to_df(cell_corners_aacgm_lons=cell_corners_aacgm_lons, cell_corners_aacgm_lats=cell_corners_aacgm_lats,
                       df=df)

    # Right now MLTs are in the range 0-24, we need to put it in the range 0-360 for circular plotting
    df['mlt'] = 15 * df['mlt']

    # Since we are using the North Pole projection we need all latitudes to be positive
    df['lat'] = df['lat'].abs()

    # Compute month, we will need it to filter the data
    month = []
    for i in range(len(df)):
        month.append(df['datetime'].iat[i].month)
    df['month'] = month

    print("     Preparing the figure...")
    fig = plt.figure(figsize=[8, 14], constrained_layout=True, dpi=300)

    lat_extreme = 59  # TODO: Adjust lat extreme based on the radar (use the most extreme point in the fan)
    month_axes, season_axes, cbar_axis = add_axes(fig=fig, radar_info=this_radars_info, lat_extreme=lat_extreme)
    fig.suptitle(str(year) + " at " + station.upper() + "; " + beam_string + "; " + freq_string +
                 "\n Data Plotted in AACGM Latitudes and " + time_units.upper() +
                 "\nProduced by " + str(os.path.basename(__file__)), fontsize=18)

    # Compute mlt edges
    deg_mlt_per_bin = 2
    n_bins_mlt = int(360 / deg_mlt_per_bin)
    mlt_edges = np.linspace(0, 360, num=(n_bins_mlt + 1))

    # Compute latitude edges
    deg_lat_per_bin = 1
    n_bins_lat = int((90 - lat_extreme) / deg_lat_per_bin)
    lat_edges = np.linspace(lat_extreme, 90, num=(n_bins_lat + 1))

    # Add the month data to the plot
    for month_str in month_axes:
        print("     Computing occurrence data for " + month_str)

        month_datetime_object = datetime.datetime.strptime(month_str, "%B")
        month = month_datetime_object.month

        df_mm = df[df['month'] == month]
        ax = month_axes[month_str]

        add_data_to_plot(df=df_mm, ax=ax, mlt_edges=mlt_edges, lat_edges=lat_edges, plot_type=plot_type)

    # Add the seasonal data to the plot
    for season in season_axes:
        print("     Computing occurrence data for " + season)

        seasonal_month_ranges = get_seasonal_month_ranges()
        restriction_range = seasonal_month_ranges[season]

        if season == "Winter":
            # Note: winter time requires a unique filter because we need months outside this range
            df_ss = df[(df['month'] >= restriction_range[0]) | (df['month'] <= restriction_range[1])]
        else:
            df_ss = df[(df['month'] >= restriction_range[0]) & (df['month'] <= restriction_range[1])]

        ax = season_axes[season]

        cbar_ref_plot, _ = add_data_to_plot(df=df_ss, ax=ax, mlt_edges=mlt_edges, lat_edges=lat_edges,
                                            plot_type=plot_type)

    add_colour_bar(fig=fig, cax=cbar_axis, plot=cbar_ref_plot)

    return df, fig


def add_data_to_plot(df, ax, mlt_edges, lat_edges, plot_type='pixel'):
    """
    Add clock data to a stereographic plot
    :param df: pandas.DataFrame: The restricted dataframe containing the occurance data
    :param ax: matplotlib.axes: The axes to draw on
    :param mlt_edges: The MLT (x data) edges
    :param lat_edges: The LAT (y data) edges
    :param plot_type: str: The type of plot - either "pixel" (default) or "contour"
    :return matplotlib.pyplot.plot, matplotlib.pyplot.plot: The ionospheric and ground scatter plots
    """

    n_bins_mlt = len(mlt_edges) - 1
    n_bins_lat = len(lat_edges) - 1

    delta_mlt = mlt_edges[1] - mlt_edges[0]
    delta_lat = lat_edges[1] - lat_edges[0]

    # Loop through the plot cells and compute occurrence rates
    data_is = np.empty(shape=(n_bins_mlt, n_bins_lat))
    data_is[:] = math.nan

    data_gs = np.empty(shape=(n_bins_mlt, n_bins_lat))
    data_gs[:] = math.nan

    for mlt_idx, start_mlt in enumerate(mlt_edges):
        if start_mlt == mlt_edges[-1]:
            continue  # The last edge is not a start
        end_mlt = start_mlt + delta_mlt
        df_mlt = df[(df['mlt'] >= start_mlt) & (df['mlt'] <= end_mlt)]

        for lat_idx, start_lat in enumerate(lat_edges):
            if start_lat == lat_edges[-1]:
                continue  # The last edge is not a start

            end_lat = start_lat + delta_lat
            df_mlt_lat = df_mlt[(df_mlt['lat'] >= start_lat) & (df_mlt['lat'] <= end_lat)]

            try:
                data_is[mlt_idx][lat_idx] = sum(df_mlt_lat['good_iono_echo']) / len(df_mlt_lat)
                data_gs[mlt_idx][lat_idx] = sum(df_mlt_lat['good_grndscat_echo']) / len(df_mlt_lat)
            except ZeroDivisionError:
                # There are no point in this interval
                data_is[mlt_idx][lat_idx] = math.nan
                data_gs[mlt_idx][lat_idx] = math.nan
            except BaseException as e:
                print("MLT index: " + str(mlt_idx))
                print("LAT index: " + str(lat_idx))
                raise e

    is_plot, gs_plot = plot_data_on_axis(axes=ax, data_is=data_is, data_gs=data_gs,
                                         mlt_edges=mlt_edges, lat_edges=lat_edges, plot_type=plot_type)

    # All plots have the same scale, we just need to return one for the colour bar to reference
    return is_plot, gs_plot


def plot_data_on_axis(axes, data_is, data_gs, mlt_edges, lat_edges, plot_type='pixel'):
    """
    Add the ionospheric and ground scatter data to the plots
    :param axes: dictionary of matplotlib.axes: The axes to draw on
    :param data_is: 2d array containing the ionospheric data
    :param data_gs: 2d array containing the ground scatter data
    :param mlt_edges: The MLT (x data) edges
    :param lat_edges: The LAT (y data) edges
    :param plot_type: str: The type of plot - either "pixel" (default) or "contour"
    :return matplotlib.pyplot.plot, matplotlib.pyplot.plot: The ionospheric and ground scatter plots
    """

    if plot_type == "contour":

        # Compute bin centers
        bin_xwidth = (mlt_edges[1] - mlt_edges[0])
        bin_ywidth = (lat_edges[1] - lat_edges[0])
        bin_xcenters = mlt_edges[1:] - bin_xwidth / 2
        bin_ycenters = lat_edges[1:] - bin_ywidth / 2

        # Adding a cyclic points is required to complete the circle
        # Without adding in this point, there will be one pie-shaped piece missing from the circle
        contour_data_is_cyclic, bin_xcenters_cyclic = add_cyclic_point(data_is.transpose(), coord=bin_xcenters)
        contour_data_gs_cyclic, bin_xcenters_cyclic = add_cyclic_point(data_gs.transpose(), coord=bin_xcenters)

        levels = 12
        levels = np.linspace(start=0, stop=1, num=(levels + 1))
        cmap = "jet"
        # cmap = modified_jet(levels=len(levels) - 1)

        is_plot = axes["is"].contourf(bin_xcenters_cyclic, bin_ycenters, contour_data_is_cyclic,
                                      cmap=cmap, levels=levels, transform=ccrs.PlateCarree())
        gs_plot = axes["gs"].contourf(bin_xcenters_cyclic, bin_ycenters, contour_data_gs_cyclic,
                                      cmap=cmap, levels=levels, transform=ccrs.PlateCarree())

    elif plot_type == "pixel":

        cmap = 'jet'
        is_plot = axes["is"].pcolormesh(mlt_edges, lat_edges, data_is.transpose(),
                                        transform=ccrs.PlateCarree(), cmap=cmap, vmin=0, vmax=1)
        gs_plot = axes["gs"].pcolormesh(mlt_edges, lat_edges, data_gs.transpose(),
                                        transform=ccrs.PlateCarree(), cmap=cmap, vmin=0, vmax=1)

    else:
        raise Exception("plot_type not recognized")

    return is_plot, gs_plot


def add_colour_bar(fig, cax, plot):
    """
    Add a colour bar to the the provided fig/cax with reference to the provided plot
    """
    cbar = fig.colorbar(plot, cax=cax, orientation="horizontal", format='%.1f')
    cbar.ax.tick_params(labelsize=18)


def get_seasonal_month_ranges():
    """
    :return: dictionary of int tuples representing suitable month ranges for each season
    """
    return {'Spring': (3, 5),
            'Summer': (6, 8),
            'Autumn': (9, 11),

            # Note: winter time requires a unique filter because we need months outside this range
            'Winter': (12, 2)
            }


def add_axes(fig, radar_info, lat_extreme=60):
    """
    :param lat_extreme: The most extreme latitude to plot. Default is 60 deg.
        Note: since we always use a NorthPoleStereo projection, this should be positive
    :param radar_info: The radars info
    :param fig: The figure to draw the axes on
    :return: numpy array of
    """

    # Regardless of hemisphere, we will always use the North Pole projection because it gives us a stereo projection
    # with the zero degree line at the bottom of the plot and East is CCW (this is what we need for a clock diagram).
    proj = ccrs.NorthPolarStereo()

    gs = fig.add_gridspec(ncols=8, nrows=13)

    month_axes = dict()
    season_axes = dict()

    # Build all of the spring axes
    # For seasonal axis, remember that splices don't include last index
    month_axes["March"] = {"is": fig.add_subplot(gs[0, 0], projection=proj),
                           "gs": fig.add_subplot(gs[0, 7], projection=proj)}
    month_axes["April"] = {"is": fig.add_subplot(gs[1, 0], projection=proj),
                           "gs": fig.add_subplot(gs[1, 7], projection=proj)}
    month_axes["May"] = {"is": fig.add_subplot(gs[2, 0], projection=proj),
                         "gs": fig.add_subplot(gs[2, 7], projection=proj)}
    season_axes["Spring"] = {"is": fig.add_subplot(gs[0:3, 1:4], projection=proj),
                             "gs": fig.add_subplot(gs[0:3, 4:7], projection=proj)}

    # Build all of the summer time axes
    month_axes["June"] = {"is": fig.add_subplot(gs[3, 0], projection=proj),
                          "gs": fig.add_subplot(gs[3, 7], projection=proj)}
    month_axes["July"] = {"is": fig.add_subplot(gs[4, 0], projection=proj),
                          "gs": fig.add_subplot(gs[4, 7], projection=proj)}
    month_axes["August"] = {"is": fig.add_subplot(gs[5, 0], projection=proj),
                            "gs": fig.add_subplot(gs[5, 7], projection=proj)}
    season_axes["Summer"] = {"is": fig.add_subplot(gs[3:6, 1:4], projection=proj),
                             "gs": fig.add_subplot(gs[3:6, 4:7], projection=proj)}

    # Build all of the fall time axes
    month_axes["September"] = {"is": fig.add_subplot(gs[6, 0], projection=proj),
                               "gs": fig.add_subplot(gs[6, 7], projection=proj)}
    month_axes["October"] = {"is": fig.add_subplot(gs[7, 0], projection=proj),
                             "gs": fig.add_subplot(gs[7, 7], projection=proj)}
    month_axes["November"] = {"is": fig.add_subplot(gs[8, 0], projection=proj),
                              "gs": fig.add_subplot(gs[8, 7], projection=proj)}
    season_axes["Autumn"] = {"is": fig.add_subplot(gs[6:9, 1:4], projection=proj),
                             "gs": fig.add_subplot(gs[6:9, 4:7], projection=proj)}

    # Build all of the winter time axes
    month_axes["December"] = {"is": fig.add_subplot(gs[9, 0], projection=proj),
                              "gs": fig.add_subplot(gs[9, 7], projection=proj)}
    month_axes["January"] = {"is": fig.add_subplot(gs[10, 0], projection=proj),
                             "gs": fig.add_subplot(gs[10, 7], projection=proj)}
    month_axes["February"] = {"is": fig.add_subplot(gs[11, 0], projection=proj),
                              "gs": fig.add_subplot(gs[11, 7], projection=proj)}
    season_axes["Winter"] = {"is": fig.add_subplot(gs[9:12, 1:4], projection=proj),
                             "gs": fig.add_subplot(gs[9:12, 4:7], projection=proj)}

    cbar_axis = fig.add_subplot(gs[12, 1:7])

    apply_common_subplot_formatting(month_axes=month_axes, season_axes=season_axes,
                                    radar_info=radar_info, lat_extreme=lat_extreme)

    return month_axes, season_axes, cbar_axis


def apply_common_subplot_formatting(month_axes, season_axes, radar_info, lat_extreme):
    """
    :param lat_extreme: The most extreme latitude to plot.
        Note: since we always use a NorthPoleStereo projection, this should be positive
    :param month_axes: matplotlib.axes: All of the smaller monthly axis to draw on
    :param season_axes: matplotlib.axes: All of the larger seasonal axis to draw on
    :param radar_info: The radars info
    """

    subplot_types = ["is", "gs"]

    hemisphere = radar_info.hemisphere

    # Compute a circle in axis coordinates which can be used as a boundary
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5  # 0.5 means middle of the circle
    vertices = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(vertices * radius + center)

    # Apply common subplot formatting
    all_axes = {**month_axes, **season_axes}
    for _, ax in all_axes.items():
        for subplot_type in subplot_types:
            ax[subplot_type].set_extent([-180, 180, 90, lat_extreme], crs=ccrs.PlateCarree())
            ax[subplot_type].set_boundary(circle, transform=ax[subplot_type].transAxes)

            # Add gridlines   # Note: Labels wont draw on a circular axis
            gl = ax[subplot_type].gridlines(draw_labels=True, linestyle='--', linewidth=0.5, zorder=5)
            gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135])

    # Add printouts to month axes
    for month in month_axes:
        for subplot_type in subplot_types:
            ax = month_axes[month][subplot_type]

            # Print out the month
            ax.text(ax.get_xlim()[0], ax.get_ylim()[1], month[:3], ha='center', va='top', fontsize=9)

            # Print out the type of scatter
            ax.text(ax.get_xlim()[1], ax.get_ylim()[1], subplot_type.upper(), ha='center', va='top', fontsize=9)

    # Add printouts to seasonal axes
    for season in season_axes:
        for subplot_type in subplot_types:
            ax = season_axes[season][subplot_type]

            if hemisphere.value == 1:
                season_str = season
            elif hemisphere.value == -1:
                # The seasons are reversed in the southern hemisphere
                season_str = south_hemi_season(north_hemi_season=season)
            else:
                raise Exception("Error: Season not recognized.")

            # Print out the season
            ax.text(ax.get_xlim()[0], ax.get_ylim()[1], season_str, ha='left', va='top', fontsize=12)

            # Print out the type of scatter
            ax.text(ax.get_xlim()[1], ax.get_ylim()[1], subplot_type.upper(), ha='right', va='top', fontsize=12)

            # Print clock numbers
            text_offset_multiplier = 1.03
            ax.text(0, text_offset_multiplier * ax.get_ylim()[1], "12", ha='center', va='bottom')
            ax.text(0, text_offset_multiplier * ax.get_ylim()[0], "00", ha='center', va='top')
            ax.text(text_offset_multiplier * ax.get_xlim()[1], 0, "06", ha='left', va='center')
            ax.text(text_offset_multiplier * ax.get_xlim()[0], 0, "18", ha='right', va='center')


def south_hemi_season(north_hemi_season):
    """
    :param north_hemi_season: str: the season in the Northern hemisphere
    :return: str: the corresponding southern hemisphere season
    """
    if north_hemi_season.lower() == "spring":
        return "Autumn"
    elif north_hemi_season.lower() == "summer":
        return "Winter"
    elif north_hemi_season.lower() == "autumn":
        return "Spring"
    elif north_hemi_season.lower() == "winter":
        return "Summer"
    else:
        raise Exception("North hemisphere season of " + str(north_hemi_season) + " not recognized.")


if __name__ == '__main__':
    """ Testing """

    local_testing = False

    if local_testing:
        station = "dce"

        _, fig = occ_clock_diagram_full_year(station=station, year=2011, day_range=None,
                                             gate_range=(0, 74), beam_range=(6, 7), freq_range=None,
                                             plot_type='pixel', local_testing=local_testing)

        plt.show()


    else:
        station = "dcn"
        freq_range = (8, 10)
        plot_type = 'pixel'

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        for year in range(2019, 2022, 1):
            _, fig = occ_clock_diagram_full_year(station=station, year=year, day_range=None,
                                                 gate_range=(0, 74), beam_range=None, freq_range=freq_range,
                                                 plot_type=plot_type, local_testing=local_testing)

            out_fig = out_dir + "/occ_clock_diagram_full_year_" + station + "_" + str(year) + "_" + \
                      str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz" + "_" + plot_type
            print("Saving plot as " + out_fig)
            fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
