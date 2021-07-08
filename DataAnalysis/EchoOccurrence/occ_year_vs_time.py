import math
import os
import pathlib
import pydarn

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import SuperDARNRadars

from lib.get_middle_of_fov import get_middle_of_fov
from lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def occ_year_vs_time(station, year_range, month_range=None, day_range=None, hour_range=None,
                     gate_range=None, beam_range=None, freq_range=None,
                     time_units='mlt', plot_type='contour', angle_contours=None,
                     local_testing=False):
    """

    Produce plots with year on the y-axis and time along the x-axis.
    There are 2 subplots, one for ionospheric scatter and one for ground scatter
    Echo occurrence rates are plotted along z (as either contours or pixels)

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only considers 45 km data
        - To check which fitACF program is being used, refer to the data readers in lib.data_getters

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year_range: (int, int):
            Inclusive. The year range to consider.
    :param month_range: (int, int) (optional):
            Inclusive. The month range to consider.  If omitted (or None), then all months will be considered.
    :param day_range: (int, int) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param hour_range: (int, int) (optional):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not quite inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param gate_range: (int, int) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (int, int) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param freq_range: (float, float) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.
    :param time_units: str:
            The time units to plot along x.  Default is 'mlt'.
                'ut' for universal time
                'mlt' for magnetic local time
                'lt' for local time (based on longitude)
                'lst' for local standard time (based on time zones)
                'ast' for apparent solar time (based on the apparent angular motion of the sun across the sky)
    :param plot_type: str (optional): default is 'contour'.
            The type of plot, either 'contour' or 'pixel'.
    :param angle_contours: str (optional):
            The type of angle contours to superimpose on top of the data.
                'zenith' for solar zenith angle (the angle between the sun’s rays and a vertical plane)
                'altitude' for solar altitude angle (the angle between the sun’s rays and a horizontal plane)
    :param local_testing: bool (optional): default is False.
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :return: pandas.DataFrame, matplotlib.pyplot.figure
            The dataframe used and the figure created.
            The figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    vmax = 0.8  # The top of the colour bar

    time_units = check_time_units(time_units)
    year_range = check_year_range(year_range)
    day_range = check_day_range(day_range)
    hour_range = check_hour_range(hour_range)
    month_range = check_month_range(month_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    print("Retrieving data...")
    df = get_data_handler(station, year_range=year_range, month_range=month_range, day_range=day_range,
                          gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                          occ_data=True, local_testing=local_testing)
    df = only_keep_45km_res_data(df)

    # Add decimal hour to df in whatever units were requested
    # Use the middle of the mid year as magnetic field estimate
    mid_year = int(year_range[0] + (year_range[1] - year_range[0]) / 2)
    date_time_est, _ = build_datetime_epoch(year=mid_year, month=6, day=15, hour=0)
    df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id, date_time_est=date_time_est)

    df = df.loc[(df[time_units] >= hour_range[0]) & (df[time_units] <= hour_range[1])]
    df.reset_index(drop=True, inplace=True)

    print("Computing decimal year...")
    months_in_a_year = 12
    days_in_a_year = 365
    hours_in_a_year = 8760

    decimal_years = []
    for i in range(len(df)):
        datetime_obj_here = df['datetime'].iat[i]
        decimal_hour_here = df[time_units].iat[i]

        decimal_year_here = datetime_obj_here.year + (datetime_obj_here.month - 1) / months_in_a_year + \
                            (datetime_obj_here.day - 1) / days_in_a_year + decimal_hour_here / hours_in_a_year

        decimal_years.append(decimal_year_here)

    df['decimal_year'] = np.asarray(decimal_years)

    if day_range == (1, 31):
        day_string = "All Days"
    else:
        day_string = "Days " + str(day_range[0]) + "-" + str(day_range[1])
    beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"

    # Prepare the plot
    fig, ax = plt.subplots(figsize=[10, 12], dpi=300, nrows=2, ncols=1, constrained_layout=True)
    fig.suptitle(station.upper() + "; " + beam_string + "; " + freq_string + "; " + day_string +
                 "\nProduced by " + str(os.path.basename(__file__)), fontsize=18)

    ax[0].set_title("Ionospheric Scatter", fontsize=20)
    ax[1].set_title("Ground Scatter", fontsize=20)

    if plot_type == "contour":
        tick_colour = "black"
    else:
        tick_colour = "white"

    # Apply common subplot formatting
    for i in range(ax.size):

        ax[i].set_xlim(hour_range)
        ax[i].set_xlabel("Time, " + time_units.upper(), fontsize=16)
        ax[i].xaxis.set_major_locator(MultipleLocator(4))

        ax[i].set_ylim((year_range[0], year_range[1] + 1))
        ax[i].set_ylabel("Year", fontsize=16)
        ax[i].yaxis.set_major_locator(MultipleLocator(1))
        ax[i].yaxis.set_minor_locator(MultipleLocator(0.25))

        ax[i].tick_params(axis='both', which='major', direction='in', color=tick_colour, labelsize=14)
        ax[i].tick_params(axis='y', which='minor', direction='in', color=tick_colour, labelsize=14)

    # Compute hour edges
    bins_per_hour = 4  # quarter hour bins
    n_bins_x = int((hour_range[1] - hour_range[0]) * bins_per_hour)
    hour_edges = np.linspace(hour_range[0], hour_range[1], num=(n_bins_x + 1))
    delta_hour = hour_edges[1] - hour_edges[0]

    # Compute year edges
    bins_per_year = 24  # Single gate bins
    n_bins_y = int(((year_range[1] + 1) - year_range[0]) * bins_per_year)
    year_edges = np.linspace(year_range[0], year_range[1] + 1, num=(n_bins_y + 1))
    delta_year = year_edges[1] - year_edges[0]

    print("Computing binned occ rates...")
    contour_data_is = np.empty(shape=(n_bins_x, n_bins_y))
    contour_data_is[:] = math.nan

    contour_data_gs = np.empty(shape=(n_bins_x, n_bins_y))
    contour_data_gs[:] = math.nan

    for hour_idx, hour_start in enumerate(hour_edges):
        if hour_start == hour_edges[-1]:
            continue  # The last edge is not a slice start
        hour_end = hour_start + delta_hour
        df_hh = df[(df[time_units] >= hour_start) & (df[time_units] <= hour_end)]

        for year_idx, year_start in enumerate(year_edges):
            if year_start == year_edges[-1]:
                continue  # The last edge is not a slice start
            year_end = year_start + delta_year
            df_hh_yy = df_hh[(df_hh['decimal_year'] >= year_start) & (df_hh['decimal_year'] <= year_end)]

            try:
                contour_data_is[hour_idx][year_idx] = sum(df_hh_yy['good_iono_echo']) / len(df_hh_yy)
                contour_data_gs[hour_idx][year_idx] = sum(df_hh_yy['good_grndscat_echo']) / len(df_hh_yy)
            except ZeroDivisionError:
                # There are no points in this interval
                contour_data_is[hour_idx][year_idx] = math.nan
                contour_data_gs[hour_idx][year_idx] = math.nan
            except BaseException as e:
                print("Hour start: " + str(hour_start))
                print("Year start: " + str(year_start))
                raise e

    # By default, Nan cells will be plotted white.  If this is not what is wanted, replace all nans with 0
    contour_data_is = np.nan_to_num(contour_data_is, nan=0)
    contour_data_gs = np.nan_to_num(contour_data_gs, nan=0)

    # And we can now plot the data
    if plot_type == "contour":

        # Compute bin centers
        bin_xcenters = hour_edges[1:] - delta_hour / 2
        bin_ycenters = year_edges[1:] - delta_year / 2

        levels = 12
        levels = np.linspace(start=0, stop=vmax, num=(levels + 1))
        cmap = 'jet'
        # cmap = modified_jet(levels=len(levels) - 1)

        plot0 = ax[0].contourf(bin_xcenters, bin_ycenters, contour_data_is.transpose(),
                               cmap=cmap, levels=levels, zorder=0)
        plot1 = ax[1].contourf(bin_xcenters, bin_ycenters, contour_data_gs.transpose(),
                               cmap=cmap, levels=levels, zorder=0)

    elif plot_type == "pixel":

        cmap = 'jet'
        plot0 = ax[0].pcolormesh(hour_edges, year_edges, contour_data_is.transpose(),
                                 cmap=cmap, vmin=0, vmax=vmax, zorder=0)
        plot1 = ax[1].pcolormesh(hour_edges, year_edges, contour_data_gs.transpose(),
                                 cmap=cmap, vmin=0, vmax=vmax, zorder=0)

    else:
        raise Exception("plot_type not recognized")

    cbar0 = fig.colorbar(plot0, ax=ax[0], orientation="vertical", format='%.1f')
    cbar0.ax.tick_params(labelsize=14)

    cbar1 = fig.colorbar(plot1, ax=ax[1], orientation="vertical", format='%.1f')
    cbar1.ax.tick_params(labelsize=14)

    if angle_contours is not None and time_units == 'lt':
        if time_units != 'lt':
            # Note: angle contours assume that the data along x is in local time
            # TODO: Update so angle contours work with any time unit
            warnings.warn("Angle contours only plotted if you choose time_units='lt'", category=Warning)
        else:
            if angle_contours != 'zenith' and angle_contours != 'altitude':
                raise Exception("Error in occ_year_vs_time(): param angle_contour not recognized, "
                                "options are 'zenith' and 'altitude'.")

            print("Drawing Angle Contours...")
            # We are good to go ahead and draw angle contours
            for i in range(ax.size):
                add_angle_contours(ax=ax[i], type_of_contours=angle_contours, year_range=year_range,
                                   hour_range=hour_range, beam_range=beam_range, gate_range=gate_range)

    # For some reason gridlines have to be turned on at the end otherwise they don't show up?
    for i in range(ax.size):
        ax[i].grid(which='minor', axis='y', linestyle='--', linewidth=0.2, color='white', zorder=4)
        ax[i].grid(which='major', axis='both', linestyle='--', linewidth=0.5, color='white', zorder=4)

    return df, fig


def add_angle_contours(ax, type_of_contours, year_range, hour_range, beam_range, gate_range):
    """

    Superimpose angle contours of the requested type onto the provided axis.

    Angle contours are based on the center of the field of view (the centroid of the sub-fan).  Therefore it is not
     a good idea to use angle contours if you are plotting data for a large gate/beam range.

    :param ax: matplotlib.pyplot.axis:
            The axis to draw on
    :param type_of_contours: str:
            The type of angle contours to superimpose on top of the data.
                'zenith' for solar zenith angle (the angle between the sun’s rays and a vertical plane)
                'altitude' for solar altitude angle (the angle between the sun’s rays and a horizontal plane)
    :param year_range: See occ_year_vs_time() docstring.
    :param hour_range: See occ_year_vs_time() docstring.
    :param beam_range: See occ_year_vs_time() docstring.
    :param gate_range: See occ_year_vs_time() docstring.

    """

    degree_sign = u'\N{DEGREE SIGN}'
    contour_colour = 'white'

    # Compute hour edges
    bins_per_hour = 4  # quarter hour bins
    n_bins_x = int((hour_range[1] - hour_range[0]) * bins_per_hour)
    hour_edges = np.linspace(hour_range[0], hour_range[1], num=(n_bins_x + 1))
    delta_hour = hour_edges[1] - hour_edges[0]

    # Compute year edges
    bins_per_year = 24  # Single gate bins
    n_bins_y = int(((year_range[1] + 1) - year_range[0]) * bins_per_year)
    year_edges = np.linspace(year_range[0], year_range[1] + 1, num=(n_bins_y + 1))
    delta_year = year_edges[1] - year_edges[0]

    # Compute bin centers
    bin_hour_centers = hour_edges[1:] - delta_hour / 2
    bin_year_centers = year_edges[1:] - delta_year / 2
    hour_centers, year_centers = np.meshgrid(bin_hour_centers, bin_year_centers)

    # We need to convert from decimal year to day of the year, n
    days_in_a_year = 365
    percent_of_year, _ = np.modf(year_centers)
    n = percent_of_year * days_in_a_year
    B = (n - 81) * 360 / 364  # (Eq. 2.2 in Solar Energy Engineering)

    # Compute declination (delta)
    # Declination (Equation 2.5 in Solar Energy Engineering)
    declination = 23.45 * np.sin(np.deg2rad(360 / 365 * (284 + n)))

    # Equation of time (Equation 2.1 in Solar Energy Engineering)
    eq_of_time_in_min = 9.87 * np.sin(2 * B) - 7.53 * np.cos(B) - 1.5 * np.sin(B)  # in minutes
    eq_of_time_in_hours = eq_of_time_in_min / 60

    # No latitude correction is required because I have defined local time based on longitude
    # If hour_centers were in local standard time then a longitude correction would be required here
    apparent_solar_times = hour_centers + eq_of_time_in_hours  # We are assuming hour centers are in LT

    # The hour angle, h, of a point on the earth’s surface is defined as the angle through which the earth
    # would turn to bring the meridian of the point directly under the sun.
    hour_angle = (apparent_solar_times - 12) * 15  # 15 deg of longitude in 1 hour

    # Get local latitude (in matrix form like everything else))
    local_lat = np.empty(shape=hour_centers.shape)
    cent_lon, cent_lat = get_middle_of_fov(station=station, beam_range=beam_range, gate_range=gate_range, coords='geo')
    local_lat[:] = cent_lat

    # Put all parameters into radians  # Convenience
    local_lat = np.deg2rad(local_lat)
    hour_angle = np.deg2rad(hour_angle)
    declination = np.deg2rad(declination)

    # Finally, solve for the angle of interest (Equation 2.12 in Solar Energy Engineering)
    term1 = np.multiply(np.sin(local_lat), np.sin(declination))
    term2 = np.multiply(np.multiply(np.cos(local_lat), np.cos(declination)), np.cos(hour_angle))

    if type_of_contours == 'zenith':
        angles = np.arccos(term1 + term2)
    else:  # 'altitude'
        angles = np.arcsin(term1 + term2)

    # Put angles in degrees
    angles = np.rad2deg(angles)

    contours = ax.contour(hour_centers, year_centers, angles, colors=contour_colour, zorder=5)

    plt.clabel(contours, inline=True, fontsize=12, colors=contour_colour, fmt='%d' + degree_sign,
               inline_spacing=3, zorder=5)

    # Add in a legend specifying what contours represent
    contour_label = "Solar " + type_of_contours.capitalize() + " Angles"
    contours.collections[0].set_label(contour_label)
    ax.legend(loc=(0, 1.02), fontsize=12, facecolor='red')


if __name__ == '__main__':
    """ Testing """

    local_testing = False

    if local_testing:
        station = "rkn"

        # Note: year, month, and day don't matter for local testing
        df, fig = occ_year_vs_time(station=station, year_range=(2011, 2012), month_range=None, day_range=None,
                                   hour_range=None, gate_range=(0, 74), beam_range=None, freq_range=(11, 13),
                                   time_units='lt', plot_type='contour', angle_contours='zenith',
                                   local_testing=local_testing)

        plt.show()


    else:
        station = "dcn"
        freq_range = (8, 10)
        year_range = (2019, 2021)
        angle_contours = 'zenith'

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        _, fig = occ_year_vs_time(station=station, year_range=year_range, day_range=None,
                                  gate_range=(10, 30), beam_range=(6, 8), freq_range=freq_range,
                                  time_units='lt', plot_type='contour', angle_contours=angle_contours,
                                  local_testing=local_testing)

        out_fig = out_dir + "/occ_yearVtime_" + station + "_" + \
                  str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz - with_" + angle_contours + "_contours"

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)

