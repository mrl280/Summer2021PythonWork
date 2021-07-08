import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

from DataAnalysis.EchoOccurrence.lib.get_middle_of_fov import get_middle_of_fov


def test_solar_zenith_angle():
    """
    Test computing solar zenith angle and plotting contours
    Computations are based on Solar Energy Engineering. Processes and Systems by Soteris A. Kalogirou

    Important links:

    Solar zenith angle:
    https://en.wikipedia.org/wiki/Solar_zenith_angle

    Declination:
    https://en.wikipedia.org/wiki/Position_of_the_Sun#Declination_of_the_Sun_as_seen_from_Earth

    Solar time:
    https://susdesign.com/popups/sunangle/time-basis.php
    https://faculty.eng.ufl.edu/jonathan-scheffe/wp-content/uploads/sites/100/2020/08/Solar-Time1419.pdf

    Other links:
    https://www.sciencedirect.com/topics/engineering/solar-declination


    """

    # station = "inv"
    # beam_range = (13, 15)

    station = "dce"
    beam_range = (10, 12)

    year_range = (2013, 2014)
    hour_range = (0, 24)
    gate_range = (10, 30)
    time_units = "lt"

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

    # Compute n, the number of days in the year
    # We need to convert from decimal year to day of the year
    days_in_a_year = 365
    percent_of_year, _ = np.modf(year_centers)
    n = percent_of_year * days_in_a_year

    # Compute declination (delta)
    # Declination (Equation 2.5 in Solar Energy Engineering)
    declination = 23.45 * np.sin(np.deg2rad(360 / 365 * (284 + n)))

    B = (n - 81) * 360 / 364  # (Eq. 2.2 in Solar Energy Engineering)
    # Compute equation of time (Equation 2.1 in Solar Energy Engineering)
    ET_in_min = 9.87 * np.sin(2 * B) - 7.53 * np.cos(B) - 1.5 * np.sin(B)  # in minutes
    ET_in_hours = ET_in_min / 60

    AST = hour_centers + ET_in_hours

    # The hour angle, h, of a point on the earth’s surface is defined as the angle through which the earth
    # would turn to bring the meridian of the point directly under the sun.
    hour_angle = (AST - 12) * 15  # 15 deg of longitude in 1 hour

    # Get local latitude
    cent_lon, cent_lat = get_middle_of_fov(station=station, beam_range=beam_range, gate_range=gate_range, coords='geo')

    # Put it into matrix form like everything else
    local_lat = np.empty(shape=hour_centers.shape)
    local_lat[:] = cent_lat

    # Put all parameters into radians
    local_lat = np.deg2rad(local_lat)
    hour_angle = np.deg2rad(hour_angle)
    declination = np.deg2rad(declination)

    # Finally, solve for solar Zenith angle (Equation 2.12 in Solar Energy Engineering)
    term1 = np.multiply(np.sin(local_lat), np.sin(declination))
    term2 = np.multiply(np.multiply(np.cos(local_lat), np.cos(declination)), np.cos(hour_angle))

    solar_zenith_angles = np.arccos(term1 + term2)
    solar_altitude_angles = np.arcsin(term1 + term2)

    # Put the zenith and altitude angles into degrees
    solar_zenith_angles = np.rad2deg(solar_zenith_angles)
    solar_altitude_angles = np.rad2deg(solar_altitude_angles)


    fig, ax = plt.subplots(figsize=[10, 12], dpi=300, nrows=2, ncols=1, constrained_layout=True)
    fig.suptitle(station.upper() + "\nProduced by " + str(os.path.basename(__file__)), fontsize=18)

    ax[0].set_title("Ionospheric Scatter", fontsize=16)
    ax[1].set_title("Ground Scatter", fontsize=16)

    # Apply common subplot formatting
    for i in range(ax.size):
        ax[i].set_xlim(hour_range)
        ax[i].set_xlabel("Time, " + time_units.upper(), fontsize=16)
        ax[i].xaxis.set_major_locator(MultipleLocator(4))

        ax[i].set_ylim((year_range[0], year_range[1] + 1))
        ax[i].set_ylabel("Year", fontsize=16)
        ax[i].yaxis.set_major_locator(MultipleLocator(1))
        ax[i].yaxis.set_minor_locator(MultipleLocator(0.25))

        ax[i].tick_params(axis='both', which='major', direction='in', color='black', labelsize=14)
        ax[i].tick_params(axis='y', which='minor', direction='in', color='black', labelsize=14)

        ax[i].grid(which='minor', axis='y', linestyle='--', linewidth=0.2, color='black', zorder=4)
        ax[i].grid(which='major', axis='both', linestyle='--', linewidth=0.5, color='black', zorder=4)

    # ax[0].scatter(hour_centers, year_centers, c='black', s=1)

    degree_sign = u'\N{DEGREE SIGN}'

    contour_colour = 'blue'
    contour_levels = [60, 75, 90, 105, 120]  # Degrees
    contours = ax[0].contour(hour_centers, year_centers, solar_zenith_angles,
                             colors=contour_colour, levels=contour_levels, zorder=5)
    plt.clabel(contours, inline=True, fontsize=12, colors=contour_colour, fmt='%d' + degree_sign,
               inline_spacing=3, zorder=5)

    return fig


if __name__ == "__main__":
    """ Testing """

    fig = test_solar_zenith_angle()

    plt.show()
    plt.close(fig)
