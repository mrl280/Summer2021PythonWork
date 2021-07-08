import os

import numpy as np
import pydarn
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import SuperDARNRadars, radar_fov

from DataAnalysis.EchoOccurrence.lib.add_decimal_hour_to_df import add_decimal_hour_to_df
from DataAnalysis.EchoOccurrence.lib.add_mlt_to_df import centroid
from DataAnalysis.EchoOccurrence.lib.build_datetime_epoch import build_datetime_epoch
from DataAnalysis.EchoOccurrence.lib.get_data_handler import get_data_handler
from DataAnalysis.EchoOccurrence.lib.only_keep_45km_res_data import only_keep_45km_res_data


def test_solar_zenith_angle(station, year_range):
    """
    Test computing solar zenith angle and plotting contours

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

    :param station:
    :param year_range:
    :return:
    """

    hour_range = (0, 24)
    time_units = "lt"

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid
    radar_lon = this_radars_info.hardware_info.geographic.lon
    radar_lat = this_radars_info.hardware_info.geographic.lat

    print("Getting some data...")
    df = get_data_handler(station, year_range=year_range, month_range=None, day_range=None,
                          gate_range=(0, 74), beam_range=(6, 8), freq_range=None,
                          occ_data=True, local_testing=True)
    df = df.head(n=200)  # To speed things up
    df = only_keep_45km_res_data(df)

    # Add in longitude based local time
    mid_year = int(year_range[0] + (year_range[1] - year_range[0]) / 2)
    date_time_est, _ = build_datetime_epoch(year=mid_year, month=6, day=15, hour=0)
    df = add_decimal_hour_to_df(df=df, time_units=time_units, stid=radar_id, date_time_est=date_time_est)

    # Compute n, the number of days in the year
    avg_days_in_a_month = 30.42
    n = []
    for i in range(len(df)):
        datetime_obj = df['datetime'].iat[i]
        # N is the day of the year beginning with N=0 at midnight Universal Time (UT) as January 1
        n.append(datetime_obj.month * avg_days_in_a_month + datetime_obj.day)
    df['n'] = np.asarray(n)

    df = add_in_solar_time(df=df)
    df = add_in_local_latitudes(df=df, radar_id=radar_id)
    df = add_in_declination(df=df)


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

        ax[i].grid(which='minor', axis='y', linestyle='--', linewidth=0.2, color='white', zorder=4)
        ax[i].grid(which='major', axis='both', linestyle='--', linewidth=0.5, color='white', zorder=4)


    return fig


def add_in_solar_time(df):
    """

    Add in the solar time, calculation based on
    https://faculty.eng.ufl.edu/jonathan-scheffe/wp-content/uploads/sites/100/2020/08/Solar-Time1419.pdf

    :param df: The dataframe, with column 'n'
    """

    df['B'] = (df['n'] - 1) * 360 / 365

    df['E'] = 299.2 * (0.000075 * 0.001868 * np.cos(df['B']) - 0.032022 * np.sin(df['B'])
                       - 0.014615 * np.cos(2 * df['B']) - 0.04089 * np.sin(2 * df['B']))

    df['solar_time'] = df['lt'] + df['E']

    return df


def add_in_local_latitudes(df, radar_id):
    """

    Add 'local_lat' to dataframe (latitudes in geographical coordinates)

    :param df: The dataframe
    """

    print("Getting local latitudes...")
    mid_year = int(year_range[0] + (year_range[1] - year_range[0]) / 2)
    date_time_est, _ = build_datetime_epoch(year=mid_year, month=6, day=15, hour=0)
    cell_corners_lats, cell_corners_lons = radar_fov(stid=radar_id, coords='geo', date=date_time_est)

    # Compute cell centroids
    fan_shape = cell_corners_lons.shape
    cell_centers_aacgm_lons = np.zeros(shape=(fan_shape[0], fan_shape[1]))
    cell_centers_aacgm_lats = np.zeros(shape=(fan_shape[0], fan_shape[1]))

    for gate_corner in range(fan_shape[0] - 1):
        for beam_corner in range(fan_shape[1] - 1):
            cent_lon, cent_lat = centroid([(cell_corners_lons[gate_corner, beam_corner],
                                            cell_corners_lats[gate_corner, beam_corner]),
                                           (cell_corners_lons[gate_corner + 1, beam_corner],
                                            cell_corners_lats[gate_corner + 1, beam_corner]),
                                           (cell_corners_lons[gate_corner, beam_corner + 1],
                                            cell_corners_lats[gate_corner, beam_corner + 1]),
                                           (cell_corners_lons[gate_corner + 1, beam_corner + 1],
                                            cell_corners_lats[gate_corner + 1, beam_corner + 1])])
            cell_centers_aacgm_lons[gate_corner, beam_corner] = cent_lon
            cell_centers_aacgm_lats[gate_corner, beam_corner] = cent_lat

    #  Loop through the dataframe, and build up local latitudes
    local_lat = []
    for i in range(len(df)):
        gate = df['slist'].iat[i]
        beam = df['bmnum'].iat[i]

        local_lat.append(cell_centers_aacgm_lons[gate, beam])

    df['local_lat'] = local_lat

    return df


def add_in_declination(df):
    """

    Add 'declination' to dataframe

    """

    print("Computing declination of the sun...")

    avg_days_in_a_month = 30.42

    # Compute declination of the sun
    # https://en.wikipedia.org/wiki/Position_of_the_Sun#Declination_of_the_Sun_as_seen_from_Earth
    N = []
    for i in range(len(df)):
        datetime_obj = df['datetime'].iat[i]
        # N is the day of the year beginning with N=0 at midnight Universal Time (UT) as January 1
        N.append(datetime_obj.month * avg_days_in_a_month + datetime_obj.day - 1)
    df['N'] = np.asarray(N)

    # TODO: Compare to formula on page 62 of the textbook
    df['declination'] = -23.44 * np.cos(np.deg2rad(360 / 365 * (df['N'] + 10)))

    return df


if __name__ == "__main__":
    """ Testing """

    station = "rkn"
    year_range = (2011, 2012)

    fig = test_solar_zenith_angle(station=station, year_range=year_range)

    plt.show()
    plt.close(fig)
