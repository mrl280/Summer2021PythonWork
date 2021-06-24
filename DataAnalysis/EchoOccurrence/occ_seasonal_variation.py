import math
import pathlib
import pydarn

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import radar_fov, SuperDARNRadars
import matplotlib.ticker as mticker

from lib.add_mlt_to_df import add_mlt_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def occ_seasonal_variation(station, year_range=None, day_range=None, hour_range=None,
                           gate_range=None, beam_range=None, freq_range=None,
                           time_units='mlt', local_testing=False):
    """

    Produces an occurrence rate versus year plot that is meant to showcase the seasonal variation in echo occurrence.

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only considers 45 km data
        - To check which fitACF program is being used, refer to the data readers in lib.data_getters

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year_range: (int, int):
            Inclusive. The year range to consider. If omitted (or None), then all years will be considered.
    :param day_range: (<int>, <int>) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param hour_range: (<int>, <int>) (optional):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not quite inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param gate_range: (<int>, <int>) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param freq_range: (<float>, <float>) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.
    :param time_units: str: 'ut' for universal time or 'mlt' for magnetic local time:
            The provided hour_range is assumed to be in these time units.
    :param local_testing: bool (optional): default is False.
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.

    :return: pandas.DataFrame, matplotlib.pyplot.figure
            The dataframe used and the figure created.
            The figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    time_units = check_time_units(time_units)
    year_range = check_year_range(year_range)
    hour_range = check_hour_range(hour_range)
    freq_range = check_freq_range(freq_range)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, this_radars_info.hardware_info)

    if year_range[1] == year_range[0]:
        year_string = str(year_range[0])
    else:
        year_string = str(year_range[0]) + " to " + str(year_range[1])
    beam_string = "Beams " + str(beam_range[0]) + "-" + str(beam_range[1])
    gate_string = "Gates " + str(gate_range[0]) + "-" + str(gate_range[1])
    freq_string = "Frequencies " + str(freq_range[0]) + "-" + str(freq_range[1]) + " MHz"
    hour_string = "Hours " + str(hour_range[0]) + "-" + str(hour_range[1]) + " " + time_units.upper()

    print("Retrieving data...")
    df = get_data_handler(station, year_range=year_range, month_range=None, day_range=day_range, hour_range=hour_range,
                          gate_range=gate_range, beam_range=beam_range, freq_range=freq_range, occ_data=True,
                          local_testing=local_testing)
    df = only_keep_45km_res_data(df)

    print("Preparing the figure...")
    y_lim = [0, 0.8]
    y_axis_major_labels = [0.0, 0.2, 0.4, 0.6, 0.8]
    x_lim = [year_range[0], year_range[1] + 1]  # Make sure we go the end of the final year

    fig, ax = plt.subplots(figsize=[8, 6], dpi=300, constrained_layout=True, nrows=1, ncols=1)

    ax.set_xlim(x_lim)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))

    ax.set_ylim(y_lim)
    ax.yaxis.set_major_locator(mticker.FixedLocator(y_axis_major_labels))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5, color='black')
    ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2, color='black')

    ax.set_xlabel("Year", fontsize=14)
    ax.set_ylabel("Echo Occurrence Rate", fontsize=14)

    fig.suptitle(year_string + " at " + station.upper() +
                 "\n" + gate_string + "; " + beam_string +
                 "\n" + freq_string + "; " + hour_string, fontsize=18)

    print("Computing UT decimal years/hours...")
    minutes_in_an_hour = 60
    seconds_in_an_hour = 3600
    months_in_a_year = 12
    days_in_a_year = 365
    hours_in_a_year = 8760

    decimal_hours, decimal_years = [], []
    for i in range(len(df)):
        datetime_obj_here = df['datetime'].iat[i]

        decimal_hour_here = datetime_obj_here.hour + datetime_obj_here.minute / minutes_in_an_hour + \
                            datetime_obj_here.second / seconds_in_an_hour

        decimal_year_here = datetime_obj_here.year + (datetime_obj_here.month - 1) / months_in_a_year + \
                            (datetime_obj_here.day - 1) / days_in_a_year + decimal_hour_here / hours_in_a_year

        decimal_hours.append(decimal_hour_here)
        decimal_years.append(decimal_year_here)

    df['decimal_hour'] = np.asarray(decimal_hours)
    df['decimal_year'] = np.asarray(decimal_years)

    if hour_range != (0, 24):
        # Then we need to  restrict ourselves to within the desired hour_range
        if time_units == "mlt":
            print("Computing MLT for hour restriction...")

            # To compute mlt we need longitudes.. use the middle of the year as magnetic field estimate
            date_time_est, _ = build_datetime_epoch(year, 6, 15, 0)
            cell_corners_aacgm_lats, cell_corners_aacgm_lons = radar_fov(stid=radar_id, coords='aacgm',
                                                                         date=date_time_est)

            df = add_mlt_to_df(cell_corners_aacgm_lons=cell_corners_aacgm_lons,
                               cell_corners_aacgm_lats=cell_corners_aacgm_lats, df=df)

            df = df.loc[(df['mlt'] >= hour_range[0]) & (df['mlt'] <= hour_range[1])]

        else:
            df = df.loc[(df['decimal_hour'] >= hour_range[0]) & (df['decimal_hour'] <= hour_range[1])]

        df.reset_index(drop=True, inplace=True)

    # Compute year_edges
    bins_per_year = 180
    n_bins_x = ((year_range[1] + 1) - year_range[0]) * bins_per_year
    year_slice_edges = np.linspace(year_range[0], (year_range[1] + 1), num=(n_bins_x + 1))
    delta_year_slice = year_slice_edges[1] - year_slice_edges[0]

    # Compute bin centers
    bin_xwidth = delta_year_slice
    bin_xcenters = year_slice_edges[1:] - bin_xwidth / 2

    # Loop through the time bins and compute occurrence rate
    occurrence_data_is = np.empty(shape=n_bins_x)
    occurrence_data_gs = np.empty(shape=n_bins_x)
    for year_idx, year_slice_start in enumerate(year_slice_edges):
        if year_slice_start == year_slice_edges[-1]:
            continue  # The last edge is not a starting edge

        year_slice_end = year_slice_start + delta_year_slice

        df_yy = df[(df['decimal_year'] >= year_slice_start) & (df['decimal_year'] <= year_slice_end)]

        try:
            occurrence_data_is[year_idx] = sum(df_yy['good_iono_echo']) / len(df_yy)
            occurrence_data_gs[year_idx] = sum(df_yy['good_grndscat_echo']) / len(df_yy)
        except ZeroDivisionError:
            # There are no points in this interval
            occurrence_data_is[year_idx] = math.nan
            occurrence_data_gs[year_idx] = math.nan
        except BaseException as e:
            print("Something bad happened while computing binned occurrence rates..")
            print("Year index: " + str(year_idx))
            print("Year slice start: " + str(year_slice_start))
            raise e

    if local_testing:
        # Plot as a point so we can see it even while testing with little data
        ax.plot(bin_xcenters, occurrence_data_is, 'ro', label='IS')
        ax.plot(bin_xcenters, occurrence_data_gs, 'bo', label='GS')
    else:
        # Plot as a line
        ax.plot(bin_xcenters, occurrence_data_is, color="blue", linestyle='-', label='IS')
        ax.plot(bin_xcenters, occurrence_data_gs, color="red", linestyle='-', label='GS')

    ax.legend(loc='upper right')

    return df, fig


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    if local_testing:
        station = "rkn"

        # Note: year, month, and day don't matter for local testing
        df, fig = occ_seasonal_variation(station=station, year_range=(2011, 2012), day_range=(12, 12),
                                         gate_range=(10, 30), beam_range=(6, 8), freq_range=(11, 13),
                                         time_units='ut', local_testing=local_testing)

        plt.show()


    else:
        station = "dcn"
        year_range = (2019, 2021)
        freq_range = (8, 10)

        datetime_now = datetime.datetime.now()
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        df, fig = occ_seasonal_variation(station=station, year_range=year_range, day_range=None,
                                         gate_range=(10, 30), beam_range=(6, 8), freq_range=freq_range,
                                         time_units='ut', local_testing=local_testing)

        out_fig = out_dir + "/occ_seasonalVariation_" + station + \
                  "_" + str(year_range[0]) + "-" + str(year_range[1]) + \
                  "_" + str(freq_range[0]) + "-" + str(freq_range[1]) + "MHz"

        print("Saving plot as " + out_fig)
        fig.savefig(out_fig + ".jpg", format='jpg', dpi=300)
