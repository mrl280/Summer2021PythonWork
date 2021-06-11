import pathlib
import pydarn

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import radar_fov, SuperDARNRadars
from scipy import stats

from lib.add_mlt_to_df import add_mlt_to_df
from lib.cm.modified_viridis import modified_viridis_2
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.z_min_max_defaults import z_min_max_defaults
from lib.build_datetime_epoch import build_datetime_epoch
from lib.data_getters.input_checkers import *


def occ_year_vs_ut(station, year, month, hour_range=None, time_units='mlt',
                   gate_range=None, beam_range=None, local_testing=False):
    """

    Produce a contour plot with gate on the y-axis and time along the x-axis.
    Plots echo occurrence rate

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Only considers 45 km data
        - Does not distinguish frequency
        - This program uses fitACF 3.0 data.  To change this, modify the source code.

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year: int:
            The year to consider
    :param month: int:
            The month to consider
    :param time_units: str: 'ut' for universal time or 'mlt' for magnetic local time:
            The time units to plot along x.  Default is 'mlt'
    :param hour_range: (<int>, <int>) (optional):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not quite inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param gate_range: (<int>, <int>) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :return: matplotlib.pyplot.figure
            The figure can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    time_units = check_time_units(time_units)
    year = check_year(year)
    hour_range = check_hour_range(hour_range)

    print("Retrieving data...")
    df = get_data_handler(station, year_range=(year, year), month_range=(month, month), hour_range=hour_range,
                          gate_range=gate_range, beam_range=beam_range, occ_data=True,
                          local_testing=local_testing)

    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid

    gate_range = check_gate_range(gate_range, this_radars_info.hardware_info)

    print("Preparing the plot...")
    fig, ax = plt.subplots(dpi=300)

    ax.set_ylim(gate_range)
    ax.set_xlim(hour_range)
    print(hour_range)

    ax.xaxis.set_major_locator(MultipleLocator(4))
    ax.tick_params(axis='y', which='major', direction='in')
    ax.tick_params(axis='x', which='major', direction='in')
    ax.grid(b=True, which='both', axis='x', linestyle='--', linewidth=0.5, zorder=4)
    ax.set_ylabel("Range Gate")
    ax.set_xlabel("Time, " + time_units.upper())


    n_bins_y = (gate_range[1] - gate_range[0])  # Single gate bins
    n_bins_x = (hour_range[1] - hour_range[0]) * 2  # half hour bins
    contour_range = [ax.get_xlim(), ax.get_ylim()]

    # TODO: Add title to plot
    # TODO: Bin occurrence data and add to plot
    print(df.head())

    return fig


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    station = "rkn"
    fig = occ_year_vs_ut(station=station, year=2014, month=2, time_units='ut',
                         gate_range=(0, 74), beam_range=None, local_testing=local_testing)

    if local_testing:
        a = 1
        plt.show()
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        out_file = out_dir + "/occ_year_vs_time_" + station
        print("Saving plot as " + out_file)
        fig.savefig(out_file + ".jpg", format='jpg', dpi=300)
