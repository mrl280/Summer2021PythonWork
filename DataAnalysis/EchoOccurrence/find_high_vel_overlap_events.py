import pathlib
import dill as pickle
import pydarn
import bz2

import numpy as np
import pandas as pd

from pydarn import SuperDARNRadars

from lib.build_datetime_epoch import build_datetime_epoch
from lib.get_data_handler import get_data_handler
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.only_keep_overlap import only_keep_overlap
from lib.data_getters.input_checkers import *


def find_high_vel_overlap_events(station1, station2, year, month_range=None, day_range=None,
                                 gate_range=None, beam_range=None, freq_range=None,
                                 event_duration_h=4, local_testing=False, even_odd_days=None):
    """

    Given two SuperDARN radars, we want to find events where both radars are reporting high velocity echoes in the
    overlapping cells.  These events will then be further analysed to determine if the two radars are measuring similar
    heights/velocities, or if each radar is in its own world.

    Data is grouped into chunks of length event_duration_h, and each chunk is assessed individually

    If there is no overlap between the radars, an exception is raised

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            we want to look thorough large amounts of data - and that data is on Maxwell
        - To check which fitACF program is being used, refer to the data readers in lib.data_getters
        - Only considers 45 km data.


    :param station1: str:
            The first radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param station2: str:
            The second radar station to consider, again as a 3 character string.
    :param year: int:
            The year to consider.
    :param day_range: (int, int) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param gate_range: (int, int) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that early gates probably are not seeing F region echoes, so you probably don't want to consider them
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (int, int) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param freq_range: (float, float) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.

    :param event_duration_h: int (optional; default is None)
            The duration of event to look for.
    :param even_odd_days: str (optional; default is None)
            'even': only even days are read in
            'odd': only odd days are read in
            None: all days are read in
    :param local_testing: bool (optional; default is False)
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.

    :return event_df: pandas.DataFrame:
            A dataframe with all of the events found.  The dataframe can then be printed out or pickled for later.
    """

    year = check_year(year=year)
    freq_range = check_freq_range(freq_range=freq_range)
    day_range = check_day_range(day_range=day_range)
    month_range = check_month_range(month_range=month_range)

    all_radars_info = SuperDARNRadars()
    first_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station1).stid]  # Grab radar info
    second_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station1).stid]  # Grab radar info

    # Check the gate and beam range against both stations hardware info
    gate_range = check_gate_range(gate_range, first_radars_info.hardware_info)
    gate_range = check_gate_range(gate_range, second_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, first_radars_info.hardware_info)
    beam_range = check_beam_range(beam_range, second_radars_info.hardware_info)

    print("Getting SuperDARN data for " + station1.upper())
    df1 = get_data_handler(station1, year_range=(year, year), month_range=month_range, day_range=day_range,
                           gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                           local_testing=local_testing, even_odd_days=even_odd_days)
    df1 = only_keep_45km_res_data(df1)

    print("Restricting " + station1.upper() + "'s data - only keep data for those cells that overlap with "
          + station2.upper())
    df1 = only_keep_overlap(station=station1, df=df1, other_station=station2, gate_min=gate_range[0])

    print("Getting SuperDARN data for " + station2.upper())
    df2 = get_data_handler(station2, year_range=(year, year), month_range=month_range, day_range=day_range,
                           gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                           local_testing=local_testing, even_odd_days=even_odd_days)
    df2 = only_keep_45km_res_data(df2)

    print("Restricting " + station2.upper() + "'s data - only keep data for those cells that overlap with "
          + station1.upper())
    df2 = only_keep_overlap(station=station2, df=df2, other_station=station1, gate_min=gate_range[0])

    # We are interested in high velocity events - sign doesn't matter.  So, remove extreme velocity outliers and make
    #  all velocities positive
    df1['v'] = np.abs(df1['v'])
    df1 = df1.loc[(df1['v'] < 1000)]
    df1.reset_index(drop=True, inplace=True)

    df2['v'] = np.abs(df2['v'])
    df2 = df2.loc[(df2['v'] < 1000)]
    df2.reset_index(drop=True, inplace=True)

    # Both dataframes are now restricted to only areas of overlap, so we can now loop though the event in intervals of
    #  event_duration_h, and look at mean |velocity|.

    # Convert event duration in hours to seconds
    seconds_in_an_hour = 3600
    delta_t = event_duration_h * seconds_in_an_hour

    starting_datetime, starting_epoch = build_datetime_epoch(year=year, month=month_range[0], day=day_range[0],
                                                             hour=0, minute=0, second=0)
    ending_datetime, ending_epoch = build_datetime_epoch(year=year, month=month_range[1], day=day_range[1],
                                                         hour=23, minute=59, second=59)

    # print("Starting datetime: " + str(starting_datetime))
    # print("Starting epoch: " + str(starting_epoch))
    #
    # print("Ending datetime: " + str(ending_datetime))
    # print("Ending epoch: " + str(ending_epoch))

    # We don't want to accidentally hit half of a good event and have it ruined by a partial quite period.
    #  So, loop though in with different starting points
    half_hour_in_s = 0.5 * seconds_in_an_hour
    starting_epochs = np.arange(start=starting_epoch, stop=starting_epoch + delta_t, step=half_hour_in_s)

    # print("Here are the starting epoch offsets:")
    # print(starting_epochs)

    event_starting_epochs, event_ending_epochs = [], []
    station1_mean_event_velocities, station2_mean_event_velocities = [], []

    # Go ahead and run through each of the time chunks for each starting epoch
    for epoch_start in starting_epochs:
        epoch_edges = np.arange(start=epoch_start, stop=ending_epoch, step=delta_t)

        print("Running through event with starting epoch " + str(epoch_start))
        # print("Epoch edges: ")
        # print(epoch_edges)

        for starting_edge_epoch in epoch_edges:
            if starting_edge_epoch == epoch_edges[-1]:
                continue  # The last edge is not a slice start

            ending_edge_epoch = starting_edge_epoch + delta_t
            df_tt_1 = df1[(df1['epoch'] >= starting_edge_epoch) & (df1['epoch'] <= ending_edge_epoch)]
            df_tt_2 = df2[(df2['epoch'] >= starting_edge_epoch) & (df2['epoch'] <= ending_edge_epoch)]

            if len(df_tt_1) <= 100 or len(df_tt_2) <= 100:
                # Then there is no need to note this event
                continue
            else:
                event_starting_epochs.append(starting_edge_epoch)
                event_ending_epochs.append(ending_edge_epoch)

                station1_mean_event_velocities.append(np.mean(df_tt_1['v']))
                station2_mean_event_velocities.append(np.mean(df_tt_2['v']))

    # Put the data into a dataframe
    print("     Putting all the events into a dataframe...")
    event_df = pd.DataFrame({'station1': [station1] * len(event_starting_epochs),
                             'station2': [station2] * len(event_starting_epochs),
                             'start_epoch': event_starting_epochs, 'end_epoch': event_ending_epochs,
                             'station1_mean_vel': station1_mean_event_velocities,
                             'station2_mean_vel': station2_mean_event_velocities
                             })

    # We want high velocity events - so sort the dataframe by the mean velocities
    event_df.sort_values(by=['station1_mean_vel', 'station2_mean_vel'], ascending=[True, True], inplace=True)
    event_df.reset_index(drop=True, inplace=True)

    return event_df


if __name__ == '__main__':
    """ Testing """

    local_testing = True

    if local_testing:
        station1 = "dcn"
        station2 = "mcm"
        year = 2011

        df = find_high_vel_overlap_events(station1=station1, station2=station2, year=year, month_range=None,
                                          day_range=None, gate_range=None, beam_range=None, freq_range=None,
                                          event_duration_h=4, local_testing=local_testing, even_odd_days=None)

        print(df.head())


    else:
        station1 = "dcn"
        station2 = "mcm"
        year = 2019
        event_duration_h = 4
        even_odd_days = None
        freq_range = (8, 11)
        gate_range = (30, 74)  # Close range gates are probably not seeing F-region echoes
        beam_range = (0, 15)

        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"

        df = find_high_vel_overlap_events(station1=station1, station2=station2, year=year, day_range=None,
                                          gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                                          event_duration_h=event_duration_h, even_odd_days=even_odd_days,
                                          local_testing=local_testing)

        out_file = out_dir + "/list_of_overlap_events-" + station1 + "_and_" + station2 + ".pbz2"

        print("     Pickling as " + out_file + "...")
        with bz2.BZ2File(out_file, "w") as file:
            pickle.dump(df, file, protocol=4)
