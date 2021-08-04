import pydarn
from pydarn import SuperDARNRadars

from lib.get_data_handler import get_data_handler
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.only_keep_overlap import only_keep_overlap
from lib.data_getters.input_checkers import *


def find_high_vel_overlap_events(station1, station2, year, month_range=None, day_range=None, freq_range=None,
                                 event_duration_h=4, local_testing=False, even_odd_days=None):
    """

    Given two SuperDARN radars, we want to find events where both radars are reporting high velocity echoes in the
    overlapping cells.  These events will then be further analysed to determine if the two radars are measuring similar
    heights/velocities, or if each radar is in it's own world.

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

    # Early gates probably are not seeing F region echoes, so don't consider them
    gate_min = 30

    year = check_year(year=year)
    freq_range = check_freq_range(freq_range=freq_range)
    day_range = check_day_range(day_range=day_range)

    all_radars_info = SuperDARNRadars()
    first_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station1).stid]  # Grab radar info
    second_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station1).stid]  # Grab radar info

    print("Getting SuperDARN data for " + station1.upper())
    df1 = get_data_handler(station1, year_range=(year, year), month_range=None, day_range=day_range,
                           gate_range=None, beam_range=None, freq_range=freq_range,
                           local_testing=local_testing, even_odd_days=even_odd_days)
    df1 = only_keep_45km_res_data(df1)

    print("Restricting " + station1.upper() + "'s data - only keep data for those cells that overlap with "
          + station2.upper())
    df1 = only_keep_overlap(station=station1, df=df1, other_station=station2, gate_min=gate_min)

    print("Getting SuperDARN data for " + station2.upper())
    df2 = get_data_handler(station2, year_range=(year, year), month_range=None, day_range=day_range,
                           gate_range=None, beam_range=None, freq_range=freq_range,
                           local_testing=local_testing, even_odd_days=even_odd_days)
    df2 = only_keep_45km_res_data(df2)

    print("Restricting " + station2.upper() + "'s data - only keep data for those cells that overlap with "
          + station1.upper())
    df2 = only_keep_overlap(station=station2, df=df2, other_station=station1, gate_min=gate_min)

    # Both dataframes are now restricted to only areas of overlap

    # TODO: Loop though the event in intervals of event_duration_h, and look at meadian or mean |velocity|.
    #  We want events where both radars see high velocities.



