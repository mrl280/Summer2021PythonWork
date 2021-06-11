import pandas as pd
import pydarn

from .data_getters.get_data import get_data
from .data_getters.get_data_occ import get_data_occ
from .data_getters.get_local_dummy_data import get_local_dummy_data
from .data_getters.input_checkers import *


def get_data_handler(station, year_range=None, month_range=None, day_range=None, hour_range=None,
                     gate_range=None, beam_range=None, occ_data=False, local_testing=False):
    """

    Get the required data, put it into a dataframe, and then return the dataframe for plotting/analysis

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year_range: (<int>, <int>):
            Inclusive. The year range to consider.
    :param month_range: (<int>, <int>) (optional):
            Inclusive. The months of the year to consider.  If omitted (or None), then all days will be considered.
    :param day_range: (<int>, <int>) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param hour_range: (<int>, <int>) (optional):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param gate_range: (<int>, <int>) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param occ_data: bool (optional):
            Set this to true if you need echo occurrence data.
            If False, you will get normal data (basically a reduced fitACF datafile),  Default if False.
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :return: pandas.DataFrame: A dataframe with select fitACF parameters.
    """

    if isinstance(station, str):
        hdw_info = pydarn.read_hdw_file(station)  # Get the hardware file
    else:
        raise Exception("Error: Please enter the station as a 3 character string (e.g. 'rkn')")

    gate_range = check_gate_range(gate_range, hdw_info)
    beam_range = check_beam_range(beam_range, hdw_info)

    if local_testing:
        # Just read in some test data
        warnings.warn("Running in local testing mode, we are just going to use local dummy data from"
                      " Nov 12, 2011, Sept 29, 2011, and Nov 10, 2015 all at RKN."
                      " Gate and Beam filtering is still applied.", category=Warning)

        # df = get_local_dummy_data(station=station, year=2011, month=9, day=29, start_hour_UT=0, end_hour_UT=23)
        df = get_local_dummy_data(station=station, year=2011, month=11, day=12, start_hour_UT=0, end_hour_UT=24,
                                  occ_data=occ_data)
        df_2 = get_local_dummy_data(station=station, year=2011, month=9, day=29, start_hour_UT=0, end_hour_UT=24,
                                    occ_data=occ_data)
        df_3 = get_local_dummy_data(station=station, year=2012, month=10, day=15, start_hour_UT=0, end_hour_UT=24,
                                    occ_data=occ_data)
        df = pd.concat([df, df_2, df_3])
        df = df.loc[(df['bmnum'] >= beam_range[0]) & (df['bmnum'] <= beam_range[1]) &
                    (df['slist'] >= gate_range[0]) & (df['slist'] <= gate_range[1])]

    else:
        # Assume are on maxwell.usask.ca, check the time parameters and use get_data()
        year_range = check_year_range(year_range)
        month_range = check_month_range(month_range)
        day_range = check_day_range(day_range)
        hour_range = check_hour_range(hour_range)

        if occ_data:
            df = get_data_occ(station=station, year_range=year_range, month_range=month_range, day_range=day_range,
                              gate_range=gate_range, beam_range=beam_range)
        else:
            df = get_data(station=station, year_range=year_range, month_range=month_range, day_range=day_range,
                          hour_range=hour_range, gate_range=gate_range, beam_range=beam_range)

    return df
