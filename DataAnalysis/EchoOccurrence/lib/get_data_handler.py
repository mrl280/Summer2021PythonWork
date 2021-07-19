import glob
import os
import pathlib
import pydarn

import pandas as pd

try:
    from .data_getters.get_data import get_data
    from .data_getters.get_data_occ import get_data_occ
    from .data_getters.get_local_dummy_data import get_local_dummy_data
    from .data_getters.input_checkers import *
except ImportError:
    # Needed for testing
    from data_getters.get_data import get_data
    from data_getters.get_data_occ import get_data_occ
    from data_getters.get_local_dummy_data import get_local_dummy_data
    from data_getters.input_checkers import *


def get_data_handler(station, year_range=None, month_range=None, day_range=None,
                     gate_range=None, beam_range=None, freq_range=None, occ_data=False,
                     local_testing=False, force_build=False, fitACF_version=2.5):
    """

    Get the required data.

    This might be reading in a pickled file, or it might mean building a dataframe.
     Either way, a dataframe with the requested data is returned for plotting/analysis


    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year_range: (<int>, <int>):
            Inclusive. The year range to consider.
    :param month_range: (<int>, <int>) (optional):
            Inclusive. The months of the year to consider.  If omitted (or None), then all days will be considered.
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
    :param occ_data: bool (optional):
            Set this to true if you need echo occurrence data.
            If False, you will get normal data (basically a reduced fitACF datafile),  Default if False.
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :param force_build: bool (optional; default is False).
            Force a dataframe build, regardless of if we find a suitable pickle file or not, we are going to build a
            dataframe from fitACF files
    :param fitACF_version: float: (optional):
            The fitACF version number.  At the time of writing the default is 2.5, but expected to move to 3.0 in the
            near future.  These are the only valid options.


    :return: pandas.DataFrame: A dataframe with select fitACF parameters.
    """

    if isinstance(station, str):
        hdw_info = pydarn.read_hdw_file(station)  # Get the hardware file
    else:
        raise Exception("Error: Please enter the station as a 3 character string (e.g. 'rkn')")

    gate_range = check_gate_range(gate_range, hdw_info)
    beam_range = check_beam_range(beam_range, hdw_info)
    freq_range = check_freq_range(freq_range)

    if local_testing:
        # Just read in some test data

        if station == "dce":
            df = get_local_dummy_data(station=station, year=2019, month=3, day=3, start_hour_UT=0, end_hour_UT=24,
                                      occ_data=occ_data)
        else:
            # Note: There is no fitACF version check with local data, you get what is there

            warnings.warn("Running in local testing mode, we are just going to use local dummy data from"
                          " Nov 12, 2011 at RKN. Gate and Beam filtering is still applied.", category=Warning)

            # df = get_local_dummy_data(station=station, year=2011, month=9, day=29, start_hour_UT=0, end_hour_UT=23)
            df = get_local_dummy_data(station=station, year=2011, month=11, day=12, start_hour_UT=0, end_hour_UT=24,
                                      occ_data=occ_data)
            # df_2 = get_local_dummy_data(station=station, year=2011, month=9, day=29, start_hour_UT=0, end_hour_UT=24,
            #                             occ_data=occ_data)
            # df_3 = get_local_dummy_data(station=station, year=2012, month=11, day=12, start_hour_UT=0, end_hour_UT=24,
            #                             occ_data=occ_data)
            # df = pd.concat([df, df_2, df_3])

    else:
        # Assume are on maxwell.usask.ca, check the time parameters and use get_data()
        year_range = check_year_range(year_range)
        month_range = check_month_range(month_range)
        day_range = check_day_range(day_range)

        if force_build:
            # Then we don't need to look for a pickled file, just go ahead and build the dataframe from scratch
            warnings.warn("We are building the dataframe ourself, this will take a while...", category=Warning)
            if occ_data:
                df = get_data_occ(station=station, year_range=year_range, month_range=month_range, day_range=day_range,
                                  gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                                  fitACF_version=fitACF_version)
            else:
                df = get_data(station=station, year_range=year_range, month_range=month_range, day_range=day_range,
                              gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                              fitACF_version=fitACF_version)

        else:
            # Try to just read in data from file
            df = look_for_pickle(station=station, year_range=year_range, occ_data=occ_data, fitACF_version=fitACF_version)

            if df is None:
                # Then we have to go ahead and build the dataframe ourselves
                warnings.warn("No pickle found, we are building the dataframe ourself, this will take a while...",
                              category=Warning)
                if occ_data:
                    df = get_data_occ(station=station, year_range=year_range, month_range=month_range,
                                      day_range=day_range, gate_range=gate_range, beam_range=beam_range,
                                      freq_range=freq_range, fitACF_version=fitACF_version)
                else:
                    df = get_data(station=station, year_range=year_range, month_range=month_range, day_range=day_range,
                                  gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                                  fitACF_version=fitACF_version)

            else:
                # We have a pickled data frame, restrict it to the day, month, and year ranges of interest
                year, month, day = [], [], []
                for i in range(len(df)):
                    datetime_obj = df['datetime'].iat[i]
                    year.append(datetime_obj.year)
                    month.append(datetime_obj.month)
                    day.append(datetime_obj.day)

                df = df.loc[((year >= year_range[0]) & (year <= year_range[1])) &
                            ((month >= month_range[0]) & (month <= month_range[1]))
                            ((day >= day_range[0]) & (day <= day_range[1]))]

    # Always make sure we are providing data within the desired range, no matter where it is coming from
    df = df.loc[(df['bmnum'] >= beam_range[0]) & (df['bmnum'] <= beam_range[1]) &
                (df['slist'] >= gate_range[0]) & (df['slist'] <= gate_range[1]) &

                # Note: freq_range is in MHz while data in 'tfreq' is in kHz
                (df['tfreq'] >= freq_range[0] * 1000) & (df['tfreq'] <= freq_range[1] * 1000)]
    df.reset_index(drop=True, inplace=True)

    return df


def look_for_pickle(station, year_range, occ_data, fitACF_version):
    """

    Try to find a pickled datafile that will do the trick.

    Unfortunately, pickle files would take up too much space if they are for all months, days, gates, and beams.
    Therefore, the pickle file might not have everything you are looking for


    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year_range: (<int>, <int>):
            Inclusive. The year range to consider.
    :param occ_data: bool:
            Set this to true if you need echo occurrence data.
            If False, you will get normal data (basically a reduced fitACF datafile)
    :param fitACF_version: float: (optional):
            The fitACF version number.  At the time of writing the default is 2.5, but expected to move to 3.0 in the
            near future.  These are the only valid options.

    :return: A pandas.DataFrame object, possibly None:
            If we find a valid pickled datafile then go ahead and return it, otherwise return None
    """

    if fitACF_version == 3.0:
        needed_fitACT_string = "fitacf_30"
    elif fitACF_version == 2.5:
        needed_fitACT_string = "fitacf_25"
    else:
        raise Exception("look_for_pickle(): fitACF_version " + str(fitACF_version) +
                        " not recognized.  Right now only versions 2.5 and 3.0 are recognized.")

    # We need to look for a pickled file before we do anything rash
    loc_root = str((pathlib.Path().parent.absolute()))
    if "lib" in loc_root:
        # We are calling from within lib
        loc_root = str((pathlib.Path().parent.absolute().parent.absolute()))

    in_dir = loc_root + "/data/" + station
    print("Looking for a pickled file in: " + in_dir)

    for in_file in glob.iglob(in_dir + "/*"):
        file_name = str(os.path.basename(in_file))

        try:
            # Try to pull the required information our of the filename
            file_station = file_name[0:3]
            file_year_start = int(file_name[4:8])
            file_year_end = int(file_name[9:13])
            file_fitACF_string = file_name[14:23]

            # Lets see if this file will do the trick
            if file_station == station and file_fitACF_string == needed_fitACT_string and \
                    file_year_start <= year_range[0] and file_year_end >= year_range[1]:
                # We are probably good, just need to confirm occ

                if file_name[24:27] == 'occ':
                    # It is an occ_file
                    if occ_data:
                        warnings.warn("Using data from " + in_file + ".  This file might not have data for all months, "
                                                                     "days, gates, and beams", category=Warning)
                        return pd.read_pickle(in_file)
                    else:
                        continue

                else:
                    # It is not an occ file
                    warnings.warn("Using data from " + in_file + ".  This file might not have data for all months, "
                                                                 "days, gates, and beams", category=Warning)
                    return pd.read_pickle(in_file)

            else:
                continue  # This file is not what we are looking for
        except BaseException as e:
            # print(e)
            pass

    # If we didn't find anything, then we can just go ahead and return None
    return None


if __name__ == "__main__":
    """ Testing """

    station = "dce"
    year_range = (2019, 2021)
    occ_data = False

    df = get_data_handler(station, year_range=year_range, month_range=None, day_range=None,
                          gate_range=None, beam_range=None, freq_range=None, occ_data=occ_data,
                          local_testing=False, fitACF_version=2.5)

    print(df.head())
