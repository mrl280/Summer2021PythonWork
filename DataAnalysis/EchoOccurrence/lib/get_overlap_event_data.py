import pathlib
import glob
import os

import pandas as pd

from adjust_velocities_los import adjust_velocities_los
from get_data_handler import get_data_handler
from only_keep_45km_res_data import only_keep_45km_res_data
from only_keep_overlap import only_keep_overlap
from data_getters.input_checkers import *


def get_overlap_event_df1_and_df2(station1, station2, year, use_only_coinciding_beams=True,
                                  station1_ref_beam=None, station2_ref_beam=None,
                                  gate_range1=None, beam_range1=None, gate_range2=None, beam_range2=None,
                                  local_testing=False):
    """

    This program reads in an overlap event csv file, and builds dataframes with the combined data for every entry
     listed in the file.

    :param station1: str:
            The first radar station to consider, as 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param station2: str:
            The second radar station to consider, again as a 3 character string.
    :param year: int:
            The year to consider.

    :param use_only_coinciding_beams: bool (optional; default is True).
            - If True; apply gate and beam range restrictions, but don't adjust velocities
            - If False; still apply gate and beam range restrictions, but also modify velocities (requires reference
            beams be passed in)
    :param station1_ref_beam: int:
            Only referenced if use_only_coinciding_beams == False.
            The beam used to adjust station1 velocity measurements so they are as if they are looking along this beam.
            We need to modify LoS velocities so that line-of-sight velocity measurements from different directions can
            be compared directly - to do this select two beam (one from each station) that point the same direction and
            adjust all velocities so it is as if all velocity measurements are looking along this same line-of-sight.
    :param station2_ref_beam: int:
            Only referenced if use_only_coinciding_beams == False.
            The beam used to adjust station2 velocity measurements so they are as if they are looking along this beam.

    :param gate_range1: (int, int) (optional):
            Inclusive. The gate range for the first station.  If omitted (or None), then all the gates will be considered.
            Note that early gates probably are not seeing F region echoes, so you probably don't want to consider them
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range1: (int, int) (optional):
            Inclusive. The beam range for the first station.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param gate_range2: (int, int) (optional):
            Inclusive. The gate range for the second station.  If omitted (or None), then all the gates will be considered.
            Note that early gates probably are not seeing F region echoes, so you probably don't want to consider them
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range2: (int, int) (optional):
            Inclusive. The beam range for the second station.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.

    :param local_testing: bool (optional; default is False):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    """

    pattern = "%Y-%m-%d %H:%M"

    loc_root = str((pathlib.Path().parent.absolute()))
    if "data_getters" in loc_root:
        # We are calling from low down
        loc_root = str((pathlib.Path().parent.absolute().parent.absolute().parent.absolute()))
    if "lib" in loc_root:
        # We are calling from within lib
        loc_root = str((pathlib.Path().parent.absolute().parent.absolute()))

    in_dir = loc_root + "/data/overlap_events"
    print("Looking for a pickled file in: " + in_dir)

    for in_file in glob.iglob(in_dir + "/*"):
        file_name = str(os.path.basename(in_file))

        if station1 in file_name and station2 in file_name and str(year) in file_name:
            print("Reading overlap events from " + file_name)
            event_df = pd.read_csv(in_file)
            event_df = event_df.loc[event_df['start'].notna()]

            # Use the first event to initialize the dataframe
            starting_datetime = datetime.datetime.strptime(event_df['start'].iat[0], pattern)
            ending_datetime = datetime.datetime.strptime(event_df['end'].iat[0], pattern)
            start_epoch = event_df['start_epoch'].iat[0]
            end_epoch = event_df['end_epoch'].iat[0]

            year_range = (starting_datetime.year, ending_datetime.year)
            year_range = check_year_range(year_range=year_range)
            month_range = (starting_datetime.month, ending_datetime.month)
            month_range = check_month_range(month_range=month_range)
            day_range = (starting_datetime.day, ending_datetime.day)
            day_range = check_day_range(day_range=day_range)

            df1 = get_data_handler(station=station1, year_range=year_range, month_range=month_range,
                                   day_range=day_range, gate_range=gate_range1, beam_range=beam_range1,
                                   local_testing=local_testing, even_odd_days=None)
            df1 = df1.loc[(df1['epoch'] >= start_epoch) & (df1['epoch'] <= end_epoch)]

            df2 = get_data_handler(station=station2, year_range=year_range, month_range=month_range,
                                   day_range=day_range, gate_range=gate_range2, beam_range=beam_range2,
                                   local_testing=local_testing, even_odd_days=None)
            df2 = df2.loc[(df2['epoch'] >= start_epoch) & (df2['epoch'] <= end_epoch)]

            # Loop through and append the rest of the events
            for i in range(len(event_df)):

                if i == 0:
                    continue

                try:

                    starting_datetime = datetime.datetime.strptime(event_df['start'].iat[i], pattern)
                    ending_datetime = datetime.datetime.strptime(event_df['end'].iat[i], pattern)
                    start_epoch = event_df['start_epoch'].iat[i]
                    end_epoch = event_df['end_epoch'].iat[i]
                except BaseException as e:
                    print(event_df['start'].iat[i])
                    raise e

                year_range = (starting_datetime.year, ending_datetime.year)
                year_range = check_year_range(year_range=year_range)
                month_range = (starting_datetime.month, ending_datetime.month)
                month_range = check_month_range(month_range=month_range)
                day_range = (starting_datetime.day, ending_datetime.day)
                day_range = check_day_range(day_range=day_range)

                df1_here = get_data_handler(station=station1, year_range=year_range, month_range=month_range,
                                            day_range=day_range, gate_range=gate_range1, beam_range=beam_range1,
                                            local_testing=local_testing, even_odd_days=None)
                df1_here = df1_here.loc[(df1_here['epoch'] >= start_epoch) & (df1_here['epoch'] <= end_epoch)]

                df2_here = get_data_handler(station=station2, year_range=year_range, month_range=month_range,
                                            day_range=day_range, gate_range=gate_range2, beam_range=beam_range2,
                                            local_testing=local_testing, even_odd_days=None)
                df2_here = df2_here.loc[(df2_here['epoch'] >= start_epoch) & (df2_here['epoch'] <= end_epoch)]

                df1 = pd.concat([df1, df1_here])
                df2 = pd.concat([df2, df2_here])


            print("Cleaning up the dataframes...")
            df1 = df1.drop_duplicates()
            df2 = df2.drop_duplicates()

            if not use_only_coinciding_beams:
                # We want to use the whole overlap range - cosine adjust velocities
                df1 = adjust_velocities_los(station=station1, df=df1, ref_beam=station1_ref_beam)
                df2 = adjust_velocities_los(station=station2, df=df2, ref_beam=station2_ref_beam)

            df1 = only_keep_45km_res_data(df1)
            df2 = only_keep_45km_res_data(df2)

            df1 = df1.loc[(df1['v'] > -1000) & (df1['v'] < 1000)]  # Remove extreme values
            df2 = df2.loc[(df2['v'] > -1000) & (df2['v'] < 1000)]  # Remove extreme values

            df1.reset_index(drop=True, inplace=True)
            df2.reset_index(drop=True, inplace=True)

            print("Restricting " + station1.upper() + "'s data - only keeping data for those cells that overlap with "
                  + station2.upper())
            df1 = only_keep_overlap(station=station1, df=df1, gate_range=gate_range1, beam_range=beam_range1,
                                    other_station=station2, other_gate_range=gate_range2, other_beam_range=beam_range2)

            print("Restricting " + station2.upper() + "'s data - only keep data for those cells that overlap with "
                  + station1.upper())
            df2 = only_keep_overlap(station=station2, df=df2, gate_range=gate_range2, beam_range=beam_range2,
                                    other_station=station1, other_gate_range=gate_range1, other_beam_range=beam_range1)

            return df1, df2

    raise Exception("get_overlap_event_df1_and_df2(): No event overlap file found for " + station1 +
                    " and " + station2 + " for " + str(year))


if __name__ == "__main__":
    """ Testing """

    local_testing = True

    station1 = "dcn"
    gate_range1 = (20, 74)
    beam_range1 = (14, 15)
    station1_ref_beam = 15

    station2 = "mcm"
    gate_range2 = (20, 74)
    beam_range2 = (8, 8)
    station2_ref_beam = 8

    year = 2019

    df1, df2 = get_overlap_event_df1_and_df2(station1=station1, year=year, gate_range1=gate_range1,
                                             beam_range1=beam_range1, station1_ref_beam=station1_ref_beam,
                                             station2=station2, gate_range2=gate_range2, beam_range2=beam_range2,
                                             local_testing=local_testing)
