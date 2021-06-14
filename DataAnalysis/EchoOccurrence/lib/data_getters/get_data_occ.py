import bz2
import glob
import os
import calendar
import time
import pydarn

import datetime as datetime
import pandas as pd
import numpy as np


def get_data_occ(station, year_range, month_range, day_range, gate_range, beam_range):
    """
    Take a fitACF file and for each possible echo, record whether or not there was a good echo there
    Note: occurrence rate is not actually computed here, but all the data required to compute it is put into the df

    From Koustov and Syd's paper:
    "The echo occurrence rate was computed as a ratio of the number of registered echoes in selected beams and gates
    to the total number of possible echo detections in the same gates over the same period."

    :param station: str:
            The radar station to consider, as a 3 character string (e.g. "rkn").
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
    :return: pandas.DataFrame: A dataframe with select fitACF parameters.
    """

    loc_root = "/data/fitacf_30"  # fitACF 3.0
    # loc_root = "/data/fitacf_25"  # fitACF 2.5

    # Create empty arrays for the parameters we need
    epoch, date_time = [], []
    slist, beam = [], []
    frang, rsep, tfreq = [], [], []
    good_echo = []

    # Data on maxwell.usask.ca is stored in a file structure that looks like this:
    # /data/fitacf_30/<year>/<month>/<year><month><day>.<start-time>.<radar-site>.fitacf.bz2
    for in_dir_parent in glob.iglob(loc_root + "/*"):
        year_here = int(os.path.basename(in_dir_parent))

        if year_range[0] <= year_here <= year_range[1]:
            for in_dir in glob.iglob(in_dir_parent + "/*"):
                month_here = int(os.path.basename(in_dir))

                if month_range[0] <= month_here <= month_range[1]:
                    for in_file in glob.iglob(in_dir + "/*" + station + "*"):
                        day_here = int(os.path.basename(in_file)[6:8])

                        if day_range[0] <= day_here <= day_range[1]:
                            # We will read in the whole day, and worry about hour restrictions later
                            print("    Reading: " + str(in_file))
                            try:
                                with bz2.open(in_file) as fp:
                                    fitacf_stream = fp.read()
                                sdarn_read = pydarn.SuperDARNRead(fitacf_stream, True)
                                fitacf_data = sdarn_read.read_fitacf()  # this is
                            except BaseException:
                                # Sometimes files are corrupted, or there is something wrong with them
                                pass

                            # Loop through every record in this file
                            for record in range(len(fitacf_data)):
                                beam_here = fitacf_data[record]['bmnum']
                                if beam_here < beam_range[0] or beam_here > beam_range[1]:
                                    continue  # We are not interested in records for this beam

                                date_time_here, epoch_here = build_datetime_epoch(
                                    year=fitacf_data[record]['time.yr'],
                                    month=fitacf_data[record]['time.mo'],
                                    day=fitacf_data[record]['time.dy'],
                                    hour=fitacf_data[record]['time.hr'],
                                    minute=fitacf_data[record]['time.mt'],
                                    second=fitacf_data[record]['time.sc'])

                                try:
                                    gates_reporting = np.ndarray.tolist(fitacf_data[record]['slist'])
                                except:
                                    # Sometimes there won't be any gates reporting, this is legacy behaviour and
                                    #  results in partial records
                                    gates_reporting = []

                                for gate in range(gate_range[0], gate_range[1], 1):
                                    epoch.append(epoch_here)
                                    date_time.append(date_time_here)
                                    slist.append(gate)
                                    beam.append(beam_here)
                                    tfreq.append(fitacf_data[record]['tfreq'])
                                    frang.append(fitacf_data[record]['frang'])
                                    rsep.append(fitacf_data[record]['rsep'])

                                    try:
                                        gate_index = gates_reporting.index(gate)
                                        if fitacf_data[record]['qflg'][gate_index] == 1 and \
                                                fitacf_data[record]['p_l'][gate_index] >= 3 and \
                                                fitacf_data[record]['gflg'][gate_index] == 0:
                                            good_echo.append(True)
                                        else:
                                            good_echo.append(False)
                                    except ValueError:
                                        good_echo.append(False)  # The gate is not in the list of reporting gates
                                    except BaseException as e:
                                        print(e)

    if len(epoch) == 0:
        # We have found no data, panic
        raise Exception("get_data_occ() found no data matching the provided criteria.")

    # Put the data into a dataframe
    print("     Building the data frame...")
    df = pd.DataFrame({'epoch': epoch, 'datetime': date_time,
                       'slist': slist, 'bmnum': beam,
                       'frang': frang, 'rsep': rsep, 'tfreq': tfreq,
                       'good_echo': good_echo
                       })

    return df


def build_datetime_epoch(year, month, day, hour, minute, second):
    """
    Build a datetime struct and compute epoch from raw date/time data.
    :param year: int: Y
    :param month: int: m
    :param day: d
    :param hour: H
    :param minute: M
    :param second: S
    :return: time.struct_time, int: The datetime and epoch
    """
    pattern = "%Y.%m.%d %H:%M:%S"  # This is the pattern we will use to convert time info to epoch
    datetime_here_str = str(year) + "." + str(month) + "." + str(day) + " " + \
                        str(hour) + ":" + str(minute) + ":" + str(second)

    date_time_struct = time.strptime(datetime_here_str, pattern)
    epoch = calendar.timegm(date_time_struct)
    date_time = datetime.datetime.fromtimestamp(time.mktime(date_time_struct))

    return date_time, epoch


if __name__ == '__main__':
    """ Testing """
    df = get_data_occ("sas", year_range=(2001, 2001), month_range=(1, 1), day_range=(1, 1),
                      gate_range=(0, 99), beam_range=(6, 7))

    print(df.head())