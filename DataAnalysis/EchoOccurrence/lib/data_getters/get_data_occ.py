import bz2
import glob
import os
import calendar
import pathlib
import time
import warnings
import pydarn

import datetime as datetime
import pandas as pd
import numpy as np
import _pickle as cPickle


def get_data_occ(station, year_range, month_range, day_range, gate_range, beam_range, freq_range,
                 fitACF_version=2.5, even_odd_days=None, print_timing_info=False):
    """
    Take a fitACF file and for each possible echo, record whether or not there was a good echo there
    Note: occurrence rate is not actually computed here, but all the data required to compute it is put into the df

    From Koustov and Syd's paper:
    "The echo occurrence rate was computed as a ratio of the number of registered echoes in selected beams and gates
    to the total number of possible echo detections in the same gates over the same period."

    :param station: str:
            The radar station to consider, as a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year_range: (int, int):
            Inclusive. The year range to consider.
    :param month_range: (int, int):
            Inclusive. The months of the year to consider.  If omitted (or None), then all days will be considered.
    :param day_range: (int, int):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param gate_range: (int, int):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (int, int):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param freq_range: (float, float):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.

    :param fitACF_version: float: (optional; default is 2.5):
            The fitACF version number.  At the time of writing the default is 2.5, but expected to move to 3.0 in the
            near future.  These are the only valid options.
    :param print_timing_info: bool (optional; default is False)
            Print out some time information for performance testing
    :param even_odd_days: (optional; default is None)
            'even': only even days are read in
            'odd': only odd days are read in
            None: all days are read in

    :return: pandas.DataFrame:
            A dataframe with select fitACF data.
    """

    if fitACF_version == 3.0:
        loc_root = "/data/fitacf_30"
    elif fitACF_version == 2.5:
        loc_root = "/data/fitacf_25"
    else:
        raise Exception("get_data_occ(): fitACF_version " + str(fitACF_version) +
                        " not recognized.  Right now only versions 2.5 and 3.0 are recognized.")

    if even_odd_days is not None:
        if even_odd_days == "even":
            warnings.warn("get_data_occ() is only considering even days of the month", category=Warning)
        elif even_odd_days == "odd":
            warnings.warn("get_data_occ() is only considering odd days of the month", category=Warning)
        else:
            raise Exception("get_data_occ(): even_odd_days: " + str(even_odd_days) +
                            " not recognized.  Options are 'even', 'odd', or None.")

    # Create empty arrays for the parameters we need
    epoch, date_time = [], []
    slist, beam = [], []
    frang, rsep, tfreq = [], [], []
    good_iono_echo, good_grndscat_echo = [], []

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
                            if even_odd_days is not None:
                                if even_odd_days == "even" and day_here % 2 == 1:
                                    # We only want even days, but this day is odd
                                    continue
                                elif even_odd_days == "odd" and day_here % 2 == 0:
                                    # We only want odd days, but this day is even
                                    continue

                            # We will read in the whole day, hour restrictions are left to the caller encase the care
                            #  to restrict based on something other than UT
                            print("    Reading: " + str(in_file))
                            try:
                                t0 = time.time()
                                with bz2.open(in_file) as fp:
                                    fitacf_stream = fp.read()
                                t1 = time.time()
                                sdarn_read = pydarn.SuperDARNRead(fitacf_stream, True)
                                t2 = time.time()
                                fitacf_data = sdarn_read.read_fitacf()  # this is
                                t3 = time.time()
                            except BaseException:
                                # Sometimes files are corrupted, or there is something wrong with them
                                pass

                            # Loop through every record in this file
                            for record in range(len(fitacf_data)):
                                beam_here = fitacf_data[record]['bmnum']
                                if beam_here < beam_range[0] or beam_here > beam_range[1]:
                                    continue  # We are not interested in records for this beam

                                date_time_here, epoch_here = build_datetime_epoch_local(
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
                                                fitacf_data[record]['p_l'][gate_index] >= 3:
                                            # We have a good echo
                                            if fitacf_data[record]['gflg'][gate_index] == 0:
                                                # We have a good ionospheric echo
                                                good_iono_echo.append(True)
                                                good_grndscat_echo.append(False)
                                            else:
                                                # We have a good ground scatter echo
                                                good_grndscat_echo.append(True)
                                                good_iono_echo.append(False)

                                        else:
                                            # Bad echo
                                            good_iono_echo.append(False)
                                            good_grndscat_echo.append(False)
                                    except ValueError:
                                        # The gate is not in the list of reporting gates
                                        good_iono_echo.append(False)
                                        good_grndscat_echo.append(False)
                                    except BaseException as e:
                                        print(e)

                            t4 = time.time()
                            if print_timing_info:
                                print("Time for unzip: " + str(t1 - t0))
                                print("Time for SuperDARN read: " + str(t2 - t1))
                                print("Time for fitacf conversion: " + str(t3 - t2))
                                print("Time for my stuff: " + str(t4 - t3))

    if len(epoch) == 0:
        # We have found no data, panic
        warnings.warn("get_data_occ() found no data matching the provided criteria.  "
                      "year_range: " + str(year_range) + ", month_range:" + str(month_range) +
                      ", day_range" + str(day_range), category=Warning)

    # Put the data into a dataframe
    print("     Building the data frame...")
    df = pd.DataFrame({'epoch': epoch, 'datetime': date_time,
                       'slist': slist, 'bmnum': beam,
                       'frang': frang, 'rsep': rsep, 'tfreq': tfreq,
                       'good_iono_echo': good_iono_echo, 'good_grndscat_echo': good_grndscat_echo
                       })

    # Ensure data is filtered for the the needed beam, gate, and freq ranges
    df = df.loc[(df['bmnum'] >= beam_range[0]) & (df['bmnum'] <= beam_range[1]) &  # Redundant - just to be safe
                (df['slist'] >= gate_range[0]) & (df['slist'] <= gate_range[1]) &  # Redundant - just to be safe

                # Note: freq_range is in MHz while data in 'tfreq' is in kHz
                (df['tfreq'] >= freq_range[0] * 1000) & (df['tfreq'] <= freq_range[1] * 1000) &

                # Only keep 45 km resolution data
                (df['frang'] == 180) & (df['rsep'] == 45)]

    # Remove everything we can
    df.drop(columns=['frang', 'rsep'], inplace=True)

    df.reset_index(drop=True, inplace=True)

    return df


def build_datetime_epoch_local(year, month, day, hour, minute, second):
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
    """ 
    Testing
    
    Also, to prebuild and pickle occurrence files for each reading
    When pre-building, it is best to use all of everything: all gates, beams, frequencies, months, and day
    """

    # testing = False
    # station = "dcn"
    # freq_range = (5, 25)
    # year_range = (2019, 2021)
    # beam_range = (6, 8)
    # gate_range = (10, 30)
    # fitACF_version = 2.5
    # even_odd_days = None
    # month_range = (1, 12)
    # day_range = (1, 31)

    # testing = False
    # station = "dce"
    # freq_range = (5, 25)
    # year_range = (2013, 2021)
    # beam_range = (6, 8)
    # gate_range = (10, 30)
    # fitACF_version = 2.5
    # even_odd_days = None
    # month_range = (1, 12)
    # day_range = (1, 31)

    testing = False
    station = "mcm"
    freq_range = (5, 25)
    year_range = (2013, 2021)
    beam_range = (6, 8)
    gate_range = (10, 30)
    fitACF_version = 2.5
    even_odd_days = None
    month_range = (1, 12)
    day_range = (1, 31)

    # testing = False
    # station = "dcn"
    # freq_range = (8, 10)
    # year_range = (2019, 2021)
    # beam_range = (6, 8)
    # gate_range = (10, 30)
    # fitACF_version = 2.5
    # even_odd_days = 'odd'
    # month_range = (1, 12)
    # day_range = (1, 31)

    # station = "dce"
    # freq_range = (5, 25)

    df = get_data_occ(station=station, year_range=year_range, month_range=month_range, day_range=day_range,
                      gate_range=gate_range, beam_range=beam_range, freq_range=freq_range,
                      fitACF_version=fitACF_version, even_odd_days=even_odd_days, print_timing_info=False)

    # df = get_data_occ("sas", year_range=(2001, 2001), month_range=(1, 1), day_range=(1, 1),
    #                   gate_range=(0, 99), beam_range=(6, 7), freq_range=(5, 25),
    #                   print_timing_info=True)

    print(df.head())

    if not testing:
        # Go ahead and pickle the dataframe for later
        if fitACF_version == 3.0:
            fitACT_string = "fitacf_30"
        elif fitACF_version == 2.5:
            fitACT_string = "fitacf_25"
        else:
            raise Exception("fitACF_version " + str(fitACF_version) + " not recognized.")

        loc_root = str((pathlib.Path().parent.absolute().parent.absolute().parent.absolute()))
        out_dir = loc_root + "/data/" + station
        out_file = out_dir + "/" + station + "_" + str(year_range[0]) + "_" + str(year_range[1]) + "_" + fitACT_string \
                   + "_occ-" + str(freq_range[0]) + "to" + str(freq_range[1]) + "MHz.pbz2"

        print("     Pickling as " + out_file + "...")
        with bz2.BZ2File(out_file, "w") as file:
            cPickle.dump(df, file)
