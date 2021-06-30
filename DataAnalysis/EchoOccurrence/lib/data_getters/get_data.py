import bz2
import glob
import math
import os
import pydarn
import time
import calendar

import pandas as pd
import datetime as datetime


def get_data(station, year_range, month_range, day_range, gate_range, beam_range, freq_range):
    """
    Get all of the SuperDARN data within the desired range/time, put it into a dataframe,
    and then return the dataframe for plotting/analysis

    Dataframe keys follow the same convention as fitACF:
        https://radar-software-toolkit-rst.readthedocs.io/en/latest/references/general/fitacf/

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
    :param freq_range: (<float>, <float>) (optional):
            Inclusive.  The frequency range to consider in MHz.
            If omitted (or None), then all frequencies are considered.
    :return: pandas.DataFrame: A dataframe with select fitACF parameters.
    """

    loc_root = "/data/fitacf_30"    # fitACF 3.0
    # loc_root = "/data/fitacf_25"  # fitACF 2.5

    # Create empty arrays for scalar parameters
    epoch, date_time = [], []
    bmnum = []
    tfreq = []
    frang, rsep = [], []

    # Create empty arrays for vector parameters
    slist = []
    qflg, gflg = [], []
    v, p_l, w_l = [], [], []
    phi0, elv = [], []

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

                            # Loop through all the scans and build up the arrays
                            # As far as I know, list extensions are the fastest way to do this
                            for scan in range(len(fitacf_data)):

                                date_time_here, epoch_here = build_datetime_epoch_local(
                                    year=fitacf_data[scan]['time.yr'],
                                    month=fitacf_data[scan]['time.mo'],
                                    day=fitacf_data[scan]['time.dy'],
                                    hour=fitacf_data[scan]['time.hr'],
                                    minute=fitacf_data[scan]['time.mt'],
                                    second=fitacf_data[scan]['time.sc'])

                                try:
                                    num_gates_reporting = len(fitacf_data[scan]['slist'])
                                except:
                                    # Sometimes there won't be any gates reporting, this is legacy behaviour and
                                    #  results in partial records
                                    num_gates_reporting = 0

                                if num_gates_reporting > 0:
                                    # Build up scalar parameters
                                    epoch.extend([epoch_here] * num_gates_reporting)
                                    date_time.extend([date_time_here] * num_gates_reporting)
                                    bmnum.extend([fitacf_data[scan]['bmnum']] * num_gates_reporting)
                                    tfreq.extend([fitacf_data[scan]['tfreq']] * num_gates_reporting)
                                    frang.extend([fitacf_data[scan]['frang']] * num_gates_reporting)
                                    rsep.extend([fitacf_data[scan]['rsep']] * num_gates_reporting)

                                    # Build up vector parameters
                                    slist.extend(fitacf_data[scan]['slist'])
                                    qflg.extend(fitacf_data[scan]['qflg'])
                                    gflg.extend(fitacf_data[scan]['gflg'])
                                    v.extend(fitacf_data[scan]['v'])
                                    p_l.extend(fitacf_data[scan]['p_l'])
                                    w_l.extend(fitacf_data[scan]['w_l'])
                                    try:
                                        phi0.extend(fitacf_data[scan]['phi0'])
                                    except KeyError:
                                        # Some files might not have this parameter
                                        phi0.extend([math.nan] * num_gates_reporting)
                                    except BaseException:
                                        raise
                                    try:
                                        elv.extend(fitacf_data[scan]['elv'])
                                    except KeyError:
                                        # Some files might not have this parameter
                                        elv.extend([math.nan] * num_gates_reporting)
                                    except BaseException:
                                        raise

    if len(epoch) == 0:
        # We have found no data, panic
        raise Exception("get_data() found no data matching the provided criteria.")

    df = pd.DataFrame({'epoch': epoch,      'datetime': date_time,
                       'bmnum': bmnum,
                       'tfreq': tfreq,
                       'frang': frang,      'rsep': rsep,

                       'slist': slist,
                       'qflg': qflg,        'gflg': gflg,
                       'v': v,              'p_l': p_l,                 'w_l': w_l,
                       'phi0': phi0,        'elv': elv
                       })

    # Filter the data for the needed time, beam, gate, and freq ranges
    df = df.loc[(df['bmnum'] >= beam_range[0]) & (df['bmnum'] <= beam_range[1]) &
                (df['slist'] >= gate_range[0]) & (df['slist'] <= gate_range[1]) &

                # Note: freq_range is in MHz while data in 'tfreq' is in kHz
                (df['tfreq'] >= freq_range[0] * 1000) & (df['tfreq'] <= freq_range[1] * 1000)]

    # Until I have an application that requires bad quality points, I will assume they always need to be filtered out
    df = df.loc[(df['qflg'] == 1) & (df['p_l'] >= 3)]

    df.drop(columns=['qflg'], inplace=True)
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
    """ Testing """
    df = get_data("rkn", year_range=(2001, 2001), month_range=(1, 1), day_range=(1, 1),
                  gate_range=(0, 99), beam_range=(0, 15), freq_range=(5, 25))

    print(df.head())
