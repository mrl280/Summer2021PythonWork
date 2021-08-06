import calendar
import pydarn
import bz2
import time
import glob
import os

import pandas as pd
import numpy as np
import datetime as datetime
import dill as pickle


def PickleFITACF_occ(station, date, beam_range):
    """
    Take a fitACF file and for each possible echo, record whether or not there was a good echo there
    Note: occurrence rate is not actually computed here, but all the data required to compute it is put into the df

    From Koustov and Syd's paper
    "The echo occurrence rate was computed as a ratio of the number of registered echoes in selected beams and gates
    to the total number of possible echo detections in the same gates over the same period. 15-min averaging
    intervals were considered."

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param date: str:
            The date as a string of the form 'yyyymmdd'.  Single digit months and days need to be zero padded.
    :param beam_range: (int, int) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.

            In the past the following beams have been used:
            - INV 13-15     - DCE 10-12
            - RKN 1-3       - MCM 6-8
            - CLY 4-6       - SPS 3-5
    """

    gate_range = (0, 74)  # TODO: Update this to read in gate range from hardware?

    epoch, date_time = [], []
    slist, beam = [], []
    frang, rsep, tfreq = [], [], []
    good_iono_echo, good_grndscat_echo = [], []

    # Loop through all the files for this station/date
    in_dir = "data/" + station + "/" + station + date
    for in_file in glob.iglob(in_dir + "/*.fitacf.bz2"):
        # Unpack and open the file
        print("     Reading " + in_file + "...")
        try:
            with bz2.open(in_file) as fp:
                fitacf_stream = fp.read()
            sdarn_read = pydarn.SuperDARNRead(fitacf_stream, True)
            fitacf_data = sdarn_read.read_fitacf()
        except BaseException as e:
            # Sometimes files are corrupted, or there is something wrong with them
            print(e)
            pass
        # print("List of available parameters: " + str(fitacf_data[0].keys()))

        # Loop through every record in this file
        for record in range(len(fitacf_data)):
            beam_here = fitacf_data[record]['bmnum']
            if beam_here < beam_range[0] or beam_here > beam_range[1]:
                continue  # We are not interested in records for this beam

            date_time_here, epoch_here = build_datetime_epoch(year=fitacf_data[record]['time.yr'],
                                                              month=fitacf_data[record]['time.mo'],
                                                              day=fitacf_data[record]['time.dy'],
                                                              hour=fitacf_data[record]['time.hr'],
                                                              minute=fitacf_data[record]['time.mt'],
                                                              second=fitacf_data[record]['time.sc'])

            try:
                gates_reporting = np.ndarray.tolist(fitacf_data[record]['slist'])
            except:
                # Sometimes there won't be any gates reporting, this is legacy behaviour and results in partial records
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
                    if fitacf_data[record]['qflg'][gate_index] == 1 and fitacf_data[record]['p_l'][gate_index] >= 3:
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

    if len(epoch) == 0:
        # We have found no data, panic
        raise Exception("PickleFITACF_occ() found no data matching the provided criteria.")

    # Put the data into a dataframe
    print("     Building the data frame...")
    df = pd.DataFrame({'epoch': epoch,      'datetime': date_time,
                       'slist': slist,      'bmnum': beam,
                       'frang': frang,      'rsep': rsep,           'tfreq': tfreq,
                       'good_iono_echo': good_iono_echo,            'good_grndscat_echo': good_grndscat_echo
                       })

    # Save to file
    out_file = in_dir + "/" + station + date + "_occ.pbz2"
    print("     Pickling as " + out_file + "...")
    with bz2.BZ2File(out_file, "w") as file:
        pickle.dump(df, file)


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
    """
    Handler to call PickleFITACF on SuperDARN data files
    """
    PICKLE_ALL = False  # To prevent accidentally pickling all data

    if PICKLE_ALL:
        print("Occ Pickling all downloaded SuperDARN data...")
        for station in os.listdir("data/"):
            for in_dir in os.listdir("data/" + station):
                print("\nStarting " + in_dir)
                PickleFITACF_occ(station, in_dir[3:], (0, 15))
    else:
        station = "dce"
        date = "20190303"
        print("Occ Pickling " + station + date + "...")
        PickleFITACF_occ(station, date, (0, 15))
