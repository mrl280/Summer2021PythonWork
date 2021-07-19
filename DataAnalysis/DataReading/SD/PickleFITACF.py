import math
import pydarn
import bz2
import glob
import os

import pandas as pd
import numpy as np
import _pickle as cPickle

from DataAnalysis.DataReading.SD.PickleFITACF_occ import build_datetime_epoch


def PickleFITACF(station, date):
    """
    Take a fitACF file, turn it into a reduced Pandas DataFrame, and pickle it for later
    More on fitACF here: https://radar-software-toolkit-rst.readthedocs.io/en/latest/references/general/fitacf/
    
    Takes a days worth of fitACF data and combine it all into one pickled data frame.
    
    Parameters of Interest:
        Scalar Parameters
            - Station id (three letter station identifier e.g. "rkn")
            - Epoch time (seconds since 1970.01.01)
                    This is a combination of year 'time.yr', month 'time.mo', day 'time.dy', 
                    hour 'time.hr', minute 'time.mt', second 'time.sc', and microsecond 'time.us'
            - Date time (human readable date time)
            - Year 'year', Month 'month', Day 'day'
            - Hour, Minute, Second (included microseconds)
            - Beam azimuth 'bmazm'
            - Integration Seconds (include microseconds) 'intt' 
            - Beam Number 'bmnum' (starting at Beam 0)
            - Transmitted Frequency 'tfreq' (in KHz)
            - Station Identifier 'stid'
            - Distance to the first range gate 'frang'
            - Range separation 'rsep' 
            - Control program name or command line, or comment buffer 'combf' 
            - Version number of the FITACF algorithm 'fitacf.revision.major'.'fitacf.revision.minor'
        
        Vector Parameters
            - Range Gate 'slist'
            - Ground Scatter Flag 'qflg'
            - Quality Flag 'gflg'
            - Velocity 'v' with Error 'v_e'
            - Power (actually SNR) from lambda fit 'p_l' with Error 'p_l_e'
            - Spectral width from lambda fit (w_l) 'w_l' with Error 'w_l_e'
            - Standard deviation of lambda fit 'sd_l' and phase fit 'sd_phi'
            - Elevation angle 'elv' with low 'elv_low' and high 'elv_high' estimates

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param date: str:
            The date as a string of the form 'yyyymmdd'.  Single digit months and days need to be zero padded.
    """

    # Create empty arrays for scalar parameters
    epoch, date_time = [], []
    year, month, day = [], [], []
    hour, minute, second = [], [], []
    date_time = []
    bmnum, bmazm = [], []
    intt = []
    tfreq = []
    frang, rsep = [], []
    combf, fitACF_rev = [], []

    # Create empty arrays for vector parameters
    slist = []
    qflg, gflg = [], []
    v, v_e = [], []
    p_l, p_l_e = [], []
    w_l, w_l_e = [], []
    sd_l, sd_phi = [], []
    phi0, phi0_e = [], []
    elv, elv_low, elv_high = [], [], []

    # Loop through all the files for this station/date
    in_dir = "data/" + station + "/" + station + date
    for in_file in glob.iglob(in_dir + "/*.fitacf.bz2"):

        # Unpack and open the file
        print("     Reading " + in_file + "...")
        with bz2.open(in_file) as fp:
            fitacf_stream = fp.read()
        sdarn_read = pydarn.SuperDARNRead(fitacf_stream, True)
        fitacf_data = sdarn_read.read_fitacf()
        # print("List of available parameters: " + str(fitacf_data[0].keys()))

        # Loop through every record in this file and add to the parallel arrays
        for record in range(len(fitacf_data)):

            date_time_here, epoch_here = build_datetime_epoch(year=fitacf_data[record]['time.yr'],
                                                              month=fitacf_data[record]['time.mo'],
                                                              day=fitacf_data[record]['time.dy'],
                                                              hour=fitacf_data[record]['time.hr'],
                                                              minute=fitacf_data[record]['time.mt'],
                                                              second=fitacf_data[record]['time.sc'])
            try:
                num_gates_reporting = len(fitacf_data[record]['slist'])
            except:
                # Sometimes there won't be any gates reporting, this is legacy behaviour and results in partial records
                num_gates_reporting = 0

            # Loop through all the reporting gates and build up the data
            for gate_idx in range(num_gates_reporting):
                # Build up vectored data
                slist.append(fitacf_data[record]['slist'][gate_idx])
                qflg.append(fitacf_data[record]['qflg'][gate_idx])
                gflg.append(fitacf_data[record]['gflg'][gate_idx])
                v.append(fitacf_data[record]['v'][gate_idx])
                v_e.append(fitacf_data[record]['v_e'][gate_idx])
                p_l.append(fitacf_data[record]['p_l'][gate_idx])
                p_l_e.append(fitacf_data[record]['p_l_e'][gate_idx])
                w_l.append(fitacf_data[record]['w_l'][gate_idx])
                w_l_e.append(fitacf_data[record]['w_l_e'][gate_idx])
                sd_l.append(fitacf_data[record]['sd_l'][gate_idx])
                sd_phi.append(fitacf_data[record]['sd_phi'][gate_idx])

                try:
                    phi0.append(fitacf_data[record]['phi0'][gate_idx])
                    phi0_e.append(fitacf_data[record]['phi0_e'][gate_idx])
                except KeyError:
                    # Some files might not have this parameter
                    phi0.append(math.nan)
                    phi0_e.append(math.nan)
                except BaseException:
                    raise

                try:
                    elv.append(fitacf_data[record]['elv'][gate_idx])
                    elv_low.append(fitacf_data[record]['elv_low'][gate_idx])
                    elv_high.append(fitacf_data[record]['elv_high'][gate_idx])
                except KeyError:
                    # Some files might not have this parameter
                    elv.append(math.nan)
                    elv_low.append(math.nan)
                    elv_high.append(math.nan)
                except BaseException:
                    raise

                # Build up scalar data, it is faster to have this in this inner loop
                epoch.append(epoch_here)
                date_time.append(date_time_here)
                year.append(fitacf_data[record]['time.yr'])
                month.append(fitacf_data[record]['time.mo'])
                day.append(fitacf_data[record]['time.dy'])
                hour.append(fitacf_data[record]['time.hr'])
                minute.append(fitacf_data[record]['time.mt'])
                second.append(fitacf_data[record]['time.sc'] + fitacf_data[record]['time.us'] / 1e6)
                bmnum.append(fitacf_data[record]['bmnum'])
                bmazm.append(fitacf_data[record]['bmazm'])
                intt.append(fitacf_data[record]['intt.sc'] + fitacf_data[record]['intt.us'] / 1e6)
                tfreq.append(fitacf_data[record]['tfreq'])
                frang.append(fitacf_data[record]['frang'])
                rsep.append(fitacf_data[record]['rsep'])
                fitACF_rev.append(str(fitacf_data[record]['fitacf.revision.major']) + "." + str(fitacf_data[record]['fitacf.revision.minor']))
                combf.append(fitacf_data[record]['combf'])

    # Put the data into a dataframe
    print("     Building the data frame...")
    df = pd.DataFrame({'station': [station] * len(epoch),
                       'datetime': date_time,       'epoch': epoch,
                       'decimalTime': np.asarray(hour) + np.asarray(minute) / 60.0 + np.asarray(second) / 3600.0,
                       'year': year,            'month': month,         'day': day,
                       'hour': hour,            'minute': minute,       'second': second,
                       'bmnum': bmnum,          'bmazm': bmazm,
                       'intt': intt,
                       'gate': slist,
                       'transFreq': tfreq,
                        'firstRang': frang,     'rangeSep': rsep,
                       'fitACFRev': fitACF_rev,
                       'CtrlPrgrm': combf,
                       'qflg': qflg,
                       'gflg': gflg,
                       'vel': v,            'velErr': v_e,
                       'pwr': p_l,          'pwrErr': p_l_e,
                       'wdt': w_l,          'wdtErr': w_l_e,
                       'phase': phi0,       'phaseErr': phi0_e,
                       'stdDevLambda': sd_l,
                       'stdDevPhase': sd_phi,
                       'elv': elv,          'elvLow': elv_low,          'elvHigh': elv_high
                       })

    # Compress and save
    out_file = in_dir + "/" + station + date + ".pbz2"
    print("     Pickling as " + out_file + "...")
    with bz2.BZ2File(out_file, "w") as file:
        cPickle.dump(df, file)


if __name__ == '__main__':
    """
    Handler to call PickleFITACF on SuperDARN data files
    """
    PICKLE_ALL = False  # To prevent accidentally pickling all data

    if PICKLE_ALL:
        print("Pickling all downloaded SuperDARN data...")
        for station in os.listdir("data/"):
            for in_dir in os.listdir("data/" + station):
                print("\nStarting " + in_dir)
                PickleFITACF(station, in_dir[3:])
    else:
        station = "rkn"
        date = "20111112"
        print("Pickling " + station + date + "...")
        PickleFITACF(station, date)
