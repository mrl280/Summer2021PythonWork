import calendar
import math
import os

import numpy as np
import time
import pathlib
import statistics

import pandas as pd

from DataAnalysis.DataReading.SD.basic_SD_df_filter import basic_SD_df_filter

if __name__ == '__main__':
    """
    Match up SuperDARN velocity points from multi-freq mode

    """
    # 2016 09 25 at RKN: gg [40, 74]
    # 4 pts, elv max 25 deg

    # 2019 09 26 at RKN: gg [10, 74]
    # 4 pts, elv max 27 deg

    # 2016 09 25 at CLY: gg [30, 74]
    # 4 pts, elv max 25 deg

    # 2016 09 26 at CLY: gg [10, 74]  # Not much here
    # 4 pts, elv max 25 deg

    # 2016 09 25 at INV: gg [10 74]
    # 4 pts, elv max 25

    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd

    station = "rkn"
    start_hour = 0  # Start and end times must be integer values for the loop
    end_hour = 4
    time_interval_s = 60  # seconds

    start_date_time = year + "-" + month + "-" + day + " " + str(start_hour) + ":00:00"
    end_date_time = year + "-" + month + "-" + day + " " + str(end_hour) + ":00:00"
    gates = [10, 74]  # We will match up data over the whole gate range now and just restrict when plotting
    Re = 6370  # Radius of the Earth, [km]

    # Compute start and end epoch
    pattern = '%Y-%m-%d %H:%M:%S'
    start_epoch = calendar.timegm(time.strptime(start_date_time, pattern))
    end_epoch = calendar.timegm(time.strptime(end_date_time, pattern))

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
    in_dir = loc_root + "/DataReading/SD/data/" + station + "/" + station + year + month + day
    in_file = in_dir + "/" + station + year + month + day + ".pkl"
    df = pd.read_pickle(in_file)

    # We are only interested in 15 km resolution data (from the multi-freq analysis)
    # This double filter should be redundant but better to be safe
    df = df.loc[(df['firstRang'] == 90) & (df['rangeSep'] == 15)]
    df = basic_SD_df_filter(df)

    # Filter for a gate range and time
    df = df.loc[(df['gate'] >= gates[0]) & (df['gate'] <= gates[1])]
    df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]

    # There should not be any nan velocities, but just to be safe
    df = df.loc[df['vel'].notna()]
    df.reset_index(drop=True, inplace=True)

    print("Number of points in the data frame: " + str(len(df['epoch'])))

    # Initialize arrays to hold matched data
    matched_times = []
    matched_gates = []
    vel10, count10, height10 = [], [], []
    vel12, count12, height12 = [], [], []
    vel13, count13, height13 = [], [], []
    vel14, count14, height14 = [], [], []

    time_step_h = time_interval_s / 3600
    # Look through the gates (we need to maintain 15 km resolution)
    for start_gate in range(gates[0], gates[1] + 1, 1):
        end_gate = start_gate
        print("     Looking at gate " + str(start_gate))

        # Restrict data to within that gate
        gate_restricted_df = df[(df['gate'] >= start_gate) & (df['gate'] <= end_gate)]
        print("         Number of points at this gate: " + str(len(gate_restricted_df['epoch'])))

        range_here = 90 + 15 * start_gate  # Slant range [km]

        # Loop through decimal time in jumps
        for start_time in np.arange(start_hour, end_hour, time_step_h):

            end_time = start_time + time_step_h
            time_restricted_df = gate_restricted_df[(gate_restricted_df['decimalTime'] >= start_time)
                                                    & (gate_restricted_df['decimalTime'] < end_time)]

            # Use the mid-time as the time marker
            matched_times.append(start_time + 0.5 * time_step_h)
            matched_gates.append(start_gate)

            df_10 = time_restricted_df[(time_restricted_df['transFreq'] >= (10 - 0.4) * 1000)
                                       & (time_restricted_df['transFreq'] <= (10 + 0.4) * 1000)]

            try:
                vel10.append(statistics.median(df_10['vel']))
                count10.append(df_10.shape[0])
                height10.append(statistics.median(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                          * np.sin(np.radians(np.asarray(df_10['elv']))))
                                                  - Re))
            except:
                vel10.append(math.nan)
                count10.append(math.nan)
                height10.append(math.nan)

            df_12 = time_restricted_df[(time_restricted_df['transFreq'] >= (12 - 0.4) * 1000)
                                       & (time_restricted_df['transFreq'] <= (12 + 0.4) * 1000)]
            try:
                vel12.append(statistics.median(df_12['vel']))
                count12.append(df_12.shape[0])
                height12.append(statistics.median(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                          * np.sin(np.radians(np.asarray(df_12['elv']))))
                                                  - Re))
            except:
                vel12.append(math.nan)
                count12.append(math.nan)
                height12.append(math.nan)

            df_13 = time_restricted_df[(time_restricted_df['transFreq'] >= (13 - 0.4) * 1000)
                                       & (time_restricted_df['transFreq'] <= (13 + 0.4) * 1000)]
            try:
                vel13.append(statistics.median(df_13['vel']))
                count13.append(df_13.shape[0])
                height13.append(statistics.median(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                          * np.sin(np.radians(np.asarray(df_13['elv']))))
                                                  - Re))
            except:
                vel13.append(math.nan)
                count13.append(math.nan)
                height13.append(math.nan)

            df_14 = time_restricted_df[(time_restricted_df['transFreq'] >= (14 - 0.4) * 1000)
                                       & (time_restricted_df['transFreq'] <= (14 + 0.4) * 1000)]
            try:
                vel14.append(statistics.median(df_14['vel']))
                count14.append(df_14.shape[0])
                height14.append(statistics.median(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                          * np.sin(np.radians(np.asarray(df_14['elv']))))
                                                  - Re))
            except:
                vel14.append(math.nan)
                count14.append(math.nan)
                height14.append(math.nan)

    # Put the data into a dataframe
    matched_data = pd.DataFrame({'decimalTime': matched_times,
                                 'gate': matched_gates,
                                 'vel10': vel10, 'count10': count10, 'height10': height10,
                                 'vel12': vel12, 'count12': count12, 'height12': height12,
                                 'vel13': vel13, 'count13': count13, 'height13': height13,
                                 'vel14': vel14, 'count14': count14, 'height14': height14
                                 })

    matched_data['10over12'] = matched_data['vel10'] / matched_data['vel12']
    matched_data['13over12'] = matched_data['vel13'] / matched_data['vel12']
    matched_data['14over12'] = matched_data['vel14'] / matched_data['vel12']

    matched_data['14over13'] = matched_data['vel14'] / matched_data['vel13']

    matched_data['13over10'] = matched_data['vel13'] / matched_data['vel10']
    matched_data['14over10'] = matched_data['vel14'] / matched_data['vel10']

    # Ensure out directory
    out_dir = loc_root + "/MultiFreqExperiment/VelocityAnalysis/data/" + station + "/" + station + year + month + day
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Save the data for later
    out_file = out_dir + "/" + station + year + month + day + ".MatchedData.1gg" + str(time_interval_s) + "s.pkl"
    print("Pickling as " + out_file + "...")
    matched_data.to_pickle(out_file)
