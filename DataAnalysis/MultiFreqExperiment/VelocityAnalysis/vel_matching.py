import calendar
import math
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
    # 2016 09 25 at RKN: gg [40, 75]
    # 4 pts, elv max 25 deg

    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd
    gates = [30, 74]
    elv_max = 25  # deg

    station = "rkn"
    start_date_time = year + "-" + month + "-" + day + " " + "00:00:00"
    end_date_time = year + "-" + month + "-" + day + " " + "08:00:00"

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

    # There should not be any nan velocities, but drop to be safe
    df = df.loc[df['vel'].notna()]

    # Lets make sure we are looking at E region echoes by filtering out high elevation angle data
    df = df.loc[df['elv'] <= elv_max]
    df.reset_index(drop=True, inplace=True)

    # All data for this experiment is in one beam
    beam = df['bmnum'][0]

    # Recover decimal time
    print("Recovering decimal time...")
    df['decimalTime'] = df['hour'] + df['minute'] / 60.0 + df['second'] / 3600.0
    print("Number of points in the data frame: " + str(len(df['epoch'])))

    # Initialize arrays to hold matched data
    matched_times = []
    starting_gates = []
    vel10, count10 = [], []
    vel12, count12 = [], []
    vel13, count13 = [], []
    vel14, count14 = [], []

    # Look through the gates in sets of 5
    for start_gate in range(gates[0], gates[1], 2):
        end_gate = start_gate + 1
        print("Looking at gates " + str(start_gate) + " to " + str(end_gate))

        # Restrict data to with that gate
        gate_restricted_df = df[(df['gate'] >= start_gate) & (df['gate'] <= end_gate)]
        print("Number of points here: " + str(len(gate_restricted_df['epoch'])))

        # Loop through decimal time in jumps of 30 seconds
        time_interval_s = 30  # seconds
        for start_time in np.arange(0, 8, time_interval_s / 3600):

            end_time = start_time + time_interval_s / 3600
            time_restricted_df = gate_restricted_df[(gate_restricted_df['decimalTime'] >= start_time)
                                                    & (gate_restricted_df['decimalTime'] < end_time)]
            # print("     Number of points here: " + str(len(time_restricted_df['epoch'])))

            # Use the mid-time as the time marker
            matched_times.append(start_time + 0.5 * time_interval_s / 3600)
            starting_gates.append(start_gate)

            df_10 = time_restricted_df[(time_restricted_df['transFreq'] >= (10 - 0.4) * 1000)
                                       & (time_restricted_df['transFreq'] <= (10 + 0.4) * 1000)]
            try:
                vel10.append(statistics.median(df_10['vel']))
                count10.append(df_10.shape[0])
            except:
                vel10.append(math.nan)
                count10.append(math.nan)

            df_12 = time_restricted_df[(time_restricted_df['transFreq'] >= (12 - 0.4) * 1000)
                                       & (time_restricted_df['transFreq'] <= (12 + 0.4) * 1000)]
            try:
                vel12.append(statistics.median(df_12['vel']))
                count12.append(df_12.shape[0])
            except:
                vel12.append(math.nan)
                count12.append(math.nan)

            df_13 = time_restricted_df[(time_restricted_df['transFreq'] >= (13 - 0.4) * 1000)
                                       & (time_restricted_df['transFreq'] <= (13 + 0.4) * 1000)]
            try:
                vel13.append(statistics.median(df_13['vel']))
                count13.append(df_13.shape[0])
            except:
                vel13.append(math.nan)
                count13.append(math.nan)

            df_14 = time_restricted_df[(time_restricted_df['transFreq'] >= (14 - 0.4) * 1000)
                                       & (time_restricted_df['transFreq'] <= (14 + 0.4) * 1000)]
            try:
                vel14.append(statistics.median(df_14['vel']))
                count14.append(df_14.shape[0])
            except:
                vel14.append(math.nan)
                count14.append(math.nan)

    # Put the data into a dataframe
    matched_data = pd.DataFrame({'decimalTimes': matched_times,
                                 'startingGates': starting_gates,
                                 'vel10': vel10,
                                 'count10': count10,
                                 'vel12': vel12,
                                 'count12': count12,
                                 'vel13': vel13,
                                 'count13': count13,
                                 'vel14': vel14,
                                 'count14': count14
                                 })

    matched_data['10over12'] = matched_data['vel10'] / matched_data['vel12']
    matched_data['13over12'] = matched_data['vel13'] / matched_data['vel12']
    matched_data['14over12'] = matched_data['vel14'] / matched_data['vel12']
    # print(matched_data[matched_data['10over12'].notna()][['vel10', 'count10', 'vel12', 'count12', '10over12']])

    # Save the data for later
    out_dir = loc_root + "/MultiFreqExperiment/VelocityAnalysis/data/" + station
    out_file = out_dir + "/" + station + year + month + day + ".MatchedData.pkl"
    print("Pickling as " + out_file + "...")
    matched_data.to_pickle(out_file)


    # matched_data = matched_data.loc[(matched_data['count10'] >= 3) & (matched_data['count12'] >= 3)
    #                                 & matched_data['10over12'].notna()]
    # print(matched_data[['vel10', 'count10', 'vel12', 'count12', '10over12']])


    # df.loc[df['vel'].notna()]

    # print(df['gate'].unique())
    #
    # print(df[['hour', 'minute', 'second', 'bmnum', 'gate', 'transFreq']])
    # print("Data start time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[0])))
    # print("Data end time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[df.shape[0] - 1])))