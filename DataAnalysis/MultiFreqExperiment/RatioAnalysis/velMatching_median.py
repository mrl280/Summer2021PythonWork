import calendar
import math
import os
import time
import pathlib
import statistics
import bz2

import dill as pickle
import numpy as np
import pandas as pd

from lib.basic_SD_df_filter import basic_SD_df_filter
from lib.elevation_v2 import elevation_v2

if __name__ == '__main__':
    """
    Match up SuperDARN velocity points from multi-freq mode
    Use median filtering - it is thought that median filtering smooths out velocity differences
    """

    station = "rkn"
    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd
    area = "3c"  # options: None, 1, 2, 3, 4, 5, "3a", "3b", "3c"
    start_hour = 0  # Start and end times must be integer values for the loop
    end_hour = 4
    time_interval_s = 60  # seconds
    t_diff = 0.003  # Elevation angle correction in microseconds

    h_ratio_limits = [0.92, 1.08]  # Height ratio limits.  If a height ratio is outside of this range, it is flagged
    elv_ratio_limits = [0.88, 1.12]  # Elv ratio limits.  If an elev ratio is outside of this range, it is flagged

    start_date_time = year + "-" + month + "-" + day + " " + str(start_hour) + ":00:00"
    end_date_time = year + "-" + month + "-" + day + " " + str(end_hour) + ":00:00"
    gates = [10, 74]  # We will match up data over the whole gate range now and just restrict when plotting
    # Note: you can not use elevation angle data before gate 5ish because elevation is resolved in the 0-40 deg range
    Re = 6370  # Radius of the Earth, [km]

    # Compute start and end epoch
    pattern = '%Y-%m-%d %H:%M:%S'
    start_epoch = calendar.timegm(time.strptime(start_date_time, pattern))
    end_epoch = calendar.timegm(time.strptime(end_date_time, pattern))

    # Read in SuperDARN data
    if area is None:
        loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
        in_dir = loc_root + "/DataReading/SD/data/" + station + "/" + station + year + month + day
        in_file = in_dir + "/" + station + year + month + day + ".pbz2"
    else:
        # We are looking at area sectioned data for Sept 26, 2016
        loc_root = str((pathlib.Path().parent.absolute()).parent.absolute())
        in_dir = loc_root + "/Sept26EventDetailedInvestigation/data"
        in_file = in_dir + "/" + station + year + month + day + "_area" + str(area) + ".pbz2"

    print("Reading in file: " + in_file)
    data_stream = bz2.BZ2File(in_file, "rb")
    df = pd.read_pickle(data_stream)

    # We are only interested in 15 km resolution data (from the multi-freq analysis)
    # This double filter should be redundant but better to be safe
    df = df.loc[(df['firstRang'] == 90) & (df['rangeSep'] == 15)]
    df = basic_SD_df_filter(df)

    # Filter for a gate range and time
    df = df.loc[(df['gate'] >= gates[0]) & (df['gate'] <= gates[1])]
    df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]

    # There should not be any nan velocities, but just to be safe
    df = df.loc[df['vel'].notna()]

    # Drop points with |velocity| < 100 m/s
    df = df.loc[(df['vel'] > 100) | (df['vel'] < -100)]

    # Remove extreme values
    df = df.loc[(df['vel'] > -1000) & (df['vel'] < 1000)]

    # The different frequencies might be clumped together - so lets sort by epoch to make sure they are all mixed up
    df.sort_values(by='epoch', ascending=True, inplace=True)
    df.reset_index(drop=True, inplace=True)

    elevation_v2(df, t_diff)  # Compute adjusted elevation angle

    print("Number of points in the data frame: " + str(len(df['epoch'])))

    # Initialize arrays to hold matched data
    matched_times = []
    matched_gates = []
    vel10, elv10, count10, height10 = [], [], [], []
    vel12, elv12, count12, height12 = [], [], [], []
    vel13, elv13, count13, height13 = [], [], [], []
    vel14, elv14, count14, height14 = [], [], [], []

    time_step_h = time_interval_s / 3600
    gate_resolution = 1  # We need to maintain 15 km resolution

    # Look through the gates
    for start_gate in range(gates[0], gates[1] + 1, gate_resolution):
        end_gate = start_gate
        print("     Looking at gate " + str(start_gate))

        # Restrict data to within that gate
        gate_restricted_df = df[(df['gate'] >= start_gate) & (df['gate'] <= end_gate)]
        print("         Number of points at this gate: " + str(len(gate_restricted_df['epoch'])))

        range_here = 90 + 15 / 2 + 15 * start_gate  # Slant range [km]

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
                elv10.append(statistics.median(df_10['adjElv']))
                count10.append(df_10.shape[0])
                height10.append(statistics.median(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                          * np.sin(np.radians(np.asarray(df_10['adjElv']))))
                                                  - Re))
            except:
                vel10.append(math.nan)
                elv10.append(math.nan)
                count10.append(math.nan)
                height10.append(math.nan)

            df_12 = time_restricted_df[(time_restricted_df['transFreq'] >= (12 - 0.4) * 1000)
                                       & (time_restricted_df['transFreq'] <= (12 + 0.4) * 1000)]
            try:
                vel12.append(statistics.median(df_12['vel']))
                elv12.append(statistics.median(df_12['adjElv']))
                count12.append(df_12.shape[0])
                height12.append(statistics.median(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                          * np.sin(np.radians(np.asarray(df_12['adjElv']))))
                                                  - Re))
            except:
                vel12.append(math.nan)
                elv12.append(math.nan)
                count12.append(math.nan)
                height12.append(math.nan)

            df_13 = time_restricted_df[(time_restricted_df['transFreq'] >= (13 - 0.4) * 1000)
                                       & (time_restricted_df['transFreq'] <= (13 + 0.4) * 1000)]
            try:
                vel13.append(statistics.median(df_13['vel']))
                elv13.append(statistics.median(df_13['adjElv']))
                count13.append(df_13.shape[0])
                height13.append(statistics.median(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                          * np.sin(np.radians(np.asarray(df_13['adjElv']))))
                                                  - Re))
            except:
                vel13.append(math.nan)
                elv13.append(math.nan)
                count13.append(math.nan)
                height13.append(math.nan)

            df_14 = time_restricted_df[(time_restricted_df['transFreq'] >= (14 - 0.4) * 1000)
                                       & (time_restricted_df['transFreq'] <= (14 + 0.4) * 1000)]
            try:
                vel14.append(statistics.median(df_14['vel']))
                elv14.append(statistics.median(df_14['adjElv']))
                count14.append(df_14.shape[0])
                height14.append(statistics.median(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                          * np.sin(np.radians(np.asarray(df_14['adjElv']))))
                                                  - Re))
            except:
                vel14.append(math.nan)
                elv14.append(math.nan)
                count14.append(math.nan)
                height14.append(math.nan)

    # Put the data into a dataframe
    matched_data = pd.DataFrame({'decimalTime': matched_times,
                                 'gate': matched_gates,
                                 'vel10': vel10, 'elv10': elv10, 'count10': count10, 'height10': height10,
                                 'vel12': vel12, 'elv12': elv12, 'count12': count12, 'height12': height12,
                                 'vel13': vel13, 'elv13': elv13, 'count13': count13, 'height13': height13,
                                 'vel14': vel14, 'elv14': elv14, 'count14': count14, 'height14': height14
                                 })

    # Compute velocity ratios
    matched_data['10over12'] = matched_data['vel10'] / matched_data['vel12']
    matched_data['13over12'] = matched_data['vel13'] / matched_data['vel12']
    matched_data['14over12'] = matched_data['vel14'] / matched_data['vel12']

    matched_data['14over13'] = matched_data['vel14'] / matched_data['vel13']

    matched_data['13over10'] = matched_data['vel13'] / matched_data['vel10']
    matched_data['14over10'] = matched_data['vel14'] / matched_data['vel10']

    # Compute height ratios
    h_ratio_10over12 = matched_data['height10'] / matched_data['height12']
    h_ratio_13over12 = matched_data['height13'] / matched_data['height12']
    h_ratio_14over12 = matched_data['height14'] / matched_data['height12']

    h_ratio_14over13 = matched_data['height14'] / matched_data['height13']

    h_ratio_13over10 = matched_data['height13'] / matched_data['height10']
    h_ratio_14over10 = matched_data['height14'] / matched_data['height10']

    # Flag points where the y-axis frequency is larger
    matched_data['diffHeightFlag_10largerthan12'] = h_ratio_10over12 > h_ratio_limits[1]
    matched_data['diffHeightFlag_13largerthan12'] = h_ratio_13over12 > h_ratio_limits[1]
    matched_data['diffHeightFlag_14largerthan12'] = h_ratio_14over12 > h_ratio_limits[1]

    matched_data['diffHeightFlag_14largerthan13'] = h_ratio_14over13 > h_ratio_limits[1]

    matched_data['diffHeightFlag_13largerthan10'] = h_ratio_13over10 > h_ratio_limits[1]
    matched_data['diffHeightFlag_14largerthan10'] = h_ratio_14over10 > h_ratio_limits[1]

    # Flag points where the base (x-axis) frequency is larger
    matched_data['diffHeightFlag_10lessthan12'] = h_ratio_10over12 < h_ratio_limits[0]
    matched_data['diffHeightFlag_13lessthan12'] = h_ratio_13over12 < h_ratio_limits[0]
    matched_data['diffHeightFlag_14lessthan12'] = h_ratio_14over12 < h_ratio_limits[0]

    matched_data['diffHeightFlag_14lessthan13'] = h_ratio_14over13 < h_ratio_limits[0]

    matched_data['diffHeightFlag_13lessthan10'] = h_ratio_13over10 < h_ratio_limits[0]
    matched_data['diffHeightFlag_14lessthan10'] = h_ratio_14over10 < h_ratio_limits[0]

    # Flag points where both heights are about the same - this simplifies plotting
    matched_data['diffHeightFlag_10about12'] = np.logical_and(h_ratio_10over12 >= h_ratio_limits[0],
                                                              h_ratio_10over12 <= h_ratio_limits[1])
    matched_data['diffHeightFlag_13about12'] = np.logical_and(h_ratio_13over12 >= h_ratio_limits[0],
                                                              h_ratio_13over12 <= h_ratio_limits[1])
    matched_data['diffHeightFlag_14about12'] = np.logical_and(h_ratio_14over12 >= h_ratio_limits[0],
                                                              h_ratio_14over12 <= h_ratio_limits[1])

    matched_data['diffHeightFlag_14about13'] = np.logical_and(h_ratio_14over13 >= h_ratio_limits[0],
                                                              h_ratio_14over13 <= h_ratio_limits[1])

    matched_data['diffHeightFlag_13about10'] = np.logical_and(h_ratio_13over10 >= h_ratio_limits[0],
                                                              h_ratio_13over10 <= h_ratio_limits[1])
    matched_data['diffHeightFlag_14about10'] = np.logical_and(h_ratio_14over10 >= h_ratio_limits[0],
                                                              h_ratio_14over10 <= h_ratio_limits[1])

    # Compute elevation angle ratios
    elv_ratio_10over12 = matched_data['elv10'] / matched_data['elv12']
    elv_ratio_13over12 = matched_data['elv13'] / matched_data['elv12']
    elv_ratio_14over12 = matched_data['elv14'] / matched_data['elv12']

    elv_ratio_14over13 = matched_data['elv14'] / matched_data['elv13']

    elv_ratio_13over10 = matched_data['elv13'] / matched_data['elv10']
    elv_ratio_14over10 = matched_data['elv14'] / matched_data['elv10']

    # Flag points where the y-axis frequency is larger
    matched_data['diffElvFlag_10largerthan12'] = elv_ratio_10over12 > elv_ratio_limits[1]
    matched_data['diffElvFlag_13largerthan12'] = elv_ratio_13over12 > elv_ratio_limits[1]
    matched_data['diffElvFlag_14largerthan12'] = elv_ratio_14over12 > elv_ratio_limits[1]

    matched_data['diffElvFlag_14largerthan13'] = elv_ratio_14over13 > elv_ratio_limits[1]

    matched_data['diffElvFlag_13largerthan10'] = elv_ratio_13over10 > elv_ratio_limits[1]
    matched_data['diffElvFlag_14largerthan10'] = elv_ratio_14over10 > elv_ratio_limits[1]

    # Flag points where the base (x-axis) frequency is larger
    matched_data['diffElvFlag_10lessthan12'] = elv_ratio_10over12 < elv_ratio_limits[0]
    matched_data['diffElvFlag_13lessthan12'] = elv_ratio_13over12 < elv_ratio_limits[0]
    matched_data['diffElvFlag_14lessthan12'] = elv_ratio_14over12 < elv_ratio_limits[0]

    matched_data['diffElvFlag_14lessthan13'] = elv_ratio_14over13 < elv_ratio_limits[0]

    matched_data['diffElvFlag_13lessthan10'] = elv_ratio_13over10 < elv_ratio_limits[0]
    matched_data['diffElvFlag_14lessthan10'] = elv_ratio_14over10 < elv_ratio_limits[0]

    # Flag points where both heights are about the same - this simplifies plotting
    matched_data['diffElvFlag_10about12'] = np.logical_and(elv_ratio_10over12 >= elv_ratio_limits[0],
                                                           elv_ratio_10over12 <= elv_ratio_limits[1])
    matched_data['diffElvFlag_13about12'] = np.logical_and(elv_ratio_13over12 >= elv_ratio_limits[0],
                                                           elv_ratio_13over12 <= elv_ratio_limits[1])
    matched_data['diffElvFlag_14about12'] = np.logical_and(elv_ratio_14over12 >= elv_ratio_limits[0],
                                                           elv_ratio_14over12 <= elv_ratio_limits[1])

    matched_data['diffElvFlag_14about13'] = np.logical_and(elv_ratio_14over13 >= elv_ratio_limits[0],
                                                           elv_ratio_14over13 <= elv_ratio_limits[1])

    matched_data['diffElvFlag_13about10'] = np.logical_and(elv_ratio_13over10 >= elv_ratio_limits[0],
                                                           elv_ratio_13over10 <= elv_ratio_limits[1])
    matched_data['diffElvFlag_14about10'] = np.logical_and(elv_ratio_14over10 >= elv_ratio_limits[0],
                                                           elv_ratio_14over10 <= elv_ratio_limits[1])

    # Ensure out directory
    loc_root = str(pathlib.Path().parent.absolute())
    out_dir = loc_root + "/data/" + station + "/" + station + year + month + day
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Save the data for later
    if area is None:
        out_file = out_dir + "/" + station + year + month + day + ".MedianMatchedData." + str(gate_resolution) + "gg" \
                   + str(time_interval_s) + "s.pbz2"
    else:
        out_file = out_dir + "/" + station + year + month + day + ".MedianMatchedData." + str(gate_resolution) + "gg" \
                   + str(time_interval_s) + "s_area" + str(area) + ".pbz2"

    print("Pickling as " + out_file + "...")

    with bz2.BZ2File(out_file, "w") as file:
        pickle.dump(matched_data, file)
