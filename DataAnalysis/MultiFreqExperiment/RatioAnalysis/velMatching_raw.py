import calendar
import math
import os
import time
import pathlib
import bz2

import _pickle as cPickle
import numpy as np
import pandas as pd

from lib.basic_SD_df_filter import basic_SD_df_filter
from DataAnalysis.DataReading.SD.elevation_v2 import elevation_v2

if __name__ == '__main__':
    """
    Match up SuperDARN velocity points from multi-freq mode
    Use raw matching - it is thought that median filtering smooths out velocity differences
    """

    station = "rkn"
    year = "2016"  # yyyy
    month = "09"  # mm
    day = "26"  # dd
    area = "3c"  # options: None, 1, 2, 3, 4, 5, "3a", "3b", "3c"
    start_hour = 0  # Start and end times must be integer values for the loop
    end_hour = 4
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

    # Put frequencies in MHz and round, this makes them easier to compare
    df['transFreq_MHz'] = round(df['transFreq'] * 1e-3, 0)

    # print(df[['epoch', 'minute', 'second', 'gate', 'transFreq', 'vel', 'adjElv']].head())
    time_interval_s = 14  # Max allowed separation between points

    # Initialize arrays to hold matched data
    matched_times = []
    matched_gates = []
    vel10, elv10, height10 = [], [], []
    vel12, elv12, height12 = [], [], []
    vel13, elv13, height13 = [], [], []
    vel14, elv14, height14 = [], [], []

    gate_resolution = 1  # We need to maintain 15 km resolution

    # Loop through all the gates
    for start_gate in range(gates[0], gates[1] + 1, gate_resolution):
        print("Looking at gate " + str(start_gate))

        # Build a gate restricted data frame
        gg_df = df[(df['gate'] == start_gate)]
        gg_df.reset_index(drop=True, inplace=True)
        print("     Number of points at this gate: " + str(len(gg_df['epoch'])))

        range_here = 90 + 15 * start_gate + 15 / 2  # Slant range [km]
        used = [False] * len(gg_df['epoch'])  # Flag points that have already been matched up

        # Loop through all the rows at this gate
        for t in range(len(gg_df['epoch'])):

            # Find the first "unused" row
            if used[t]:
                continue

            used[t] = True  # We have found a point to use
            matched_times.append(gg_df['decimalTime'][t] + 4.5 / 3600)  # Best guess
            matched_gates.append(start_gate)
            found10, found12, found13, found14 = False, False, False, False

            # Find out what frequency our point is, and put it in it's place
            if gg_df['transFreq_MHz'][t] == 10:
                vel10.append(gg_df['vel'][t])
                elv10.append(gg_df['adjElv'][t])
                height10.append(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                        * np.sin(np.radians(np.asarray(gg_df['adjElv'][t])))) - Re)
                found10 = True
            elif gg_df['transFreq_MHz'][t] == 12:
                vel12.append(gg_df['vel'][t])
                elv12.append(gg_df['adjElv'][t])
                height12.append(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                        * np.sin(np.radians(np.asarray(gg_df['adjElv'][t])))) - Re)
                found12 = True
            elif gg_df['transFreq_MHz'][t] == 13:
                vel13.append(gg_df['vel'][t])
                elv13.append(gg_df['adjElv'][t])
                height13.append(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                        * np.sin(np.radians(np.asarray(gg_df['adjElv'][t])))) - Re)
                found13 = True
            elif gg_df['transFreq_MHz'][t] == 14:
                vel14.append(gg_df['vel'][t])
                elv14.append(gg_df['adjElv'][t])
                height14.append(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                        * np.sin(np.radians(np.asarray(gg_df['adjElv'][t])))) - Re)
                found14 = True
            else:
                raise Exception("Error: found a point with an unrecognized frequency: " + str(gg_df['transFreq'][t]))

            # Look at the next three points (dt=1, 2, 3), we might use them
            for dt in range(1, 4):
                if t + dt < len(gg_df['epoch']):  # So as to not run off the edge

                    if (abs(gg_df['epoch'][t + dt] - gg_df['epoch'][t]) < time_interval_s) \
                            and (gg_df['transFreq_MHz'][t + dt] != gg_df['transFreq_MHz'][t]):

                        used[t + dt] = True  # The point is close by and of a different frequency, it is a match, use it

                        # Find out what frequency it is, and put it in it's place
                        if gg_df['transFreq_MHz'][t + dt] == 10 and not found10:
                            vel10.append(gg_df['vel'][t + dt])
                            elv10.append(gg_df['adjElv'][t + dt])
                            height10.append(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                    * np.sin(np.radians(np.asarray(gg_df['adjElv'][t + dt])))) - Re)
                            found10 = True
                        elif gg_df['transFreq_MHz'][t + dt] == 12 and not found12:
                            vel12.append(gg_df['vel'][t + dt])
                            elv12.append(gg_df['adjElv'][t + dt])
                            height12.append(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                    * np.sin(np.radians(np.asarray(gg_df['adjElv'][t + dt])))) - Re)
                            found12 = True
                        elif gg_df['transFreq_MHz'][t + dt] == 13 and not found13:
                            vel13.append(gg_df['vel'][t + dt])
                            elv13.append(gg_df['adjElv'][t + dt])
                            height13.append(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                    * np.sin(np.radians(np.asarray(gg_df['adjElv'][t + dt])))) - Re)
                            found13 = True
                        elif gg_df['transFreq_MHz'][t + dt] == 14 and not found14:
                            vel14.append(gg_df['vel'][t + dt])
                            elv14.append(gg_df['adjElv'][t + dt])
                            height14.append(np.sqrt(Re * Re + range_here * range_here + 2 * Re * range_here
                                                    * np.sin(np.radians(np.asarray(gg_df['adjElv'][t + dt])))) - Re)
                            found14 = True
                        else:
                            # There were 2 measurements at the same frequency, ignore the second one
                            pass

            # Any frequencies that were not found, just fill in with nan
            if not found10:
                vel10.append(math.nan)
                elv10.append(math.nan)
                height10.append(math.nan)
            if not found12:
                vel12.append(math.nan)
                elv12.append(math.nan)
                height12.append(math.nan)
            if not found13:
                vel13.append(math.nan)
                elv13.append(math.nan)
                height13.append(math.nan)
            if not found14:
                vel14.append(math.nan)
                elv14.append(math.nan)
                height14.append(math.nan)

    # Put the data into a dataframe
    matched_data = pd.DataFrame({'decimalTime': matched_times,
                                 'gate': matched_gates,
                                 'vel10': vel10, 'elv10': elv10, 'height10': height10,
                                 'vel12': vel12, 'elv12': elv12, 'height12': height12,
                                 'vel13': vel13, 'elv13': elv13, 'height13': height13,
                                 'vel14': vel14, 'elv14': elv14, 'height14': height14,
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
        out_file = out_dir + "/" + station + year + month + day + ".RawMatchedData." + str(gate_resolution) + "gg" \
                   + str(time_interval_s) + "s.pbz2"
    else:
        out_file = out_dir + "/" + station + year + month + day + ".RawMatchedData." + str(gate_resolution) + "gg" \
                   + str(time_interval_s) + "s_area" + str(area) + ".pbz2"

    print("Pickling as " + out_file + "...")

    with bz2.BZ2File(out_file, "w") as file:
        cPickle.dump(matched_data, file)
