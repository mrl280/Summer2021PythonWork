import calendar
import glob
import math

import numpy as np
import time
import os
import pathlib
import statistics
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy import stats
from PyPDF2 import PdfFileMerger
from matplotlib.ticker import MultipleLocator

import pandas as pd

from DataAnalysis.DataReading.SD.basic_SD_df_filter import basic_SD_df_filter

if __name__ == '__main__':
    """
    Match up SuperDARN velocity points from multi-freq mode

    """

    year = "2016"  # yyyy
    month = "09"  # mm
    day = "25"  # dd
    gates = [40, 74]

    station = "rkn"
    start_date_time = year + "-" + month + "-" + day + " " + "00:00:00"
    end_date_time = year + "-" + month + "-" + day + " " + "04:00:00"

    # Compute start and end epoch
    pattern = '%Y-%m-%d %H:%M:%S'
    start_epoch = calendar.timegm(time.strptime(start_date_time, pattern))
    end_epoch = calendar.timegm(time.strptime(end_date_time, pattern))

    # Read in SuperDARN data
    loc_root = str(((pathlib.Path().parent.absolute()).parent.absolute()).parent.absolute())
    SD_in_dir = loc_root + "/DataReading/SD/data/" + station + "/" + station + year + month + day
    SD_in_file = SD_in_dir + "/" + station + year + month + day + ".pkl"
    df = pd.read_pickle(SD_in_file)

    # We are only interested in 15 km resolution data
    # This double filter should be redundant but better to be safe
    df = df.loc[(df['firstRang'] == 90) & (df['rangeSep'] == 15)]
    df = basic_SD_df_filter(df)  # This re-indexes

    # All data for this experiment is in one beam
    beam = df['bmnum'][0]

    # Filter for a gate range and time
    df = df.loc[(df['gate'] >= gates[0]) & (df['gate'] <= gates[1])]
    df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]
    df.reset_index(drop=True, inplace=True)

    # Recover decimal time
    print("Recovering decimal time...")
    decimal_time = []
    for t in range(len(df['epoch'])):
        decimal_time.append(df['hour'][t] + df['minute'][t] / 60.0 + df['second'][t] / 3600.0)
    df['decimalTime'] = decimal_time
    print("Number of points in the data frame: " + str(len(df['epoch'])))

    # Initialize arrays
    start_decimal_times = []
    vel10 = []
    vel12 = []
    vel13 = []
    vel14 = []

    # Look through the gates in sets of 5
    for start_gate in range(gates[0], gates[1], 5):
        end_gate = start_gate + 4
        print("Looking at gates " + str(start_gate) + " to " + str(end_gate))

        # Restrict data to with that gate
        restricted_df = df[(df['gate'] >= start_gate) & (df['gate'] <= end_gate)]
        print("Number of points here: " + str(len(restricted_df['epoch'])))

        # Loop through decimal time in jumps of 30 seconds
        time_interval_s = 120  # seconds
        for start_time in np.arange(0, 4, time_interval_s / 3600):
            # TODO: Match up points.  There are lots of points in each gate but not each time.
            #  Might have to remove .loc
            end_time = start_time + time_interval_s / 3600
            restricted_df = restricted_df.loc[(restricted_df['decimalTime'] >= start_time)
                                              & (restricted_df['decimalTime'] < end_time)]
            start_decimal_times.append(start_time)

            df_10 = restricted_df.loc[(restricted_df['transFreq'] >= (10 - 0.4) * 1000)
                                      & (restricted_df['transFreq'] <= (10 + 0.4) * 1000)]
            try:
                vel10.append(statistics.median(df_10['vel']))
            except:
                vel10.append(math.nan)

            df_12 = restricted_df.loc[(restricted_df['transFreq'] >= (12 - 0.4) * 1000)
                                      & (restricted_df['transFreq'] <= (12 + 0.4) * 1000)]
            try:
                vel12.append(statistics.median(df_12['vel']))
            except:
                vel12.append(math.nan)

            df_13 = restricted_df.loc[(restricted_df['transFreq'] >= (13 - 0.4) * 1000)
                                      & (restricted_df['transFreq'] <= (13 + 0.4) * 1000)]
            try:
                vel13.append(statistics.median(df_13['vel']))
            except:
                vel13.append(math.nan)

            df_14 = restricted_df.loc[(restricted_df['transFreq'] >= (14 - 0.4) * 1000)
                                      & (restricted_df['transFreq'] <= (14 + 0.4) * 1000)]
            try:
                vel14.append(statistics.median(df_14['vel']))
            except:
                vel14.append(math.nan)

    # Put the data into a dataframe
    matched_data = pd.DataFrame({'startTimes': start_decimal_times,
                                 'vel10': vel10,
                                 'vel12': vel12,
                                 'vel13': vel13,
                                 'vel14': vel14})

    print(matched_data[matched_data['vel12'].notnull()])

    # print(df['gate'].unique())
    #
    # print(df[['hour', 'minute', 'second', 'bmnum', 'gate', 'transFreq']])
    # print("Data start time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[0])))
    # print("Data end time: " + time.strftime(pattern, time.gmtime(df['epoch'].iloc[df.shape[0] - 1])))
