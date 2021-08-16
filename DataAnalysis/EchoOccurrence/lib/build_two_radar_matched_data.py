import numpy as np
import pandas as pd

try:
    # Program is being imported into a main level program
    from lib.compute_sd_radar_overlap import compute_df_radar_overlap
except:
    from DataAnalysis.EchoOccurrence.lib.compute_sd_radar_overlap import compute_df_radar_overlap


def build_two_radar_matched_data(station1, df1, gate_range, beam_range,
                                 station2, df2, other_gate_range, other_beam_range,
                                 time_interval_s):
    """

    Given two SuperDARN fit dataframes, build and return a dataframe holding matched data.

    Data is median matched - that is, we separate the data into temporal/spatial chunks, and take the median of all the
     data in a chunk
        - Temporal resolution - time_interval_s
        - Spatial resolution - kind of hard to say, but probably varies from 0 to 50 km?

    Note: This is just a straight chunk match - so we assume velocity (a line-of-sight vector measurement) has been
     adjusted so both radars are seeing velocities in the same direction

    # TODO: Sometimes a datapoint from one dataframe will be used multiple times because it will be the closest match
       to multiple different cells within the other dataframe.  Not sure if this is a valid approach, and should maybe
       be depicted on an overlap map.

    :param station1: str:
            The first radar station to consider, as 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param df1: pandas.DataFrame:
            The SuperDARN fit dataframe for station1.
    :param gate_range: (int, int):
            The range of allowable gates for station.  Gate matches from outside this range will not be considered.
            Often the near gates will not see F region echoes. So, you probably want to start
             around gate 20.
    :param beam_range: (int, int):
            The range of allowable beams for station.  Beam matches from outside this range will not be considered.

    :param station2: str:
            The second radar station to consider, again as a 3 character string.
    :param df2: pandas.DataFrame:
            The SuperDARN fit dataframe for station2.
            Needs to have temporal and spatial overlap with the first dataframe otherwise no matches will be found.
    :param other_gate_range: (int, int):
            The range of allowable gates for other_station.  Gate matches from outside this range will not be considered.
            Often the near gates will not see F region echoes. So, you probably want to start
             around gate 20.
    :param other_beam_range: (int, int):
            The range of allowable beams for other_station.  Beam matches from outside this range will not be considered.

    :param time_interval_s: int:
            The time resolution of the data.  Data is clumped into intervals of time_interval_s.

    :return matched_df: pandas.DataFrame:
            A dataframe with one row for each temporal/spatial chunk
            Depending on the event there might be lots of Nans
    """

    # Find the earliest and latest epoch values - these will be our starting and ending epochs
    df1_starting_epoch = np.min(df1['epoch'])
    df2_starting_epoch = np.min(df2['epoch'])

    df1_ending_epoch = np.max(df1['epoch'])
    df2_ending_epoch = np.max(df2['epoch'])

    if df1_starting_epoch < df2_starting_epoch:
        start_epoch = df1_starting_epoch
    else:
        start_epoch = df2_starting_epoch

    if df1_ending_epoch > df2_ending_epoch:
        end_epoch = df1_ending_epoch
    else:
        end_epoch = df2_ending_epoch

    epoch_edges = np.arange(start=start_epoch, stop=end_epoch, step=time_interval_s)
    # print("Epoch edges: " + str(epoch_edges))

    matched_epoch = []
    v1, v2 = [], []
    height1, height2 = [], []
    count1, count2 = [], []

    # Loop through every cell in the first stations dataframe - and match it up
    station1_gates = df1['slist'].unique()
    station1_beams = df1['bmnum'].unique()
    for station1_gate in station1_gates:
        for station1_beam in station1_beams:

            # Spatially restrict the data
            df1_ss = df1[(df1['bmnum'] == station1_beam) & (df1['slist'] == station1_gate)]

            station2_beam, station2_gate = compute_df_radar_overlap(station1=station1, station1_beam=station1_beam,
                                                                    station1_gate=station1_gate,
                                                                    gate_range1=gate_range,
                                                                    beam_range1=beam_range,
                                                                    station2=station2,
                                                                    gate_range2=other_gate_range,
                                                                    beam_range2=other_beam_range)

            df2_ss = df2[(df2['bmnum'] == station2_beam) & (df2['slist'] == station2_gate)]

            for starting_edge_epoch in epoch_edges:
                if starting_edge_epoch == epoch_edges[-1]:
                    continue  # The last edge is not a slice start

                ending_edge_epoch = starting_edge_epoch + time_interval_s

                df1_ss_tt = df1_ss[(df1_ss['epoch'] >= starting_edge_epoch) & (df1_ss['epoch'] <= ending_edge_epoch)]
                df2_ss_tt = df2_ss[(df2_ss['epoch'] >= starting_edge_epoch) & (df2_ss['epoch'] <= ending_edge_epoch)]

                matched_epoch.append(int(starting_edge_epoch + time_interval_s / 2))  # Middle of the time interval

                count1.append(df1_ss_tt.shape[0])
                count2.append(df2_ss_tt.shape[0])

                v1.append(np.median(df1_ss_tt['v']))
                v2.append(np.median(df2_ss_tt['v']))

                try:
                    height1.append(np.median(df1_ss_tt['height']))
                    height2.append(np.median(df2_ss_tt['height']))
                except KeyError:
                    # Guess we don't have this parameter
                    height1.append(np.nan)
                    height2.append(np.nan)
                except BaseException as e:
                    raise e

    matched_df = pd.DataFrame({'matched_epoch': matched_epoch,
                               'v1': v1, 'v2': v2,
                               'count1': count1, 'count2': count2,
                               'height1': height1, 'height2': height2,
                               })

    return matched_df
