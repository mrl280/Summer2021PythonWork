import pandas as pd
import pathlib
import calendar
import time

def get_df_multi_event(file_name=None, flag=None):
    """
    Get a dataframe with data for multiple events
    :param file_name: str: The name of the pickled data frame you would like to use.  e.g. "event_summary"
    :param flag: str: optional: only events that have this flag will be considered
    :return: a dataframe with data for multiple events
    """

    # Check the file_name
    if file_name is None or not isinstance(file_name, str):
        print("No file name given, so we are get_df_multi_event() is using event_summary.pkl")
        file_name = "event_summary.pkl"

    # Read in the data
    try:
        # This works if it is called from the same level
        loc_root = pathlib.Path().absolute()
        event_summary_dir = str(loc_root) + "/data"
        event_summary_path = event_summary_dir + "/" + file_name + ".pkl"
        event_summary = pd.read_pickle(event_summary_path)
    except:
        # This works if we are calling from a sub-folder
        loc_root = pathlib.Path().absolute().parent
        event_summary_dir = str(loc_root) + "/data"
        event_summary_path = event_summary_dir + "/" + file_name + ".pkl"
        event_summary = pd.read_pickle(event_summary_path)

    if flag is not None:
        # Only keep the rows that are flagged
        try:
            event_summary = event_summary.loc[event_summary[flag] == 1]
            event_summary.reset_index(drop=True, inplace=True)
        except:
            print("Error: " + str(flag) + " is not a recognized key.  "
                  "Your options are: " + str([i for i in event_summary.keys()]))

    pattern = '%Y-%m-%d %H:%M:%S'

    number_of_events = len(event_summary)
    if number_of_events <= 0:
        print("There are no events to consider.")
        return pd.DataFrame()
    else:
        # There is at least one event, use the first event to initialize the dataframe
        station = event_summary['SD_station'][0]
        year = str(event_summary['year'][0])

        if event_summary['month'][0] >= 10:
            month = str(event_summary['month'][0])
        else:
            month = "0" + str(event_summary['month'][0])

        if event_summary['day'][0] >= 10:
            day = str(event_summary['day'][0])
        else:
            day = "0" + str(event_summary['day'][0])

        start_hour_UT = event_summary['start_hour_UT'][0]
        end_hour_UT = event_summary['end_hour_UT'][0]

        # Compute start and end epoch, they are needed to filter the df
        start_date_time = year + "-" + month + "-" + day + " " + str(start_hour_UT) + ":00:00"
        end_date_time = year + "-" + month + "-" + day + " " + str(end_hour_UT) + ":00:00"
        start_epoch = calendar.timegm(time.strptime(start_date_time, pattern))
        end_epoch = calendar.timegm(time.strptime(end_date_time, pattern))

        in_dir = str(loc_root.parent) + "/DataReading/SD/data/" + station + "/" + station + year + month + day
        in_file = in_dir + "/" + station + year + month + day + ".pkl"
        df = pd.read_pickle(in_file)
        df = df.loc[(df['epoch'] >= start_epoch) & (df['epoch'] <= end_epoch)]
        df.reset_index(drop=True, inplace=True)

        # Loop through every other event and add the data to the dataframe
        for event in range(number_of_events):

            station = event_summary['SD_station'][event]
            year = str(event_summary['year'][event])

            if event_summary['month'][0] >= 10:
                month = str(event_summary['month'][0])
            else:
                month = "0" + str(event_summary['month'][0])

            if event_summary['day'][0] >= 10:
                day = str(event_summary['day'][0])
            else:
                day = "0" + str(event_summary['day'][0])

            start_hour_UT = event_summary['start_hour_UT'][event]
            end_hour_UT = event_summary['end_hour_UT'][event]

            # Compute start and end epoch, they are needed to filter the df
            start_date_time = year + "-" + month + "-" + day + " " + str(start_hour_UT) + ":00:00"
            end_date_time = year + "-" + month + "-" + day + " " + str(end_hour_UT) + ":00:00"
            start_epoch = calendar.timegm(time.strptime(start_date_time, pattern))
            end_epoch = calendar.timegm(time.strptime(end_date_time, pattern))

            in_dir = str(loc_root.parent) + "/DataReading/SD/data/" + station + "/" + station + year + month + day
            in_file = in_dir + "/" + station + year + month + day + ".pkl"
            next_df = pd.read_pickle(in_file)
            next_df = next_df.loc[(next_df['epoch'] >= start_epoch) & (next_df['epoch'] <= end_epoch)]
            next_df.reset_index(drop=True, inplace=True)

            df = pd.concat([df, next_df])  # Warning: this is costly if there are a large number of events

        return df


if __name__ == '__main__':
    """
    """

    df = get_df_multi_event(file_name="event_summary", flag='for_vel_hist')
    print(df)

