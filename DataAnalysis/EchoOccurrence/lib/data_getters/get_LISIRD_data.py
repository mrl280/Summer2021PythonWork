import pathlib
import datetime
import time
import calendar

import pandas as pd


def get_LISIRD_data(dataset_name):
    """
    Get LISIRD data

    Please see https://lasp.colorado.edu/lisird/ for a complete listing of lizard datasets
    When downloading LISIRD data, use '&format_time(yyyy-MM-dd'T'HH:mm:ss.SSS)', this is required to ensure a
    time (yyyy-MM-dd'T'HH:mm:ss.SSS) data field

    :precond: The dataset must exists as a csv file in the data folder
    :param dataset_name: The name of the LISIRD dataset, as a string (without the file extension)
    :return: pandas.Date_Frame: a dataframe with the requested data
    """

    # Read in the data from csv
    root = str((pathlib.Path().parent.absolute()))
    in_dir = root + "/data"
    in_file = in_dir + "/" + dataset_name + ".csv"
    try:
        df = pd.read_csv(in_file)
    except FileNotFoundError:
        # We might have imported from somewhere else, try again
        root = str((pathlib.Path().parent.absolute().parent.absolute()))
        in_dir = root + "/data"
        in_file = in_dir + "/" + dataset_name + ".csv"
        try:
            df = pd.read_csv(in_file)
        except FileNotFoundError:
            # We might have imported from somewhere else, try again
            root = str((pathlib.Path().parent.absolute().parent.absolute().parent.absolute()))
            in_dir = root + "/data"
            in_file = in_dir + "/" + dataset_name + ".csv"
            df = pd.read_csv(in_file)

    # Rename some the columns in common datasets
    df.rename(columns={"observed_flux (solar flux unit (SFU))": "observed_flux",
                       "adjusted_flux (solar flux unit (SFU))": "adjusted_flux",
                       "ssn (count)": "sunspot_number",
                       "time (yyyy-MM-dd'T'HH:mm:ss.SSS)": "time"},
              inplace=True)

    if 'time' not in df.columns:
        raise Exception("Error in get_LISIRD_data(), the expected time column 'time (yyyy-MM-dd'T'HH:mm:ss.SSS)' "
                        "was not found in " + str(dataset_name))

    # Add datetime and epoch to file
    date_time, epoch = [], []
    for i in range(len(df)):
        time_string = df['time'].iat[i]

        year = int(time_string[0:4])
        month = int(time_string[5:7])
        day = int(time_string[8:10])
        hour = int(time_string[11:13])

        date_time_here, epoch_here = build_datetime_epoch_local(year=year, month=month, day=day, hour=hour)

        date_time.append(date_time_here)
        epoch.append(epoch_here)

    df['datetime'] = date_time
    df['epoch'] = epoch

    df.drop(columns=["time"], inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df


def build_datetime_epoch_local(year, month, day, hour):
    """
    :param year: int: the year to consider
    :param month: int: the month to consider
    :param day: int: the day to consider
    :param hour: int: The hour to consider
    :return: time.struct_time, int: The datetime and epoch
    """
    pattern = '%Y.%m.%d %H:%M:%S'

    if month < 10:
        month = "0" + str(month)
    else:
        month = str(month)
    if day < 10:
        day = "0" + str(day)
    else:
        day = str(day)
    if hour < 10:
        hour = "0" + str(hour)
    else:
        hour = str(hour)

    if hour == "24":
        date_time_str = str(year) + "." + str(month) + "." + str(day) \
            + " " + "23:59:59"
    else:
        date_time_str = str(year) + "." + str(month) + "." + str(day) \
            + " " + str(hour) + ":00:00"
    date_time_struct = time.strptime(date_time_str, pattern)
    epoch = calendar.timegm(date_time_struct)

    date_time = datetime.datetime.fromtimestamp(time.mktime(date_time_struct))

    return date_time, epoch


if __name__ == "__main__":
    """ Testing """

    dataset_1 = "penticton_radio_flux_nearest_noon"
    dataset_2 = "international_sunspot_number"

    df = get_LISIRD_data(dataset_name=dataset_1)
    print(df.keys())
    print(df.head())

    df = get_LISIRD_data(dataset_name=dataset_2)
    print(df.keys())
    print(df.head())
