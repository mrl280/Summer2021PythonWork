import datetime
import time
import calendar


def build_date_epoch(year, month, day, hour):
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
    date_time, epoch = build_date_epoch(2018, 1, 1, 1)

    print(type(date_time))
