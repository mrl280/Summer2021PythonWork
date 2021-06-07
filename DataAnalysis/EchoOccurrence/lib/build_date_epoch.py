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

    if hour == 24:
        datetime = str(year) + "." + str(month) + "." + str(day) \
            + " " + "23:59:59"
    else:
        datetime = str(year) + "." + str(month) + "." + str(day) \
            + " " + str(hour) + ":00:00"
    datetime = time.strptime(datetime, pattern)
    epoch = calendar.timegm(datetime)

    return datetime, epoch