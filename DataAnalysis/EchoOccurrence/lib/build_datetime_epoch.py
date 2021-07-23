import datetime
import time
import calendar


def build_datetime_epoch(year, month, day, hour, minute=None, second=None):
    """
    Build a datetime struct and compute epoch from raw date/time data.

    # TODO: Allow vector inputs

    :param year: int: Y
    :param month: int: m
    :param day: d
    :param hour: H
    :param minute: M (optional)
    :param second: S (optional)

    :return: time.struct_time, int: The datetime and epoch
    """

    pattern = "%Y.%m.%d %H:%M:%S"  # This is the pattern we will use to convert time info to epoch

    if minute is None:
        minute = "00"
    if second is None:
        second = "00"

    datetime_here_str = str(year) + "." + str(month) + "." + str(day) + " " + \
                        str(hour) + ":" + str(minute) + ":" + str(second)

    date_time_struct = time.strptime(datetime_here_str, pattern)
    epoch = calendar.timegm(date_time_struct)
    date_time = datetime.datetime.fromtimestamp(time.mktime(date_time_struct))

    return date_time, epoch


if __name__ == "__main__":
    """ Testing """
    date_time, epoch = build_datetime_epoch(year=2018, month=1, day=1, hour=1)
    print("Datetime: " + str(date_time))
    print("epoch: " + str(epoch))

    date_time, epoch = build_datetime_epoch(year=2018, month=1, day=1, hour=1, minute=6, second=3)
    print("Datetime: " + str(date_time))
    print("epoch: " + str(epoch))

    print(type(date_time))
