import datetime
import warnings


def check_year_range(year_range):
    """
    :param year_range: (<int>, <int>) or None: The year range to check
    :return: (<int>, <int>): A suitable year range
    """
    now = datetime.datetime.now()
    if year_range is None:
        year_range = (1993, now.year)  # Assume we all data
        warnings.warn("No year range was provided, it has defaulted to (1993, " + str(now.year) + ").",
                      category=Warning)
    if year_range[0] < 1993:
        year_range = (1993, year_range[1])
        warnings.warn("There is no SuperDARN data before 1993, year_range[0] has defaulted to 1993.",
                      category=Warning)
    if year_range[1] > now.year:
        year_range = (year_range[0], now.year)
        warnings.warn("It is only " + str(now.year) + ", year_range[1] has defaulted to " + str(str(now.year)) + ".",
                      category=Warning)
    return year_range


def check_year(year):
    """
    :param year: int or None: The year to check
    :return: int: A suitable year
    """
    now = datetime.datetime.now()
    if year is None:
        year = now.year
        warnings.warn("No year was provided, it has defaulted to " + str(now.year) + ".",
                      category=Warning)
    if year < 1993:
        year = 1993
        warnings.warn("There is no SuperDARN data before 1993, year has defaulted to 1993.",
                      category=Warning)
    if year > now.year:
        year = now.year
        warnings.warn("It is only " + str(now.year) + ", year has defaulted to " + str(year) + ".",
                      category=Warning)
    return year


def check_month_range(month_range):
    """
    :param month_range: (<int>, <int>) or None: The month range to check
    :return: (<int>, <int>): A suitable month range
    """
    if month_range is None:
        month_range = (1, 12)  # Assume we want all months
        warnings.warn("No month range was provided, it has defaulted to (1, 12).",
                      category=Warning)
    if month_range[0] < 1:
        month_range = (1, month_range[1])
        warnings.warn("month_range[0] has defaulted to 1 (Jan).",
                      category=Warning)
    if month_range[1] > 12:
        month_range = (month_range[0], 12)
        warnings.warn("month_range[1] has defaulted to 12 (Dec).",
                      category=Warning)
    return month_range


def check_day_range(day_range):
    """
    :param day_range: (<int>, <int>) or None: The day range to check
    :return: (<int>, <int>): A suitable day range
    """
    if day_range is None:
        day_range = (1, 31)  # Assume we want all days
        warnings.warn("No day range was provided, it has defaulted to (1, 31).",
                      category=Warning)
    if day_range[0] < 1:
        day_range = (1, day_range[1])
        warnings.warn("day_range[0] has defaulted to 1 (the first of the month).",
                      category=Warning)
    if day_range[1] > 31:
        day_range = (day_range[0], 31)
        warnings.warn("day_range[1] has defaulted to 31 (the assumed last day of the month).",
                      category=Warning)
    return day_range


def check_hour_range(hour_range):
    """
    :param hour_range: (<int>, <int>) or None: The hour range to check
    :return: (<int>, <int>): A suitable hour range
    """
    if hour_range is None:
        hour_range = (0, 24)  # Assume we want all hours
        warnings.warn("No hour range was provided, it has defaulted to (0, 24).",
                      category=Warning)
    if hour_range[0] < 0:
        hour_range = (0, hour_range[1])
        warnings.warn("hour_range[0] has defaulted to 0 (the very start of the day).",
                      category=Warning)
    if hour_range[1] > 24:
        hour_range = (hour_range[0], 24)
        warnings.warn("hour_range[1] has defaulted to 24 (the very end of the day).",
                      category=Warning)
    return hour_range


def check_gate_range(gate_range, hdw_info):
    """
    :param gate_range: (<int>, <int>) or None: The gate range to check
    :param hdw_info: _HdwInfo: A hardware info object
    :return: (<int>, <int>): A suitable gate range
    """
    if gate_range is None:
        gate_range = (0, 99)  # This is 100 gates
        warnings.warn("No gate range was provided, it has defaulted to (0, 99) - that is the first 100 gates.",
                      category=Warning)
    if gate_range[0] < 0:
        gate_range = (0, gate_range[1])
        warnings.warn("gate_range[0] has defaulted to 0 (the very first gate).",
                      category=Warning)
    if gate_range[1] > hdw_info.gates:
        gate_range = (gate_range[0], hdw_info.gates - 1)
        warnings.warn("The hardware file suggests there are only " + str(hdw_info.gates) +
                      " gates, so gate_range[1] has defaulted to " + str(gate_range[1]),
                      category=Warning)
    return gate_range


def check_beam_range(beam_range, hdw_info):
    """
    :param beam_range: (<int>, <int>) or None: The beam range to check
    :param hdw_info: _HdwInfo: A hardware info object
    :return: (<int>, <int>): A suitable beam range
    """
    if beam_range is None:
        beam_range = (0, 15)  # This is 16 beams
        warnings.warn("No beam range was provided, it has defaulted to (0, 15) - that is the first 16 beams.",
                      category=Warning)
    if beam_range[0] < 0:
        beam_range = (0, beam_range[1])
        warnings.warn("beam_range[0] has defaulted to 0 (the very first gate).",
                      category=Warning)
    if beam_range[1] > hdw_info.beams:
        beam_range = (beam_range[0], hdw_info.beams - 1)
        warnings.warn("The hardware file suggests there are only " + str(hdw_info.beams) +
                      " beams, so beam_range[1] has defaulted to " + str(beam_range[1]),
                      category=Warning)
    return beam_range


def check_time_units(time_units):
    """
    :param time_units: str: 'ut' for universal time or 'mlt' for magnetic local time
    :return: str: suitable time units
    """

    if time_units is None:
        time_units = "mlt"
        warnings.warn("No time_units provided, they have defaulted to mlt.",
                      category=Warning)

    time_units = time_units.lower()
    if time_units != "mlt" and time_units != "ut":
        warnings.warn("The provided time units were neither 'mlt' or 'ut', they have defaulted to mlt.",
                      category=Warning)
        time_units = "mlt"

    return time_units


