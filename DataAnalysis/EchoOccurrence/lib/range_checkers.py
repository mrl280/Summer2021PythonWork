def check_month_range(month_range):
    """
    :param month_range: (<int>, <int>) or None: The month range to check
    :return: (<int>, <int>): A suitable month range
    """
    if month_range is None:
        month_range = (1, 12)  # Assume we want all months
    if month_range[0] < 1:
        month_range = (1, month_range[1])
    if month_range[1] > 12:
        month_range = (month_range[0], 12)
    return month_range


def check_day_range(day_range):
    """
    :param day_range: (<int>, <int>) or None: The day range to check
    :return: (<int>, <int>): A suitable day range
    """
    if day_range is None:
        day_range = (1, 31)  # Assume we want all days
    if day_range[0] < 1:
        day_range = (1, day_range[1])
    if day_range[1] > 31:
        day_range = (day_range[0], 31)
    return day_range


def check_hour_range(hour_range):
    """
    :param hour_range: (<int>, <int>) or None: The hour range to check
    :return: (<int>, <int>): A suitable hour range
    """
    if hour_range is None:
        hour_range = (0, 24)  # Assume we want all hours
    if hour_range[0] < 0:
        hour_range = (0, hour_range[1])
    if hour_range[1] > 24:
        hour_range = (hour_range[0], 24)
    return hour_range


def check_gate_range(gate_range, hdw_info):
    """
    :param gate_range: (<int>, <int>) or None: The gate range to check
    :param hdw_info: _HdwInfo: A hardware info object
    :return: (<int>, <int>): A suitable gate range
    """
    if gate_range is None:
        gate_range = (0, 99)  # This is 100 gates
    if gate_range[0] < 0:
        gate_range = (0, gate_range[1])
    if gate_range[1] > hdw_info.gates:
        gate_range = (gate_range[0], hdw_info.gates - 1)
    return gate_range


def check_beam_range(beam_range, hdw_info):
    """
    :param beam_range: (<int>, <int>) or None: The beam range to check
    :param hdw_info: _HdwInfo: A hardware info object
    :return: (<int>, <int>): A suitable beam range
    """
    if beam_range is None:
        beam_range = (0, 15)  # This is 16 beams
    if beam_range[0] < 0:
        beam_range = (0, beam_range[1])
    if beam_range[1] > hdw_info.beams:
        beam_range = (beam_range[0], hdw_info.beams - 1)
    return beam_range
