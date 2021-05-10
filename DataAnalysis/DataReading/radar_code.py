

def radar_code(id):
    """
    :param: int id: The station id
    :return: String: 3 character long radar code
    """
    # SuperDARN Radars
    if id == 64:
        return 'inv'
    elif id == 65:
        return 'rkn'
    elif id == 66:
        return 'cly'
    elif id == 5:
        return 'sas'

    # RISR Incoherent Scatter Radars
    elif id == 91:
        return 'ran'
    elif id == 92:
        return 'ras'

    else:
        raise Exception("Error: radar_code_from_station_id() doesn't recognize station id " + str(id)
                        + ".  See https://madrigal.phys.ucalgary.ca/instMetadata for more radar codes.")

