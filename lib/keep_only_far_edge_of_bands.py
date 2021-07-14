

def keep_only_far_edge_of_bands(df):
    """

    Echoes commonly exist in bands.  This function builds a restricted dataframe containing only echoes near the far
    ends of the band.

    These echoes near the far edge of the bands might be important to isolate for determining elevation angle correction
    based on percent difference approach

    :param df: pandas.Date_Frame:
        A SuperDARD dataframe.
    :return: pandas.Date_Frame:
        The provided dataframe, except now with only far-edge-of-band echoes.
    """

    # TODO: Figure out if this is required, and if so what is the best way to do this
    #  Note: might need to bin the data?  Otherwise we might miss less dense bands?

    return df

