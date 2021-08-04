

def only_keep_overlap(station, df, other_station, gate_min=30):
    """

    Given a SuperDARN dataframe (df) with fit data for station, remove all data for cells that don't overlap with the
    other stations FoV.  The original dataframe is modified and also returned.

    :param station: str:
            The station whose df we are restricting, as 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param df: pandas.DataFrame
            station's dataframe.  Must contain 'gate' and 'beam' columns
    :param other_station: str
            The other station, also as a 3 character string.
    :param gate_min:
            Often the near gates will not see F region echoes.
            So, restrict both stations valid region of overlap to be only gates greater than this gate_min

    :return restricted_df:
            The provided dataframe, expect restricted to only those cells for which there is valid overlap with
            other_station
    """

    restricted_df = df.copy()
    return restricted_df
