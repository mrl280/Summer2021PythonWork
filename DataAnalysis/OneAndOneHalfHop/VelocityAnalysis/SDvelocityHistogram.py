import pandas as pd

import pathlib
from DataAnalysis.OneAndOneHalfHop.get_df_multi_event import get_df_multi_event

def SDvelocityHistogram(file_name, flag):
    """

    :param file_name: str: The name of the pickled data frame you would like to use.  e.g. "event_summary"
    :param flag: str: All events in the event summary that have this flag will be considered
    :return:
    """

    # Check to make sure inputs are okay
    if flag is None or file_name is None:
        raise Exception("Error: you must pass SDvelocityHistogram() a file_name and flag "
                        "so it knows which events to consider.")
    if not isinstance(flag, str) or not isinstance(file_name, str):
        raise Exception("Error: both file_name and flag must be strings.")

    df = get_df_multi_event(file_name, flag)

    print(df)


if __name__ == '__main__':
    """
    Create a velocity histogram
    """

    # Create a velocity histogram for a set of flagged events
    SDvelocityHistogram(file_name="event_summary", flag='for_vel_hist')
