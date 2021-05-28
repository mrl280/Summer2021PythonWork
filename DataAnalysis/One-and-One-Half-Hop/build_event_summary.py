import pandas as pd

import pathlib

if __name__ == '__main__':
    """
    This script builds a dataframe summarizing the relevant events
    The dataframe is pickled, so it can later be used to loop through the relevant events
    
    If you only want to process a sub-set of the events in the database, use a flag
    """

    file_name = "event_summary"

    column_names = ["SD_station",     # 3 char radar identifier  e.g. "rkn"
                    "RISR_station",   # 3 char radar identifier  e.g. "ran"  # Note: RISR-N = "ran"
                    "year",
                    "month",
                    "day",
                    "start_hour_UT",  # int
                    "end_hour_UT",    # int
                    "quality_flag"
                    # I am sure we will track more flags in the future
                    ]
    df = pd.DataFrame(columns=column_names)

    """ # Insert events below # """

    # Notes:
    df = df.append({'SD_station': "rkn",
                    'RISR_station': "ran",
                    'year': 2011,
                    'month': 11,
                    'day': 11,
                    'start_hour_UT': 16,
                    'end_hour_UT': 22,
                    'quality_flag': 1,
                    }, ignore_index=True)

    # Notes:
    df = df.append({'SD_station': "rkn",
                    'RISR_station': "ran",
                    'year': 2011,
                    'month': 11,
                    'day': 12,
                    'start_hour_UT': 18,
                    'end_hour_UT': 21,
                    'quality_flag': 1
                    }, ignore_index=True)

    # Notes:
    df = df.append({'SD_station': "rkn",
                    'RISR_station': "ran",
                    'year': 2011,
                    'month': 11,
                    'day': 14,
                    'start_hour_UT': 17,
                    'end_hour_UT': 21,
                    'quality_flag': 1
                    }, ignore_index=True)

    # Notes:
    df = df.append({'SD_station': "rkn",
                    'RISR_station': "ran",
                    'year': 2012,
                    'month': 10,
                    'day': 15,
                    'start_hour_UT': 18,
                    'end_hour_UT': 20,
                    'quality_flag': 1
                    }, ignore_index=True)


    # Save the dataframe to file
    loc_root = str(pathlib.Path().absolute())
    out_dir = loc_root + "/VelocityAnalysis/data"
    out_file = out_dir + "/" + file_name + ".pkl"
    print("Saving event summary as " + out_file)
    df.to_pickle(out_file)

