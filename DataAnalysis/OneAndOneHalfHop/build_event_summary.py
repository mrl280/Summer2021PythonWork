import bz2
import pathlib

import pickle
import pandas as pd


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

                    # flags
                    "for_vel_hist"
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

                    'for_vel_hist': 1,
                    }, ignore_index=True)

    # Notes:
    df = df.append({'SD_station': "rkn",
                    'RISR_station': "ran",
                    'year': 2011,
                    'month': 11,
                    'day': 12,
                    'start_hour_UT': 18,
                    'end_hour_UT': 21,

                    'for_vel_hist': 1
                    }, ignore_index=True)

    # Notes:
    df = df.append({'SD_station': "rkn",
                    'RISR_station': "ran",
                    'year': 2011,
                    'month': 11,
                    'day': 14,
                    'start_hour_UT': 17,
                    'end_hour_UT': 21,

                    'for_vel_hist': 1
                    }, ignore_index=True)

    # Notes:
    df = df.append({'SD_station': "rkn",
                    'RISR_station': "ran",
                    'year': 2012,
                    'month': 10,
                    'day': 15,
                    'start_hour_UT': 18,
                    'end_hour_UT': 20,

                    'for_vel_hist': 0
                    }, ignore_index=True)


    # Save the dataframe to file
    loc_root = str(pathlib.Path().absolute())
    out_dir = loc_root + "/data"
    out_file = out_dir + "/" + file_name + ".pbz2"
    print("Saving event summary as " + out_file)
    with bz2.BZ2File(out_file, "w") as file:
        pickle.dump(df, file)

