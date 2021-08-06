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

    # 2016 09 25 at RKN: gg[40, 74]
    # 2016 09 25 at CLY: gg[30, 74]
    # 2016 09 26 at CLY: gg[10, 74]
    # 2016 09 25 at INV: gg[10, 74]

    file_name = "event_summary"

    column_names = ["station",     # 3 char radar identifier  e.g. "rkn"
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

    # Notes: One of our main events
    df = df.append({'station': "rkn",
                    'year': 2016,
                    'month': 9,
                    'day': 26,
                    'start_hour_UT': 1,
                    'end_hour_UT': 4,
                    't_diff': 0.003,
                    'for_vel_hist': 1,
                    }, ignore_index=True)

    # Notes: One of our main events
    df = df.append({'station': "rkn",
                    'year': 2017,
                    'month': 10,
                    'day': 23,
                    'start_hour_UT': 4,
                    'end_hour_UT': 8,
                    't_diff': 0.003,
                    'for_vel_hist': 1,
                    }, ignore_index=True)

    # Notes: Some frequencies measure positive velocities while other measure negative
    df = df.append({'station': "rkn",
                    'year': 2017,
                    'month': 10,
                    'day': 23,
                    'start_hour_UT': 1,
                    'end_hour_UT': 4,
                    't_diff': 0.003,
                    'for_vel_hist': 0,
                    }, ignore_index=True)

    # Notes:
    df = df.append({'station': "rkn",
                    'year': 2016,
                    'month': 9,
                    'day': 5,
                    'start_hour_UT': 1,
                    'end_hour_UT': 8,
                    't_diff': 0,
                    'for_vel_hist': 1,
                    }, ignore_index=True)

    # Notes:
    df = df.append({'station': "rkn",
                    'year': 2017,
                    'month': 2,
                    'day': 4,
                    'start_hour_UT': 4,
                    'end_hour_UT': 7,
                    't_diff': 0.003,
                    'for_vel_hist': 1,
                    }, ignore_index=True)

    # Notes:
    df = df.append({'station': "rkn",
                    'year': 2016,
                    'month': 11,
                    'day': 5,
                    'start_hour_UT': 0,
                    'end_hour_UT': 8,
                    't_diff': 0.003,
                    'for_vel_hist': 1,
                    }, ignore_index=True)


    # Save the dataframe to file
    loc_root = str(pathlib.Path().absolute())
    out_dir = loc_root + "/data"
    out_file = out_dir + "/" + file_name + ".pbz2"
    print("Saving event summary as " + out_file)
    print(df)

    with bz2.BZ2File(out_file, "w") as file:
        pickle.dump(df, file)

