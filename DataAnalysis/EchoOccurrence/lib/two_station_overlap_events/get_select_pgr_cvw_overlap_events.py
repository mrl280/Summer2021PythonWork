import pandas as pd


def get_select_prg_cvw_overlap_events():
    """"
    Build and return an event dataframe to hold all the two radar overlap events we want to look at
    """

    station1 = "pgr"
    gate_range1 = (20, 74)
    beam_range1 = (5, 6)
    station1_ref_beam = 6

    station2 = "cvw"
    gate_range2 = (20, 74)
    beam_range2 = (15, 15)
    station2_ref_beam = 15

    freq_range = (8, 11)
    plot_type = 'contour'

    column_names = ["station1",  # 3 char radar identifier
                    "station2",  # 3 char radar identifier
                    "start_epoch",
                    "end_epoch",
                    "station1_ref_beam", "station2_ref_beam",
                    "gate_range1", "beam_range1", "gate_range2", "beam_range2", "freq_range",
                    "plot_type"
                    ]
    event_df = pd.DataFrame(columns=column_names)

    # # Notes:
    # event_df = event_df.append({"station1": station1,
    #                             "station2": station2,
    #                             "start_epoch": 1570647600,  # 2019-10-09 19:00
    #                             "end_epoch": 1570662000,  # 2019-10-09 23:00
    #                             "station1_ref_beam": station1_ref_beam, "station2_ref_beam": station2_ref_beam,
    #                             "gate_range1": gate_range1, "beam_range1": beam_range1,
    #                             "gate_range2": gate_range2, "beam_range2": beam_range2, "freq_range": freq_range,
    #                             "plot_type": plot_type
    #                             }, ignore_index=True)

    return event_df