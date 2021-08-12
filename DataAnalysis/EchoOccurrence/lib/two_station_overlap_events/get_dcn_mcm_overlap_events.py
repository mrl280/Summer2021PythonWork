import pandas as pd


def get_dcn_mcm_overlap_events():
    """"
    Build and return an event dataframe to hold all the two radar overlap events we want to look at
    """

    station1 = "dcn"
    station1_ref_beam = 15

    station2 = "mcm"
    station2_ref_beam = 8

    gate_range = (20, 74)
    beam_range = (0, 15)
    freq_range = (8, 11)
    plot_type = 'contour'

    column_names = ["station1",  # 3 char radar identifier
                    "station2",  # 3 char radar identifier
                    "start_epoch",
                    "end_epoch",
                    ]
    event_df = pd.DataFrame(columns=column_names)

    # Notes:
    event_df = event_df.append({"station1": station1,
                                "station2": station2,
                                "start_epoch": 1570647600,  # 2019-10-09 19:00
                                "end_epoch": 1570662000,  # 2019-10-09 23:00
                                "station1_ref_beam": station1_ref_beam,     "station2_ref_beam": station2_ref_beam,
                                "gate_range": gate_range,   "beam_range": beam_range,   "freq_range": freq_range,
                                "plot_type": plot_type
                                }, ignore_index=True)

    # Notes:
    event_df = event_df.append({"station1": station1,
                                "station2": station2,
                                "start_epoch": 1572463800,  # 2019-10-30 19:30
                                "end_epoch": 1572478200,  # 2019-10-30 23:30
                                "station1_ref_beam": station1_ref_beam, "station2_ref_beam": station2_ref_beam,
                                "gate_range": gate_range, "beam_range": beam_range, "freq_range": freq_range,
                                "plot_type": plot_type
                                }, ignore_index=True)

    # Notes:
    event_df = event_df.append({"station1": station1,
                                "station2": station2,
                                "start_epoch": 1551369600,  # 2019-02-28 16:00
                                "end_epoch": 1551384000,  # 2019-02-28 20:00
                                "station1_ref_beam": station1_ref_beam, "station2_ref_beam": station2_ref_beam,
                                "gate_range": gate_range, "beam_range": beam_range, "freq_range": freq_range,
                                "plot_type": plot_type
                                }, ignore_index=True)

    # Notes:
    event_df = event_df.append({"station1": station1,
                                "station2": station2,
                                "start_epoch": 1567305000,  # 2019-09-01 2:30
                                "end_epoch": 1567319400,  # 2019-09-01 6:30
                                "station1_ref_beam": station1_ref_beam, "station2_ref_beam": station2_ref_beam,
                                "gate_range": gate_range, "beam_range": beam_range, "freq_range": freq_range,
                                "plot_type": plot_type
                                }, ignore_index=True)

    # Notes:
    event_df = event_df.append({"station1": station1,
                                "station2": station2,
                                "start_epoch": 1571943600,  # 2019-10-24 19:00
                                "end_epoch": 1571958000,  # 2019-10-24 23:00
                                "station1_ref_beam": station1_ref_beam, "station2_ref_beam": station2_ref_beam,
                                "gate_range": gate_range, "beam_range": beam_range, "freq_range": freq_range,
                                "plot_type": plot_type
                                }, ignore_index=True)

    # Notes:
    event_df = event_df.append({"station1": station1,
                                "station2": station2,
                                "start_epoch": 1549416600,  # 2019-02-06 1:30
                                "end_epoch": 1549431000,  # 2019-02-06 5:30
                                "station1_ref_beam": station1_ref_beam, "station2_ref_beam": station2_ref_beam,
                                "gate_range": gate_range, "beam_range": beam_range, "freq_range": freq_range,
                                "plot_type": plot_type
                                }, ignore_index=True)

    # Notes:
    event_df = event_df.append({"station1": station1,
                                "station2": station2,
                                "start_epoch": 1549416600,  # 2019-02-06 1:30
                                "end_epoch": 1549431000,  # 2019-02-06 5:30
                                "station1_ref_beam": station1_ref_beam, "station2_ref_beam": station2_ref_beam,
                                "gate_range": gate_range, "beam_range": beam_range, "freq_range": freq_range,
                                "plot_type": plot_type
                                }, ignore_index=True)

    # Notes:
    event_df = event_df.append({"station1": station1,
                                "station2": station2,
                                "start_epoch": 1573486200,  # 2019-11-11 15:30
                                "end_epoch": 1573500600,  # 2019-11-11 19:30
                                "station1_ref_beam": station1_ref_beam, "station2_ref_beam": station2_ref_beam,
                                "gate_range": gate_range, "beam_range": beam_range, "freq_range": freq_range,
                                "plot_type": plot_type
                                }, ignore_index=True)

    # Notes:
    event_df = event_df.append({"station1": station1,
                                "station2": station2,
                                "start_epoch": 1570761000,  # 2019-10-11 2:30
                                "end_epoch": 1570775400,  # 2019-10-11 6:30
                                "station1_ref_beam": station1_ref_beam, "station2_ref_beam": station2_ref_beam,
                                "gate_range": gate_range, "beam_range": beam_range, "freq_range": freq_range,
                                "plot_type": plot_type
                                }, ignore_index=True)

    # Notes:
    event_df = event_df.append({"station1": station1,
                                "station2": station2,
                                "start_epoch": 1573497000,  # 2019-11-11 18:30
                                "end_epoch": 1573511400,  # 2019-11-11 22:30
                                "station1_ref_beam": station1_ref_beam, "station2_ref_beam": station2_ref_beam,
                                "gate_range": gate_range, "beam_range": beam_range, "freq_range": freq_range,
                                "plot_type": plot_type
                                }, ignore_index=True)

    return event_df
