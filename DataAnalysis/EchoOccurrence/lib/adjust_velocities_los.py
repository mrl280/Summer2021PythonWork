import numpy as np
import pydarn
from pydarn import SuperDARNRadars


def adjust_velocities_los(station, df, ref_beam):
    """

    We need to modify LoS velocities so that line-of-sight velocity measurements from different directions can
        be compared directly - to do this select a reference beam, and we will adjust all velocity measurements so they
        are all as if looking along this same reference beam (same line-of-sight).

    :param station: str:
        The first radar station to consider, as 3 character string (e.g. "rkn").
    :param df: pandas.DataFrame:
        A SuperDARN fit dataframe with data for station.  Must include velocity ('v') column.
    :param ref_beam: int:
        We will adjust all velocity measurements so they are as if they were looking along this reference beam.

        There is nothing stopping you from choosing reference beams (including negative beams) outside of the actual
         fan - it will just adjust the number of bmsep's used in the adjustment

    :return df: pandas.DataFrame:
        The same dataframe, except not with the velocities 'v' adjusted as if they were looking along new_beam
    """

    all_radars_info = SuperDARNRadars()
    station_stid = pydarn.read_hdw_file(station).stid
    radar_info = all_radars_info.radars[station_stid]

    bmsep = radar_info.hardware_info.beam_separation  # Angular separation in degrees between adjacent beams [deg]
    # bmsep might be negative but that should be okay

    # Compute the angular separation in degrees between every beam and the ref_beam
    separation_from_ref_beam = bmsep * (df['bmnum'] - ref_beam)

    df['v'] = np.cos(np.radians(separation_from_ref_beam)) * df['v']

    return df

