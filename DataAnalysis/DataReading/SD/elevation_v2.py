import math
import pydarn

import pandas as pd
import numpy as np


def elevation_v2(df, t_diff=0.0):
    """
    This function is based on the rst code found here:
        https://github.com/SuperDARN/rst/blob/master/codebase/superdarn/src.lib/tk/elevation.1.0/src/elevation_v2.c
    The method is based on 'Elevation angle determination for SuperDARN HF radar layouts' by S. G. Shepherd

    Note:
    - This program is a tad slow, so it is best to filter your dataframe before calling it.
    - In the fitACF file, a phase offset of 0 results in an elevation angle of 0.  However, that is not the case
     for this program.  Rather, you will get the appropriate non-zero elevation angle.

    :param df: pandas.DataFrame:
        SuperDARN data frame.  Must contain the following standard parameters: 'transFreq', 'phase', 'station', 'bmnum'
    :param t_diff: float: (Optional, default is 0.0)
        The extra time delay to add in, in microseconds.
    :return: pandas.DataFrame:
        The input dataframe, but with the added parameters: 'adjPhase', 'adjElv'
    """

    df.reset_index(drop=True, inplace=True)  # Just encase
    if len(df) <= 0:
        df['adjPhase'] = []
        df['adjElv'] = []
        return

    c = 2.998e+8  # Speed of light in vacuum [m/s]

    # Pull the required information out of the hardware file
    station = df['station'].iat[0]
    hdw_info = pydarn.read_hdw_file(station)
    x = hdw_info.interferometer_offset.x  # Metres
    y = hdw_info.interferometer_offset.y  # Metres
    z = hdw_info.interferometer_offset.z  # Metres
    maxbeam = hdw_info.beams
    bmsep = hdw_info.beam_separation  # Angular separation in degrees between adjacent beams [deg]
    t_diff_hdw = hdw_info.tdiff  # Microseconds

    if y < 0:
        # Then the interferometer array is behind the main antenna
        y_sign = -1
    else:
        # The interferometer array is in front of the main antenna
        y_sign = 1

    # Find the total t_diff.  That from the hardware file plus the extra we are adding in
    t_diff += t_diff_hdw

    # Find beam offset
    boff = maxbeam / 2.0 - 0.5  # middle of the beams?
    phi0 = np.radians(bmsep * (df['bmnum'] - boff))
    cp0 = np.cos(phi0)
    sp0 = np.sin(phi0)

    # Find the phase delay due to the electrical path difference
    psi_ele = -2 * math.pi * df['transFreq'] * t_diff * 1e-3

    # Find the elevation angle where phase difference is max
    # For most sites: z=0 and so a0=0
    a0 = np.arcsin(y_sign * z * cp0 / math.sqrt(y * y + z * z))

    # Assume that negative elevation angles are unphysical
    for i in a0:
        if a0[i] < 0:
            a0[i] = 0

    # ca0 = np.cos(a0)
    sa0 = np.sin(a0)

    # Find the maximum phase.  This is phi_ele + phi_geo(a0)
    # TODO: In rst it is ca0 * ca0 - sp0 * sp0, I am following what is in the paper
    psi_max = psi_ele + \
              2 * math.pi * df['transFreq'] * 1e3 / c * (x * sp0 + y * np.sqrt(cp0 * cp0 - sa0 * sa0) + z * sa0)

    # "Unwrap" phase
    # Find the number of 2pi factors required to map phase to the correct region
    dpsi = psi_max - df['phase']  # Max phase minus observed phase
    if y > 0:
        n = np.floor(dpsi / (2 * math.pi))
    else:
        n = np.ceil(dpsi / (2 * math.pi))

    # map the observed phase to the correct region
    df['adjPhase'] = df['phase'] + (n * 2 * math.pi)

    # And we can solve for elevation angle alpha
    E = (df['adjPhase'] / (2 * math.pi * df['transFreq'] * 1e3) + t_diff * 1e-6) * c - x * sp0
    df['adjElv'] = np.degrees(
        np.arcsin((E * z + np.sqrt(E * E * z * z - (y * y + z * z) * (E * E - y * y * cp0 * cp0))) /
                  (y * y + z * z)))


if __name__ == '__main__':
    """
    Compare fitACF elevation angles to elevation angles computed by elevation_v2()
    """

    # station = "rkn"
    # date = "20111011"
    station = "rkn"
    date = "20160926"

    in_dir = "data/" + station + "/" + station + date
    in_file = in_dir + "/" + station + date + ".pkl"
    df = pd.read_pickle(in_file)

    t_diff = 0.0003  # Time delay in microseconds

    elevation_v2(df, t_diff)
    print("df keys:")
    print(df.keys())
    print(df[['transFreq', 'vel', 'phase', 'adjPhase', 'elv', 'adjElv']])
    print(1e-3)
