import pandas as pd
import numpy as np
import math


def elevation_v2(df, t_diff):
    """
    This function is based on the rst code found here:
        https://github.com/SuperDARN/rst/blob/master/codebase/superdarn/src.lib/tk/elevation.1.0/src/elevation_v2.c
    The method is based on 'Elevation angle determination for SuperDARN HF radar layouts' by S. G. Shepherd

    :param df: SuperDARN data frame.  Must contain the following standard parameters: 'tfreq', 'phase'
    :param t_diff: The extra time delay to add in, in microseconds.
    :return: The input dataframe, but with the added parameters:
    """

    c = 2.998e+8  # Speed of light in vacuum

    if df['stationId'].iloc[0] == "rkn":
        x = 0  # From hdw.  Interferometer offset in metres
        y = 0
        z = -100
        maxbeam = 16  # From hdw.  Maximum number of beams to be used at a particular radar site.
        bmsep = 3.4  # From hdw.  Beam separation (Angular separation in degrees between adjacent beams).
        t_diff_hdw = 0.0  # From hdw.  Propagation time difference between the main array and the interferometer array.
    else:
        raise Exception("Error in elevation_v2(), " + df['stationId'].iloc[0] + " station not recognized")

    if y < 0:
        # Then the interferometer array is behind the main antenna
        y_sign = -1
    else:
        # The interferometer array is in front of the main antenna
        y_sign = 1

    # Find the total t_diff.  That from the hardware file plus the extra we are adding in
    t_diff += t_diff_hdw

    # Find distance from the main array to the interferometer array
    d = math.sqrt(x * x + y * y + z * z)

    # Find beam offset
    boff = maxbeam / 2.0 - 0.5  # middle of the beams?
    phi0 = math.radians(bmsep * (maxbeam - boff))
    cp0 = math.cos(phi0)
    sp0 = math.sin(phi0)

    # Find the phase delay due to the electrical path difference
    psi_ele = -2 * math.pi * df['transFreq'] * t_diff * 1e-3  # TODO: Not sure about the scale

    # Find the elevation angle where phase difference is max
    # For most sites: z=0 and so a0=0
    a0 = math.asin(y_sign * z * cp0 / math.sqrt(y * y + z * z))

    if a0 < 0:
        # Assume that negative elevation angles are unphysical
        a0 = 0

    ca0 = math.cos(a0)
    sa0 = math.sin(a0)

    # Find the maximum phase.  This is phi_ele + phi_geo(a0)
    psi_max = psi_ele + \
                    2 * math.pi * df['transFreq'] * 1e3 / c * (x * sp0 + y * math.sqrt(ca0 * ca0 - sp0 * sp0) + z * sa0)

    # "Unwrap" phase
    # Find the number of 2pi factors required to map phase to the correct region
    dpsi = psi_max - df['phase']  # Max phase minus observed phase
    if y > 0:
        n2pi = np.floor(dpsi / (2 * math.pi))
    else:
        n2pi = np.ceil(dpsi / (2 * math.pi))
    d2pi = n2pi * 2 * math.pi

    # map the observed phase to the correct region
    df['adjPhase'] = df['phase'] + d2pi

    # And we can solve for elevation angle alpha
    E = (df['adjPhase'] / (2 * math.pi * df['transFreq'] * 1e3) + t_diff * 1e-6) * c - x * sp0
    df['adjElv'] = np.arcsin((E * z + np.sqrt(E * E * z * z - (y * y + z * z) * (E * E - y * y * cp0 * cp0))) /
                             (y * y + z * z))


if __name__ == '__main__':
    """
    Read in and look at pickled SuperDARN data
    Doesn't produce anything, just for looking at the file structure
    """

    station = "rkn"
    date = "20160925"

    in_dir = "data/" + station + "/" + station + date
    in_file = in_dir + "/" + station + date + ".pkl"
    df = pd.read_pickle(in_file)

    print(df.keys())
    elevation_v2(df, 0)
    print("After call:")
    print(df.keys())
