import bz2
import glob
import os
import pathlib
import warnings

import numpy as np
import pandas as pd

from scipy import stats
from matplotlib import pyplot as plt


def get_IMF_data(year_range):
    """
    Get OMNI IMF data

    Note: IMF data should be pickled with pickle_imf() before calling this.
     An error is raised if no pickled data is found.

    Please see pickle_imf() for a more complete description of IMF data and its origins

    :param year_range: (int, int):
            Inclusive. The year range to consider.

    :return: pandas.Date_Frame:
            A dataframe with the requested data.
    """

    loc_root = str((pathlib.Path().parent.absolute()))
    if "data_getters" in loc_root:
        # We are calling from within lib/data_getters
        loc_root = str((pathlib.Path().parent.absolute().parent.absolute().parent.absolute()))
    if "lib" in loc_root:
        # We are calling from within lib
        loc_root = str((pathlib.Path().parent.absolute().parent.absolute()))

    in_dir = loc_root + "/data/omni"
    print("Looking for a pickled file in: " + in_dir)

    for in_file in glob.iglob(in_dir + "/*.pbz2"):
        file_name = str(os.path.basename(in_file))

        try:
            # Try to pull the required information our of the filename
            file_year_start = int(file_name[9:13])
            file_year_end = int(file_name[14:18])

            # Lets see if this file will do the trick
            if file_year_start <= year_range[0] and file_year_end >= year_range[1]:
                # We are good, go ahead and read in the data

                warnings.warn("Using IMF data from " + in_file, category=Warning)
                data_stream = bz2.BZ2File(in_file, "rb")
                return pd.read_pickle(data_stream)

            else:
                continue  # This file is not what we are looking for
        except BaseException as e:
            # print(e)
            pass

    # If we didn't find anything, then we can just go ahead and return None
    raise FileNotFoundError("get_IMF_data() could not find a pickled IMF datafile with information from "
                            + str(year_range[0]) + " to " + str(year_range[1])
                            + ".  Please pickle data for the desired range using pickle_imf()")


if __name__ == "__main__":
    """ Testing """

    # Most of the IMF data in both by and bz in the range of -5 to 5 nT

    year_range = (2014, 2021)

    df = get_IMF_data(year_range=year_range)

    print(df.keys())

    print("")
    print("Minimum Bx value: " + str(np.min(df['Bx_nT'])) + " nT (GSM)")
    print("Maximum Bx value: " + str(np.max(df['Bx_nT'])) + " nT (GSM)")

    print("")
    print("Minimum By value: " + str(np.min(df['By_nT_GSM'])) + " nT (GSM)")
    print("Maximum By value: " + str(np.max(df['By_nT_GSM'])) + " nT (GSM)")

    print("")
    print("Minimum Bz value: " + str(np.min(df['Bz_nT_GSM'])) + " nT (GSM)")
    print("Maximum Bz value: " + str(np.max(df['Bz_nT_GSM'])) + " nT (GSM)")

    # Quickly plot the data so we can look at the spread
    fig, ax = plt.subplots(figsize=[6, 6], nrows=1, ncols=1, constrained_layout=True, dpi=300)

    # Compute By edges
    n_bins_x = 40
    by_edges = np.linspace(-10, 10, num=(n_bins_x + 1))

    # Compute By edges
    n_bins_y = 40
    bz_edges = np.linspace(-10, 10, num=(n_bins_y + 1))

    result, _, _, _ = stats.binned_statistic_2d(df['By_nT_GSM'], df['Bz_nT_GSM'], values=None, statistic='count',
                                                bins=[by_edges, bz_edges])

    plot = ax.pcolormesh(by_edges, bz_edges, result.transpose(), cmap='jet', zorder=0)
    cbar = fig.colorbar(plot, ax=ax, orientation="vertical")

    # ax.scatter(df['By_nT_GSM'], df['Bz_nT_GSM'])
    ax.set_xlabel("By [nT] (GSM)")
    ax.set_ylabel("Bz [nT] (GSM)")
    ax.set_title("IMF Data Spread")

    plt.show()
    plt.close(fig)
