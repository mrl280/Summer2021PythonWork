import numpy as np
import pandas as pd


def boxcar_smooth(raw_data, window_size=10):
    """

    Use a boxcar filter to smooth the data

    :param raw_data: numpy.array:
            The raw data that you would like to smooth
    :param window_size: int: (default is 10)
            The size of the rolling window
    :return: numpy.array:
            The smoothed data
    """

    # We already know how to do this with Pandas
    df = pd.DataFrame({'raw_data': raw_data})

    smoothed_data = df['raw_data'].rolling(window=window_size).sum() / window_size
    smoothed_data_shifted = np.roll(smoothed_data, -1 * int(window_size / 2))

    return smoothed_data_shifted
