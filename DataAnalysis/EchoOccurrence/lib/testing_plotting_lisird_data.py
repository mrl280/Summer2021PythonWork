import os

import numpy as np
from matplotlib import pyplot as plt

from matplotlib.ticker import MultipleLocator

from DataAnalysis.EchoOccurrence.lib.data_getters.get_LISIRD_data import get_LISIRD_data


def plot_LISIRD_data(dataset_name):
    """

    TODO: Info about the solar flux data can be found here:
     https://www.spaceweather.gc.ca/forecast-prevision/solar-solaire/solarflux/sx-3-en.php

    Here is a plot of what the data should look like:
     https://lasp.colorado.edu/lisird/data/penticton_radio_flux/

    :param dataset_name:
    :return:
    """
    df = get_LISIRD_data(dataset_name=dataset_name)

    print("Computing UT decimal years/hours...")
    minutes_in_an_hour = 60
    seconds_in_an_hour = 3600
    months_in_a_year = 12
    days_in_a_year = 365
    hours_in_a_year = 8760

    # Compute decimal years to plot along the x axis
    decimal_hours, decimal_years = [], []
    for i in range(len(df)):
        datetime_obj_here = df['datetime'].iat[i]

        decimal_hour_here = datetime_obj_here.hour + datetime_obj_here.minute / minutes_in_an_hour + \
                            datetime_obj_here.second / seconds_in_an_hour

        decimal_year_here = datetime_obj_here.year + (datetime_obj_here.month - 1) / months_in_a_year + \
                            (datetime_obj_here.day - 1) / days_in_a_year + decimal_hour_here / hours_in_a_year

        decimal_hours.append(decimal_hour_here)
        decimal_years.append(decimal_year_here)

    df['decimal_hour'] = np.asarray(decimal_hours)
    df['decimal_year'] = np.asarray(decimal_years)

    fig, ax = plt.subplots(figsize=[8, 6], dpi=300, constrained_layout=True, nrows=1, ncols=1)
    fig.suptitle(dataset_name + "\nProduced by " + str(os.path.basename(__file__)), fontsize=18)

    ax.set_xlim([2005, 2021])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    # ax.xaxis.set_minor_locator(MultipleLocator(0.25))

    y_lim = [0, 200]
    ax.set_ylim(y_lim)
    # ax.yaxis.set_major_locator(mticker.FixedLocator(y_axis_major_labels))
    # ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5, color='black')
    ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2, color='black')

    ax.set_xlabel("Year", fontsize=14)
    ax.set_ylabel("Solar Flux Units (SFU)", fontsize=14)

    ax.plot(df['decimal_year'], df['observed_flux'], color="cornflowerblue", linestyle='-', linewidth=0.5)

    # Smooth the data by applying a boxcar filter
    smoothing_window = 120
    smoothed_data = df['observed_flux'].rolling(window=smoothing_window).sum() / smoothing_window
    smoothed_data_shifted = np.roll(smoothed_data, -1 * int(smoothing_window / 2))

    ax.plot(df['decimal_year'], smoothed_data_shifted, 'b-')

    return fig


if __name__ == "__main__":
    """ Testing """

    dataset_1 = "penticton_radio_flux_nearest_noon"

    fig = plot_LISIRD_data(dataset_name=dataset_1)

    plt.show()
    plt.close(fig)

