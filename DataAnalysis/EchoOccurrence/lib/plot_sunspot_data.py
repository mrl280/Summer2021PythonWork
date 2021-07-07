from matplotlib.ticker import MultipleLocator
import numpy as np

from .boxcar_smooth import boxcar_smooth
from .data_getters.get_LISIRD_data import get_LISIRD_data
from .data_getters.input_checkers import check_year_range


def plot_sunspot_data(ax, year_range, smoothing_window_size_in_days):
    """

    Plot sunspot on the provided axes.

    :param ax: matplotlib.pyplot.axis:
            The axis to draw on
    :param year_range: (int, int):
            Inclusive. The year range to consider. If omitted (or None), then all years will be considered.
    :param smoothing_window_size_in_days:
            This size of the window to use when boxcar smoothing the data
    :return: pandas.Date_Frame:
            A dataframe containing the solar flux data, it can then be further investigated as required.
    """

    sunspot_dataset = "international_sunspot_number"
    sm_in_days = smoothing_window_size_in_days  # Convenience

    # Format the axes
    ax.set_xlim(year_range)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))

    # y_lim = [50, 200]
    # ax.set_ylim(y_lim)
    # ax.yaxis.set_major_locator(mticker.FixedLocator(y_axis_major_labels))
    # ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    ax.grid(b=True, which='major', axis='both', linestyle='--', linewidth=0.5, color='black')
    ax.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.2, color='black')

    ax.set_xlabel("Year", fontsize=14)
    ax.set_ylabel("Sunspot Number [count]", fontsize=14)

    # Get the data
    df = get_LISIRD_data(dataset_name=sunspot_dataset)

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

    # Plot the raw data
    ax.plot(df['decimal_year'], df['sunspot_number'], color="cornflowerblue", linestyle='-', linewidth=0.5)

    # Smooth the data with a boxcar filter, and plot that on top of the raw data
    smoothed_data = boxcar_smooth(df['sunspot_number'], window_size=sm_in_days)
    ax.plot(df['decimal_year'], smoothed_data, color="blue", linewidth="1", linestyle="-")

    return df
