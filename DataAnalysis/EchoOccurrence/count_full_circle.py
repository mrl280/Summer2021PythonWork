import pathlib
import pydarn

import cartopy.crs as ccrs
import matplotlib.path as mpath
import matplotlib.ticker as mticker

from matplotlib import pyplot as plt
from pydarn import SuperDARNRadars, radar_fov
from scipy import stats

from lib.build_datetime_epoch import build_datetime_epoch
from lib.cm.modified_viridis import modified_viridis
from lib.add_mlt_to_df import add_mlt_to_df
from lib.only_keep_45km_res_data import only_keep_45km_res_data
from lib.get_data_handler import get_data_handler
from lib.z_min_max_defaults import z_min_max_defaults
from lib.data_getters.input_checkers import check_year


def occ_full_circle(station, year, time_units='mlt', month_range=None, day_range=None, hour_range=None, gate_range=None, beam_range=None,
                    local_testing=False, parameter=None, plot_ground_scat=False):
    """

    TODO: There are still concerns around what to do here.. not all cells cover the same area,
      and the lower latitude regions are scanned much less often then an the higher lattitude points
    TODO: Also the data is in aagmc but the plot and gridlines are in geographic coordinates - which doesn't make sense

    Produce a full circle stereographic plot in either ut or mlt.
    Can plot a simple echo count, ground scatter count, or average a fitACF parameter over the provided time range.

    Notes:
        - This program was originally written to be run on maxwell.usask.ca.  This decision was made because
            occurrence investigations often require chewing large amounts of data.
        - Does not distinguish frequency
        - Only considers 45 km data.
            (a warning will be printed if other spatial resolution data is stripped from the dataset)
        - This program uses fitACF 3.0 data.  To change this, modify the source code.
        - All times and dates are assumed UT

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param year: int:
            The year to consider.
    :param time_units: str: 'ut' for universal time or 'mlt' for magnetic local time:
            The time units to plot on the circle, 12 is always at the top.  Default is 'mlt'
    :param month_range: (<int>, <int>) (optional):
            Inclusive. The months of the year to consider.  If omitted (or None), then all days will be considered.
    :param day_range: (<int>, <int>) (optional):
            Inclusive. The days of the month to consider.  If omitted (or None), then all days will be considered.
    :param hour_range: (<int>, <int>) (optional):
            The hour range to consider.  If omitted (or None), then all hours will be considered.
            Not inclusive: if you pass in (0, 5) you will get from 0:00-4:59 UT
    :param gate_range: (<int>, <int>) (optional):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (<int>, <int>) (optional):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param local_testing: bool (optional):
            Set this to true if you are testing on your local machine.  Program will then use local dummy data.
    :param parameter: str (optional):
            Parameter to be averaged (e.g. 'v' or 'p_l')
            If omitted, then a simple echo count will be plotted.
    :param plot_ground_scat: bool (optional)
            Set this to true if you would like to plot ground scatter counts.  Default is False
            If plot_ground_scat is set to True, then :param parameter is ignored.
    :return: matplotlib.pyplot.figure: The figure.
            It can then be modified, added to, printed out, or saved in whichever file format is desired.
    """

    year = check_year(year)

    print("Retrieving data...")
    df = get_data_handler(station, year_range=(year, year), month_range=month_range, day_range=day_range,
                          hour_range=hour_range, gate_range=gate_range, beam_range=beam_range,
                          local_testing=local_testing)

    print("Getting some hardware info...")
    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    hemisphere = this_radars_info.hemisphere
    radar_id = this_radars_info.hardware_info.stid

    print("Filtering data...")
    df = df.loc[(df['p_l'] >= 3)]  # Restrict to points with at least 3 dB

    if not plot_ground_scat and parameter is not None:
        zmin, zmax = z_min_max_defaults(parameter)
        df = df.loc[(df[parameter] >= zmin) & (df[parameter] <= zmax)]

    df.reset_index(drop=True, inplace=True)
    df = only_keep_45km_res_data(df)

    print("Preparing the plot...")
    fig = plt.figure(figsize=(5, 5), dpi=300)
    if hemisphere.value == 1:
        # Northern hemisphere
        min_lat = 37  # deg
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, 90, min_lat], crs=ccrs.PlateCarree())
    elif hemisphere.value == -1:
        # Southern hemisphere
        max_lat = -34  # deg
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90, max_lat], crs=ccrs.PlateCarree())
    else:
        raise Exception("occ_full_circle(): hemisphere not recognized")

    # Compute a circle in axis coordinates which can be used as a boundary
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    vertices = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(vertices * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    # Add gridlines and mlt labels
    text_offset_multiplier = 1.03
    gl = ax.gridlines(draw_labels=False)
    gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135])
    ax.text(0, text_offset_multiplier * ax.get_ylim()[1], "12", ha='center', va='bottom')
    ax.text(0, text_offset_multiplier * ax.get_ylim()[0], "00", ha='center', va='top')
    ax.text(text_offset_multiplier * ax.get_xlim()[1], 0, "06", ha='left', va='center')
    ax.text(text_offset_multiplier * ax.get_xlim()[0], 0, "18", ha='right', va='center')

    # TODO: Only compute mlt if needed, make ut option
    print("Computing MLT...")
    # To compute mlt we need longitudes..
    # we will use the start of the year and assume magnetic longitudes don't change much over the observation period
    start_datetime, start_epoch = build_datetime_epoch(year=year, month=1, day=1, hour=0)
    cell_corners_aacgm_lats, cell_corners_aacgm_lons = radar_fov(stid=radar_id, coords='aacgm', date=start_datetime)

    df = add_mlt_to_df(cell_corners_aacgm_lons=cell_corners_aacgm_lons,
                          cell_corners_aacgm_lats=cell_corners_aacgm_lats, df=df)

    # Right now MLT is in the range 0-24, we need to put it in the range 0-360 for circular plotting
    df['mlt'] = 15 * df['mlt']

    n_bins_x = 180  # 1 bin per every 2 degrees
    n_bins_y = 40  # About 1 bin per degree
    # TODO: This depends on hemisphere, need to make compatible with southern hemisphere
    contour_range = [[0, 360], [37, 90]]
    if parameter is None:
        # We just want a simple echo count
        binned_counts, bin_xedges, bin_yedges, bin_numbers = stats.binned_statistic_2d(
            df['mlt'], df['lat'], values=None,
            statistic='count', bins=[n_bins_x, n_bins_y], range=contour_range)
    else:
        # We want to the median of the chosen parameter
        binned_counts, bin_xedges, bin_yedges, bin_numbers = stats.binned_statistic_2d(
            df['mlt'], df['lat'], values=df[parameter],
            statistic='median', bins=[n_bins_x, n_bins_y], range=contour_range)
        # TODO: Make so no data cells plot like nothing is there
        binned_counts = np.nan_to_num(binned_counts)

    # Compute bin centers
    bin_xwidth = (bin_xedges[1] - bin_xedges[0])
    bin_ywidth = (bin_yedges[1] - bin_yedges[0])
    bin_xcenters = bin_xedges[1:] - bin_xwidth / 2
    bin_ycenters = bin_yedges[1:] - bin_ywidth / 2

    # Plot the data
    if plot_ground_scat:
        raise Exception("Ground scatter not yet supported!")  # TODO: Ground scatter not currently supported
    else:
        if parameter == 'v':
            cmap = 'seismic_r'
        else:
            cmap = modified_viridis()
        cont = ax.contourf(bin_xcenters, bin_ycenters, binned_counts.transpose(),
                           levels=6, cmap=cmap, transform=ccrs.PlateCarree())

    # ax.plot([radar_geo_lon, radar_geo_lon], [radar__geo_lat, radar__geo_lat], 'ro', transform=ccrs.PlateCarree())
    # cont = ax.contourf(bin_xcenters, bin_ycenters, binned_counts.transpose(), 5, transform=ccrs.PlateCarree())
    fig.colorbar(cont, fraction=0.046, orientation="horizontal", ax=ax)

    return fig


if __name__ == '__main__':
    """ Testing """

    local_testing = True
    station = "rkn"

    fig = occ_full_circle(station=station, year=2011, month_range=(2, 2), day_range=None,
                          gate_range=(0, 74), beam_range=(0, 15),
                          plot_ground_scat=False, parameter=None, local_testing=local_testing)

    if local_testing:
        plt.show()
    else:
        loc_root = str((pathlib.Path().parent.absolute()))
        out_dir = loc_root + "/out"
        out_file = out_dir + "/count_full_circle_" + station
        print("Saving plot as " + out_file)
        fig.savefig(out_file + ".jpg", format='jpg', dpi=300)
