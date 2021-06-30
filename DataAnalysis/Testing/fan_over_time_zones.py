import pathlib
import warnings
import aacgmv2
import pydarn
import pytz

import cartopy.crs as ccrs
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import numpy as np
import cartopy.feature as cfeature
import datetime as datetime

from tzwhere import tzwhere
from aacgmv2 import convert_latlon_arr
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import radar_fov, SuperDARNRadars
from cartopy.io import shapereader as shpreader
from cartopy.feature import ShapelyFeature
import shapefile
from timezonefinder import TimezoneFinder

from DataAnalysis.EchoOccurrence.lib.add_mlt_to_df import centroid

loc_root = str((pathlib.Path().parent.absolute()))

all_radars_info = SuperDARNRadars()

# Note: Radars not in the same hemisphere will be omitted
# stations = ["dce", "dcn"]
stations = ["dce"]
reference_hemisphere = all_radars_info.radars[pydarn.read_hdw_file(stations[0]).stid].hemisphere

date = datetime.datetime.now()

print("Preparing the plot...")
fig = plt.figure(figsize=(10, 5), dpi=300, constrained_layout=True)
lat_extreme = 40 * reference_hemisphere.value  # deg
if reference_hemisphere.value == 1:
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
elif reference_hemisphere.value == -1:
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())
else:
    raise Exception("hemisphere not recognized")

ax.set_extent([-180, 180, reference_hemisphere.value * 90, lat_extreme], crs=ccrs.PlateCarree())

# Compute a circle in axis coordinates which can be used as a boundary
theta = np.linspace(0, 2 * np.pi, 100)
center, radius = [0.5, 0.5], 0.5
vertices = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(vertices * radius + center)
# ax.set_boundary(circle, transform=ax.transAxes)

# Add gridlines and mlt labels
text_offset_multiplier = 1.03
gl = ax.gridlines(draw_labels=True, linestyle='--', color='white')

gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135])

# Move the longitude labels off the 135 deg line and onto the 45 degree line
plt.draw()
for ea in gl.label_artists:
    if ea[1] == False:
        tx = ea[2]
        xy = tx.get_position()
        # print(xy)
        if xy[0] == 135:
            tx.set_position([45, xy[1]])

# Find and plot aacgm latitudes
in_lon = np.arange(0, 360, 1)
heights = np.asarray([0] * len(in_lon))
lats = [reference_hemisphere.value * x for x in [80, 70, 60, 50, 40]]
for lat in lats:
    in_lat = np.asarray([lat] * len(in_lon))
    out_lats, out_lons, out_rs = convert_latlon_arr(in_lat=in_lat, in_lon=in_lon,
                                                    height=heights, dtime=date, method_code="A2G")
    # plt.plot(out_lons, out_lats, 'k', markersize=3, transform=ccrs.Geodetic())

# Find the geodetic coordinates for the geomagnetic pole, and plot it as a black dot
pole_geodetic_lat_N, pole_geodetic_lon_E, _ = aacgmv2.convert_latlon(in_lat=reference_hemisphere.value * 90, in_lon=-90,
                                                                     height=0,
                                                                     dtime=date, method_code="A2G")
plt.plot([pole_geodetic_lon_E, pole_geodetic_lon_E], [pole_geodetic_lat_N, pole_geodetic_lat_N],
         'ko', markersize=3, transform=ccrs.Geodetic(), label="Geomagnetic Pole")

ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)

# Plot a blue dot at the geographic pole
plt.plot([-90, -90], [90, 90], 'bo', markersize=3, transform=ccrs.Geodetic(),
         label="Geographic Pole")

beam_range = (0, 15)
gate_range = (0, 74)

for station in stations:
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    hemisphere = this_radars_info.hemisphere

    if hemisphere != reference_hemisphere:
        warnings.warn(station + " is in a different hemisphere, it is being omitted.")
        continue

    radar_id = this_radars_info.hardware_info.stid
    radar_lon = this_radars_info.hardware_info.geographic.lon
    radar_lat = this_radars_info.hardware_info.geographic.lat
    print("Looking at " + station)
    print("     Radar lat: " + str(radar_lat))
    print("     Radar lon: " + str(radar_lon))

    # Plot the radar as a dot
    plt.plot([radar_lon, radar_lon], [radar_lat, radar_lat], "o", markersize=3, transform=ccrs.Geodetic(),
             label=this_radars_info.name)

    cell_corners_lats, cell_corners_lons = radar_fov(stid=radar_id, coords='geo', date=date)

    # plot all the beam boundary lines
    for beam_line in range(beam_range[0], beam_range[1] + 2):
        plt.plot(cell_corners_lons[gate_range[0]:gate_range[1] + 2, beam_line],
                 cell_corners_lats[gate_range[0]:gate_range[1] + 2, beam_line],
                 color='black', linewidth=0.1, transform=ccrs.Geodetic(), zorder=4)

    # plot the arcs boundary lines
    for range_ in range(gate_range[0], gate_range[1] + 2):
        plt.plot(cell_corners_lons[range_, beam_range[0]:beam_range[1] + 2],
                 cell_corners_lats[range_, beam_range[0]:beam_range[1] + 2],
                 color='black', linewidth=0.1, transform=ccrs.Geodetic(), zorder=4)

    # Compute cell centroids
    tf = TimezoneFinder()
    fan_shape = cell_corners_lons.shape
    cell_centers_lons = np.empty(shape=(fan_shape[0] - 1, fan_shape[1] - 1))
    cell_centers_lats = np.empty(shape=(fan_shape[0] - 1, fan_shape[1] - 1))
    cell_time_offsets = np.empty(shape=(fan_shape[0] - 1, fan_shape[1] - 1), dtype='timedelta64[h]')

    for gate_corner in range(fan_shape[0] - 1):
        for beam_corner in range(fan_shape[1] - 1):
            cent_lon, cent_lat = centroid([(cell_corners_lons[gate_corner, beam_corner],
                                            cell_corners_lats[gate_corner, beam_corner]),
                                           (cell_corners_lons[gate_corner + 1, beam_corner],
                                            cell_corners_lats[gate_corner + 1, beam_corner]),
                                           (cell_corners_lons[gate_corner, beam_corner + 1],
                                            cell_corners_lats[gate_corner, beam_corner + 1]),
                                           (cell_corners_lons[gate_corner + 1, beam_corner + 1],
                                            cell_corners_lats[gate_corner + 1, beam_corner + 1])])
            cell_centers_lons[gate_corner, beam_corner] = cent_lon
            cell_centers_lats[gate_corner, beam_corner] = cent_lat

            # Find the time zone of this cell
            try:
                cell_timezone_str = tf.timezone_at(lng=cent_lon, lat=cent_lat)
            except BaseException as e:
                print("Error on lon=" + str(cent_lon) + " and lat=" + str(cent_lat))
                print("This is for gate: " + str(gate_corner) + " and beam " + str(beam_corner))
                raise e

            cell_timezone = pytz.timezone(cell_timezone_str)
            dt = datetime.datetime.now()
            cell_time_offsets[gate_corner, beam_corner] = np.timedelta64(cell_timezone.utcoffset(dt))

    # Find the time zone of the radar itself
    tf = TimezoneFinder()
    radar_timezone_str = tf.timezone_at(lng=radar_lon, lat=radar_lat)

    radar_timezone = pytz.timezone(radar_timezone_str)
    dt = datetime.datetime.now()
    utc_offset = radar_timezone.utcoffset(dt)
    print("\nThe radar is located in the following timezone:")
    print(radar_timezone_str)
    print(utc_offset)

    print("Here are all the unique values in the fan:")
    print(np.unique(cell_time_offsets))

    np.set_printoptions(threshold=np.inf)
    print(np.unique(cell_time_offsets))


# Plot time zone lines
time_zones = cfeature.NaturalEarthFeature(
    category='cultural',
    name='time_zones',
    scale='10m',
    facecolor='none'
)
ax.add_feature(time_zones, edgecolor='red', label="TZ Lines", alpha=0.3, linewidth=0.5, zorder=5)



# Add in 15 degree lines
time_zone_est_lats = np.linspace(reference_hemisphere.value * 40, reference_hemisphere.value * 90, num=100)
time_zone_est_lons = np.linspace(7.5, 352.5, num=24)

print(time_zone_est_lons)
for lon in time_zone_est_lons:
    in_lons = np.asarray([lon] * len(time_zone_est_lats))

    if lon == time_zone_est_lons[0]:
        plt.plot(in_lons, time_zone_est_lats, 'g-', linewidth=0.5, transform=ccrs.Geodetic(), label="15 deg Lines", zorder=4)
    else:
        plt.plot(in_lons, time_zone_est_lats, 'g-', linewidth=0.5, transform=ccrs.Geodetic(), zorder=4)



plt.legend(bbox_to_anchor=(1.05, 0.8))
plt.show()

# # Reading the timezone using ShapelyFeature of Cartopy
# local_time_zones_feature = ShapelyFeature(shpreader.Reader('tz_world.shp').geometries(),
#                                ccrs.PlateCarree(), edgecolor='black')
# ax.add_feature(local_time_zones_feature)

out_dir = loc_root + "/out"
if reference_hemisphere.value == 1:
    out_file = out_dir + "/localTimeZones_" + "NH_" + str(stations[0])
else:
    out_file = out_dir + "/localTimeZones_" + "SH_" + str(stations[0])
fig.savefig(out_file + ".jpg", format='jpg', dpi=300)
plt.close(fig)

# # Convert aacgm lats to geo for plotting
# heights = np.asarray([250] * len(df['lon']))  # TODO: Figure out what to do about heights
# in_lon = df['lon']
# in_lat = df['lat']
# out_lats, _, _ = convert_latlon_arr(in_lat, in_lon, heights, date_time_est, method_code="A2G")


# tz = tzwhere.tzwhere()
# timezone_str = tz.tzNameAt(latitude=radar_lat, longitude=radar_lon)
# print(timezone_str)

# timezone = pytz.timezone(timezone_str)
# dt = datetime.datetime.now()
# utc_offset = timezone.utcoffset(dt)
# print(utc_offset)
