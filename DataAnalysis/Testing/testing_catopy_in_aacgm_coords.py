import math
import pathlib

import aacgmv2
import pydarn
import calendar

import cartopy.crs as ccrs
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import numpy as np
import cartopy.feature as cfeature
import datetime as datetime

from aacgmv2 import convert_latlon_arr
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pydarn import radar_fov, SuperDARNRadars

from DataAnalysis.EchoOccurrence.lib.data_getters.get_local_dummy_data import get_local_dummy_data
from DataAnalysis.EchoOccurrence.lib.only_keep_45km_res_data import only_keep_45km_res_data

station = "dce"

# Get radar specific hardware information
all_radars_info = SuperDARNRadars()
this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
hemisphere = this_radars_info.hemisphere
radar_id = this_radars_info.hardware_info.stid
radar_lon = this_radars_info.hardware_info.geographic.lon
radar_lat = this_radars_info.hardware_info.geographic.lat
print("Radar lat: " + str(radar_lat))
print("Radar lon: " + str(radar_lon))

# # Get some dummy data
# df = get_local_dummy_data(station="dce", year=2012, month=11, day=12, start_hour_UT=0, end_hour_UT=24, occ_data=False)
# df = df.loc[(df['p_l'] >= 3)]  # Restrict to points with at least 3 dB
# df.reset_index(drop=True, inplace=True)
# df = only_keep_45km_res_data(df)

date = datetime.datetime.now()
# date = datetime.datetime(1968, 5, 17)



print("Preparing the plot...")
fig = plt.figure(figsize=(5, 5), dpi=300)
lat_extreme = 40 * hemisphere.value  # deg
if hemisphere.value == 1:
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
elif hemisphere.value == -1:
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())
else:
    raise Exception("hemisphere not recognized")

ax.set_extent([-180, 180, hemisphere.value * 90, lat_extreme], crs=ccrs.PlateCarree())

# Find and plot aacgm latitudes
in_lon = np.arange(0, 360, 1)
heights = np.asarray([0] * len(in_lon))
lats = [hemisphere.value * x for x in [80, 70, 60, 50]]
for lat in lats:
    in_lat = np.asarray([lat] * len(in_lon))
    out_lats, out_lons, out_rs = convert_latlon_arr(in_lat=in_lat, in_lon=in_lon,
                                                    height=heights, dtime=date, method_code="A2G")
    plt.plot(out_lons, out_lats, 'k', markersize=3, transform=ccrs.Geodetic(), label=str(lat) + " deg")

# Find the geodetic coordinates for the geomagnetic pole, and plot it as a black dot
pole_geodetic_lat_N, pole_geodetic_lon_E, _ = aacgmv2.convert_latlon(in_lat=hemisphere.value * 90, in_lon=-90, height=0,
                                                                     dtime=date, method_code="A2G")
plt.plot([pole_geodetic_lon_E, pole_geodetic_lon_E], [pole_geodetic_lat_N, pole_geodetic_lat_N],
         'ko', markersize=3, transform=ccrs.Geodetic(), label="Geomagnetic North Pole")


ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)

# Plot a blue dot at the geographic north pole
plt.plot([-90, -90], [90, 90], 'bo', markersize=3, transform=ccrs.Geodetic(),
         label="North Pole")


# Plot the radar as a red dot
plt.plot([radar_lon, radar_lon], [radar_lat, radar_lat], 'ro', markersize=3, transform=ccrs.Geodetic(),
         label=this_radars_info.name)

# Compute a circle in axis coordinates which can be used as a boundary
theta = np.linspace(0, 2 * np.pi, 100)
center, radius = [0.5, 0.5], 0.5
vertices = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(vertices * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)

# Add gridlines and mlt labels
text_offset_multiplier = 1.03
gl = ax.gridlines(draw_labels=True, linestyle='--', color='white')
gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135])
ax.text(0, text_offset_multiplier * ax.get_ylim()[1], "12", ha='center', va='bottom')
ax.text(0, text_offset_multiplier * ax.get_ylim()[0], "00", ha='center', va='top')
ax.text(text_offset_multiplier * ax.get_xlim()[1], 0, "06", ha='left', va='center')
ax.text(text_offset_multiplier * ax.get_xlim()[0], 0, "18", ha='right', va='center')

beam_range = (0, 15)
gate_range = (0, 74)

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

# Plot the midnight line
ax.plot([0, 0], [hemisphere.value * 40, hemisphere.value * 90], 'g-', transform=ccrs.Geodetic())

# Plot the 6 o-clock line
ax.plot([6 * 15, 6 * 15], [hemisphere.value * 40, hemisphere.value * 90], 'm-', transform=ccrs.Geodetic())

# Plot the zero degree of aacgm lon line
lats_aagcm = np.linspace(hemisphere.value * 40, hemisphere.value * 90, num=100)
lons_aagcm = np.asarray([0] * len(lats_aagcm))  # The zero degree line in aacgm
heights = np.asarray([250] * len(lats_aagcm))
out_lats, out_lons, out_rs = convert_latlon_arr(lats_aagcm, lons_aagcm, heights, date, method_code="A2G")

ax.plot(out_lons, out_lats, 'k-', transform=ccrs.Geodetic())

# Plot the 90 degree of aacgm lon line
lons_aagcm = np.asarray([90] * len(lats_aagcm))
out_lats, out_lons, out_rs = convert_latlon_arr(lats_aagcm, lons_aagcm, heights, date, method_code="A2G")

ax.plot(out_lons, out_lats, 'k-', linewidth=0.5, transform=ccrs.Geodetic())

plt.show()


loc_root = str((pathlib.Path().parent.absolute()))
out_dir = loc_root + "/out"
out_file = out_dir + "/" + station + "_geometry_w_aacgm"
# fig.savefig(out_file + ".jpg", format='jpg', dpi=300)
plt.close(fig)


# # Convert aacgm lats to geo for plotting
# heights = np.asarray([250] * len(df['lon']))  # TODO: Figure out what to do about heights
# in_lon = df['lon']
# in_lat = df['lat']
# out_lats, _, _ = convert_latlon_arr(in_lat, in_lon, heights, date_time_est, method_code="A2G")




