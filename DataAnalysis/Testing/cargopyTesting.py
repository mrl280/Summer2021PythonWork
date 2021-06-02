import pathlib
import cartopy
import aacgmv2

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
import datetime as dt

from datetime import datetime
from cartopy.feature.nightshade import Nightshade
from pydarn import radar_fov

rkn_lon, rkn_lat = -92.113, 62.82  # RKN
inv_lon, inv_lat = -133.77, 68.41  # INV
risr_lon, risr_lat = -94.91, 74.73  # RISR

fig = plt.figure(figsize=(8, 5), dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(rkn_lon, rkn_lat))
ax.set_global()
ax.set_extent([-37, -143, 35, 90])
ax.stock_img()
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle="-", linewidth=0.5)
radar_marker_size = 4

# Add RKN to the map
plt.plot([rkn_lon, rkn_lon], [rkn_lat, rkn_lat], 'go', markersize=radar_marker_size, transform=ccrs.Geodetic())

# Add INV to the map
plt.plot([inv_lon, inv_lon], [inv_lat, inv_lat], 'bo', markersize=radar_marker_size, transform=ccrs.Geodetic())

# Add RISR to the map
plt.plot([risr_lon, risr_lon], [risr_lat, risr_lat], 'ro', markersize=radar_marker_size, transform=ccrs.Geodetic(),
         zorder=3)

date_time = datetime.utcnow()
ax.add_feature(Nightshade(date_time, alpha=0.25))
ax.set_title("Night time shading for " + str(date_time) + " UT")

# Try to plot RKN's FoV  # RKN is id 65
ranges = [0, 100]

rkn_beam_corners_lats, rkn_beam_corners_lons = radar_fov(65, coords='geo')
rkn_fan_shape = rkn_beam_corners_lons.shape

inv_beam_corners_lats, inv_beam_corners_lons = radar_fov(64, coords='geo')
inv_fan_shape = inv_beam_corners_lons.shape

# plot all the RKN beam boundary lines
for beam_line in range(rkn_fan_shape[1]):
    plt.plot(rkn_beam_corners_lons[0:ranges[1] + 1, beam_line], rkn_beam_corners_lats[0:ranges[1] + 1, beam_line],
             color='black', linewidth=0.1, transform=ccrs.Geodetic())

# # left boundary line
# plt.plot(rkn_beam_corners_lons[0:ranges[1] + 1, 0], rkn_beam_corners_lats[0:ranges[1] + 1, 0],
#          color='black', linewidth=0.5, transform=ccrs.Geodetic())
#
# # right boundary line
# plt.plot(rkn_beam_corners_lons[0:ranges[1] + 1, rkn_fan_shape[1] - 1], rkn_beam_corners_lats[0:ranges[1] + 1, rkn_fan_shape[1] - 1],
#          color='black', linewidth=0.5, transform=ccrs.Geodetic())

# plot all the RKN arcs  # have to loop through ranges
for range_ in range(ranges[1] + 1):
    plt.plot(rkn_beam_corners_lons[range_, 0:rkn_fan_shape[1]], rkn_beam_corners_lats[range_, 0:rkn_fan_shape[1]],
             color='black', linewidth=0.1, transform=ccrs.Geodetic())

# # bottom arc
# plt.plot(rkn_beam_corners_lons[0, 0:rkn_fan_shape[1] - 1], rkn_beam_corners_lats[0, 0:rkn_fan_shape[1] - 1],
#          color='black', linewidth=0.5, transform=ccrs.Geodetic())
#
# # top radar arc
# # Note: pydarn is plotting the top arc one arc too soon
# plt.plot(rkn_beam_corners_lons[ranges[1], 0:rkn_fan_shape[1]], rkn_beam_corners_lats[ranges[1], 0:rkn_fan_shape[1]],
#          color='black', linewidth=0.5, transform=ccrs.Geodetic())


# plot all the INV beam boundary lines
for beam_line in range(inv_fan_shape[1]):
    plt.plot(inv_beam_corners_lons[0:ranges[1] + 1, beam_line], inv_beam_corners_lats[0:ranges[1] + 1, beam_line],
             color='black', linewidth=0.1, transform=ccrs.Geodetic())

# plot all the INV arcs  # have to loop through ranges
for range_ in range(ranges[1] + 1):
    plt.plot(inv_beam_corners_lons[range_, 0:inv_fan_shape[1]], inv_beam_corners_lats[range_, 0:inv_fan_shape[1]],
             color='black', linewidth=0.1, transform=ccrs.Geodetic())

# Read in a RISR world day file and plot all of the points - hoping this will show where beams are
station = "ran"
year = 2011
month = 11
day = 11
resolution = 1  # min

loc_root = str(pathlib.Path().absolute().parent.absolute())
in_dir = loc_root + "/DataReading/RISR/data/" + station + "/" + station + str(year) + str(month) + str(day)
in_file = in_dir + "/" + station + str(year) + str(month) + str(day) + "." + str(resolution) + "min.pkl"
df = pd.read_pickle(in_file)

print(df.keys())

# At RKN the overlap is beam 5 and gates 25-32
# At INV the overlap is beam 12 and gates 30-38

plt.scatter(df.query('wdBmnum == 2')['gdlon'], df.query('wdBmnum == 2')['gdlat'],
            s=1, c='r', marker='.', transform=ccrs.Geodetic(), zorder=3)

plt.scatter(df.query('wdBmnum == 10')['gdlon'], df.query('wdBmnum == 10')['gdlat'],
            s=1, c='r', marker='.', transform=ccrs.Geodetic(), zorder=3)

# Highlight RKN beam 5
for beam_line in [5, 6]:
    plt.plot(rkn_beam_corners_lons[0:ranges[1] + 1, beam_line], rkn_beam_corners_lats[0:ranges[1] + 1, beam_line],
             color='black', linewidth=0.5, transform=ccrs.Geodetic())

# Highlight INV beam 12
for beam_line in [12, 13]:
    plt.plot(inv_beam_corners_lons[0:ranges[1] + 1, beam_line], inv_beam_corners_lats[0:ranges[1] + 1, beam_line],
             color='black', linewidth=0.5, transform=ccrs.Geodetic())

plt.show()
