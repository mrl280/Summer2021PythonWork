import numpy as np

from aacgmv2 import convert_latlon_arr
from matplotlib import pyplot as plt
import datetime as datetime

date = datetime.datetime.now()

# The aacgm zero degree longitude line varies with latitude

# For each lat, we need to find the required rotation


lats_aagcm = np.linspace(40, 90, num=100)
lons_aagcm = np.asarray([0] * len(lats_aagcm))  # The zero degree line in aacgm
heights = np.asarray([0] * len(lats_aagcm))
out_lats, out_lons, out_rs = convert_latlon_arr(lats_aagcm, lons_aagcm, heights, date, method_code="A2G")

heights = np.asarray([250] * len(lats_aagcm))

out_lats_w_height, out_lons_w_height, out_rs_w_height = convert_latlon_arr(lats_aagcm, lons_aagcm, heights, date,
                                                                           method_code="A2G")

fig2, ax2 = plt.subplots(nrows=1, ncols=1, dpi=300)
ax2.set_xlabel("Lat")
ax2.set_ylabel("Lon")

ax2.plot(out_lats, out_lons, 'r-', label='0 km')
ax2.plot(out_lats_w_height, out_lons_w_height, 'b-', label="250 km")


ax2.legend(loc='best')
plt.show()
