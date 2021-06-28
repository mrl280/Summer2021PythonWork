from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(ncols=1, nrows=2)

projection = ccrs.NorthPolarStereo()

ax = fig.add_subplot(gs[0, 0], projection=projection)

ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)

plt.show()

