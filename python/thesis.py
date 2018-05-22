from aacgmv2 import convert
from aacgmv2 import convert_mlt
import datetime as dt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import gitm_create_coordinate as gcc

plt.figure(figsize=(6.88,6.74))
#geographic coordinates
ax1, projection1 = gcc.create_map(
        2, 2, 1, 'pol', 90, 50, 0, coastlines=False,useLT=True,
        dlat=10, lonticklabel=(1, 1, 1, 1))
ax1.scatter(0, 90, color='r', transform=ccrs.PlateCarree())
ax2, projection2 = gcc.create_map(
        2, 2, 2, 'pol', -50, -90, 0, coastlines=False,useLT=True,
        dlat=10, lonticklabel=(1, 1, 1, 1))
ax2.scatter(0, -90, color='b', transform=ccrs.PlateCarree())

#geomagnetic coordinates
ax3, projection3 = gcc.create_map(
        2, 2, 3, 'pol', 90, 50, 0, coastlines=False,useLT=True,
        dlat=10, lonticklabel=(1, 1, 1, 1))
mlatn,mlonn = convert(90,0,0,date=dt.date(2002,3,21))
for k in range(24):
    mltn = convert_mlt(mlonn[0],dtime=dt.datetime(2003,3,21,k))
    ax3.scatter(mltn*15,mlatn[0],color='r',transform=ccrs.PlateCarree())
ax3.scatter(180,75,s=50,c='k',marker='x',transform=ccrs.PlateCarree())

ax4, projection4 = gcc.create_map(
        2, 2, 4, 'pol', -50, -90, 0, coastlines=False,useLT=True,
        dlat=10, lonticklabel=(1, 1, 1, 1))
mlats,mlons = convert(-90,0,0,date=dt.date(2002,3,21))
for k in range(24):
    mlts = convert_mlt(mlons[0],dtime=dt.datetime(2003,3,21,k))
    ax4.scatter(mlts*15,mlats[0],color='b',transform=ccrs.PlateCarree())
ax4.scatter(180,-75,s=50,c='k',marker='x',transform=ccrs.PlateCarree())
plt.savefig('/Users/guod/Documents/Thesis/006SciFig/work3_001.eps')
#plt.show()
