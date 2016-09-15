import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as tf
import matplotlib.dates as dates
from pylab import cm
from gitm_routines_new import *
from math import pi


def rmap(r):
    r = 90 - r
    return r

rtod = 180.0/pi

datapath = 'F:\\data\\gitm\\substorm\\bin\\'
file = datapath + '3DALL_t000321_120001.bin'
header = read_gitm_file_header(file)

Var = 'Rho'
var = 3
vars = [0,1,2]
vars.append(var)

alt = 400
maxLon = 360/rtod
minLon = 0
maxLat = 90
minLat = 50

data = read_gitm_one_file(file, vars)

[nLons, nLats, nAlts] = data[0].shape
Alts = data[2][0][0] / 1000.0
Lons = data[0][:, 0, 0]
Lats = data[1][0, :, 0] * rtod

iLon = 0
while Lons[iLon] <= maxLon:
    if Lons[iLon] < minLon:
        sLon = iLon
    iLon += 1
sLon = sLon
eLon = iLon + 1

iLat = 0
while Lats[iLat] <= maxLat:
    if Lats[iLat] < minLat:
        sLat = iLat
    iLat += 1
sLat = sLat
eLat = iLat + 1

if alt < 50:
    iAlt = alt
else:
    if alt > Alts[nAlts - 1]:
        iAlt = nAlts - 3
    else:
        iAlt = 2
        while Alts[iAlt] < alt:
            iAlt += 1
Alt = Alts[iAlt]

time = data["time"]
d2d = np.array(data[var][sLon:eLon, sLat:eLat, iAlt])

maxi = 2.e-13  # np.max(AllData2D)*1.05
mini = 9.e-15  # np.min(AllData2D)*0.95
#dr = (maxi-mini)/31
#levels = np.arange(mini, maxi, dr)

i = 0
fig = plt.figure()
# norm = cm.colors.Normalize(vmax=maxi, vmin=mini)
norm = cm.colors.LogNorm(vmax=maxi, vmin=mini)
cmap = cm.get_cmap('Blues')

ax = fig.add_subplot(111, projection='polar')


# outfile = datapath + Var + time.strftime('_%y%m%d_%H%M%S.png')
outfile = 'test' + time.strftime('_%y%m%d_%H%M%S.png')

ax.set_ylim(0,40)
ax.set_xlim([0,360])
title = time.strftime('%b %d, %Y %H:%M:%S')+'; Alt : '+"%.2f" % Alt + ' km'
ax.set_title(title)

ax.set_xticks(np.radians(range(0,360,30)))
ax.set_xticklabels(map(str,range(0,360,30)))
ax.set_yticks(range(0,40,10))
ax.set_yticklabels(map(str,range(90,50,-10)))

cax = ax.pcolormesh(Lons[sLon:eLon], rmap(Lats[sLat:eLat]), np.transpose(d2d), norm=norm, cmap=cmap)
ax.grid(True)

cbar = fig.colorbar(cax)
cbar.set_label(Var,rotation=90)

fig.savefig(outfile)
fig.clear()

i=i+1

print('ok')
