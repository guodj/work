import matplotlib.pyplot as plt
from gitm_routines_new import *
import re
import sys
from math import pi

rtod = 180.0/pi

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


def get_args(argv):

    filelist = []
    var = 18
    alt = 400
    lon = 120.0
    lat = 75.0
    out = '.\\'

    for arg in argv:

        bfound = 0

        if not bfound:

            m = re.match(r'-var=(.*)', arg)
            if m:
                var = int(m.group(1))
                bfound = 1

            m = re.match(r'-alt=(.*)', arg)
            if m:
                alt = int(m.group(1))
                bfound = 1

            m = re.match(r'-lon=(.*)', arg)
            if m:
                lon = int(m.group(1))
                bfound = 1

            m = re.match(r'-lat=(.*)', arg)
            if m:
                lat = int(m.group(1))
                bfound = 1

            m = re.match(r'-out=(.*)', arg)
            if m:
                out = m.group(1)
                bfound = 1

            if bfound == 0 and not(arg == argv[0]):
                m = re.search(r'\*', arg)
                if m:
                    filelist += glob(arg)
                else:
                    filelist.append(arg)

    args = {'filelist': filelist,
            'var': var,
            'alt': alt,
            'lon': lon,
            'lat': lat,
            'out': out}

    return args

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------

args = get_args(sys.argv)
filelist = args['filelist']
var = args['var']
alt = args['alt']
lon = args['lon']
lat = args['lat']

outpath = args['out']
outpath.replace('/', '\\')
if outpath[len(outpath)-1] != '\\':
    outpath += '\\'

vars = [0, 1, 2]
data = read_gitm_one_file(filelist[0], vars)

[nLons, nLats, nAlts] = data[0].shape
Alts = data[2][0, 0, :] / 1000.0
Lons = data[0][:, 0, 0] * rtod
Lats = data[1][0, :, 0] * rtod

iLat = 0
while Lats[iLat] < lat:
    iLat += 1
Lat = Lats[iLat]

iLon = 0
while Lons[iLon] < lon:
    iLon += 1
Lon = Lons[iLon]

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

dln = {}
il = 0
lons = Lons[::2]
for l in lons:
    dln[round(l)] = []

i = 0
time = []
for file in filelist:
    print(file)
    data = read_gitm_one_file(file, [var])
    time.append(data['time'])
    il = 0
    for l in lons:
        dln[round(l)].append(data[var][il, iLat, iAlt])
        il += 2

for l in lons:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid(True)
    ax.plot(np.array(dln[round(l)]))
    title = '; Lat : ' + "%d" % Lat + '; Lon : ' + "%d" % round(l) + '; Alt: %dkm' % Alt
    ax.set_title(title)
    outfile = outpath + 'Vn(up)' + 'Lon_' + '%d' % round(l)
    fig.savefig(outfile)
    fig.clear()

