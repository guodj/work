import matplotlib.pyplot as plt
from pylab import cm
from gitm_routines_new import *
import re
import sys
from math import pi

rtod = 180.0/pi

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


def get_var_minmax(var):
    minmax = {}
    minmax[3] = {'min':9e-15, 'max':2e-13}
    minmax[15] = {'min':450, 'max':700}
    minmax[18] = {'min':-52, 'max':52}
    return [minmax[var]['min'],minmax[var]['max']]

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


def get_args(argv):

    filelist = []
    IsLog = 0
    var = 18
    alt = 250
    lon = 120.0
    out = '.\\'
    onePlot = 0

    for arg in argv:

        IsFound = 0

        if not IsFound:

            m = re.match(r'-var=(.*)',arg)
            if m:
                var = int(m.group(1))
                IsFound = 1

            m = re.match(r'-alt=(.*)', arg)
            if m:
                alt = int(m.group(1))
                IsFound = 1

            m = re.match(r'-lon=(.*)',arg)
            if m:
                lon = int(m.group(1))
                IsFound = 1

            m = re.match(r'-out=(.*)', arg)
            if m:
                out = m.group(1)
                IsFound = 1

            m = re.match(r'-alog',arg)
            if m:
                IsLog = 1
                IsFound = 1

            m = re.match(r'-oneplot',arg)
            if m:
                onePlot = 1
                IsFound = 1

            if IsFound == 0 and not(arg == argv[0]):
                m = re.search(r'\*',arg)
                if m:
                    filelist += glob(arg)
                else:
                    filelist.append(arg)

    args = {'filelist': filelist,
            'var': var,
            'alt': alt,
            'lon': lon,
            'out': out,
            'IsLog': IsLog,
            'IsOnePlot': onePlot}

    return args

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------

args = get_args(sys.argv)
filelist = args['filelist']
header = read_gitm_file_header(filelist[0])

var = args['var']
Var = header['vars'][var].decode().replace(' ', '')
if var == 18:
    Var = 'Vn(Up)'

vars = [0, 1, 2]
vars.append(var)

alt = args['alt']
lon = args['lon']

oneplot = args['IsOnePlot']

outpath = args['out']
outpath.replace('/', '\\')
if outpath[len(outpath)-1] != '\\':
    outpath += '\\'

maxLat = 90
minLat = 50

if oneplot:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(minLat, maxLat)
    ax.grid(True)

i = 0
for file in filelist:
    if i % 3 != 0:
        i += 1
        continue

    print(file)
    data = read_gitm_one_file(file, vars)

    if i == 0:
        [nLons, nLats, nAlts] = data[0].shape
        Alts = data[2][0][0] / 1000.0
        Lons = data[0][:, 0, 0] * rtod
        Lats = data[1][0, :, 0] * rtod

        iLat = 0
        while Lats[iLat] <= maxLat:
            if Lats[iLat] < minLat:
                sLat = iLat
            iLat += 1
        eLat = iLat + 1

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

    time = data["time"]
    d2d = np.array(data[var][iLon, sLat:eLat, iAlt])

    if not oneplot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(minLat, maxLat)
        ax.grid(True)

    # ax.set_ylim(minAlt, maxAlt)

    # ax.set_xticks(range(minLat, maxLat, 10))
    # ax.set_xticklabels(map(str, range(0, 360, 30)))
    # ax.set_yticks(range(0,40,10))
    # ax.set_yticklabels(map(str,range(90,50,-10)))

    cax = ax.plot(Lats[sLat:eLat], np.transpose(d2d))

    # cbar = fig.colorbar(cax)
    # cbar.set_label(Var, rotation=90)

    # fig.clear()

    if not oneplot:
        title = time.strftime('%b %d, %Y %H:%M:%S') + '; Lon : ' + "%d" % round(Lon) + '; Alt: %dkm' % Alt
        ax.set_title(title)
        outfile = outpath + Var + time.strftime('_%y%m%d_%H%M%S.png')
        fig.savefig(outfile)

    data.clear()

    i += 1

if oneplot:
    title = time.strftime('%b %d, %Y %H:%M:%S') + '; Lon : ' + "%d" % round(Lon) + '; Alt: %dkm' % Alt
    ax.set_title(title)
    outfile = outpath + Var + time.strftime('_%y%m%d_%H%M.png')
    fig.savefig(outfile)

print('ok')
