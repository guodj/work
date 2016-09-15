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
    minmax[3] = {'min': 9e-15, 'max': 2e-13}
    minmax[15] = {'min': 450, 'max': 700}
    minmax[18] = {'min': -52, 'max': 52}
    minmax[33] = {'min': 0., 'max': 5e11}
    minmax[36] = {'min': -800., 'max': 800}
    return [minmax[var]['min'], minmax[var]['max']]

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


def get_args(argv):

    filelist = []
    IsLog = 0
    var = 18
    lon = 120.0
    out = '.\\'

    for arg in argv:

        IsFound = 0

        if not IsFound:

            m = re.match(r'-var=(.*)',arg)
            if m:
                var = int(m.group(1))
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

            if IsFound == 0 and not(arg == argv[0]):
                m = re.search(r'\*',arg)
                if m:
                    filelist += glob(arg)
                else:
                    filelist.append(arg)

    args = {'filelist': filelist,
            'var': var,
            'lon': lon,
            'out': out,
            'IsLog': IsLog}

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

lon = args['lon']

outpath = args['out']
outpath.replace('/', '\\')
if outpath[len(outpath)-1] != '\\':
    outpath += '\\'

maxAlt = 500
minAlt = 100
maxLat = 90
minLat = 50

[vmin, vmax] = get_var_minmax(var)

if args['IsLog']:
    norm = cm.colors.LogNorm(vmax=vmax, vmin=vmin)
else:
    norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)

if var == 18 or var == 36:
    cmap = cm.get_cmap('seismic')
elif var == 15:
    cmap = cm.get_cmap('YlOrRd')
elif var == 3 or var == 33:
    cmap = cm.get_cmap('Blues')
else:
    cmap = cm.get_cmap('jet')

fig = plt.figure()

i = 0
for file in filelist:
    print(file)
    data = read_gitm_one_file(file, vars)

    if i == 0:
        [nLons, nLats, nAlts] = data[0].shape
        Alts = data[2][0][0] / 1000.0
        Lons = data[0][:, 0, 0] * rtod
        Lats = data[1][0, :, 0] * rtod

        iAlt = 0
        while Alts[iAlt] <= maxAlt:
            if Alts[iAlt] < minAlt:
                sAlt = iAlt
            iAlt += 1
        eAlt = iAlt + 1

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

    time = data["time"]
    d2d = np.array(data[var][iLon, sLat:eLat, sAlt:eAlt])

    ax = fig.add_subplot(111)

    outfile = outpath + Var + time.strftime('_%y%m%d_%H%M%S.png')

    ax.set_xlim(minLat, maxLat)
    ax.set_ylim(minAlt, maxAlt)
    title = time.strftime('%b %d, %Y %H:%M:%S')+'; Lon : '+"%d" % round(Lon)
    ax.set_title(title)

    # ax.set_xticks(range(minLat, maxLat, 10))
    # ax.set_xticklabels(map(str, range(0, 360, 30)))
    # ax.set_yticks(range(0,40,10))
    # ax.set_yticklabels(map(str,range(90,50,-10)))

    levels = [-800, -700, -600, -500, -400, -300, -200, -100, 0, 100, 200, 300, 400, 500, 600, 700, 800]
    cax = ax.contour(Lats[sLat:eLat], Alts[sAlt:eAlt], np.transpose(d2d), norm=norm, cmap=cmap, levels=levels)
    ax.grid(True)

    cbar = fig.colorbar(cax)
    cbar.set_label(Var, rotation=90)

    fig.savefig(outfile)
    fig.clear()

    data.clear()

    i += 1

print('ok')
