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


def rmap(r):
    r = 90 - r
    return r

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
    var = 15
    alt = 400.0
    out = '.\\'

    for arg in argv:

        IsFound = 0

        if not IsFound:

            m = re.match(r'-var=(.*)',arg)
            if m:
                var = int(m.group(1))
                IsFound = 1

            m = re.match(r'-alt=(.*)',arg)
            if m:
                alt = int(m.group(1))
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
            'alt': alt,
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
Var = header['vars'][var].decode().replace(' ','')
if var == 18:
    Var = 'Vn(Up)'

vars = [0, 1, 2]
vars.append(var)

alt = args['alt']

outpath = args['out']
outpath.replace('/','\\')
if outpath[len(outpath)-1] != '\\':
    outpath += '\\'

maxLon = 360/rtod
minLon = 0
maxLat = 90
minLat = 50

[vmin, vmax] = get_var_minmax(var)

if args['IsLog']:
    norm = cm.colors.LogNorm(vmax=vmax, vmin=vmin)
else:
    norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)

if var == 18:
    cmap = cm.get_cmap('seismic')
elif var == 15:
    cmap = cm.get_cmap('YlOrRd')
elif var == 3:
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
        Lons = data[0][:, 0, 0]
        Lats = data[1][0, :, 0] * rtod

        iLon = 0
        while Lons[iLon] <= maxLon:
            if Lons[iLon] < minLon:
                sLon = iLon
            iLon += 1
        eLon = iLon + 1

        iLat = 0
        while Lats[iLat] <= maxLat:
            if Lats[iLat] < minLat:
                sLat = iLat
            iLat += 1
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

    ax = fig.add_subplot(111, projection='polar')

    outfile = outpath + Var + time.strftime('_%y%m%d_%H%M%S.png')

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

    data.clear()

    i += 1

print('ok')
