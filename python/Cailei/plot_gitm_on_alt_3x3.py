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


# def get_var_minmax(var):
#     minmax = {}
#     minmax[3] = {'min': 9e-15, 'max': 2e-13}
#     minmax[15] = {'min': 450, 'max': 700}
#     minmax[18] = {'min': -52, 'max': 52}
#     return [minmax[var]['min'], minmax[var]['max']]


def get_var_minmax(var, alt):
    minmax = {}
    minmax[3] = {250: {'min': 9e-15, 'max': 2e-13}, 400: {'min': 9e-15, 'max': 2e-13}}
    minmax[15] = {250: {'min': 450, 'max': 700}, 400: {'min': 450, 'max': 700}}
    minmax[18] = {250: {'min': -20, 'max': 20}, 400: {'min': -50, 'max': 50}}
    minmax[33] = {250: {'min': 0., 'max': 5e11}, 400: {'min': 0., 'max': 5e11}}
    minmax[36] = {250: {'min': -800., 'max': 800}, 400: {'min': -800., 'max': 800}}
    return [minmax[var][alt]['min'], minmax[var][alt]['max']]

# -----------------------------[------------------------------------------------
#
# -----------------------------------------------------------------------------


def get_args(argv):

    filelist = []
    blog = 0
    var = 15
    alt = -1
    out = '.\\'

    for arg in argv:

        bfound = 0

        if not bfound:

            m = re.match(r'-var=(.*)',arg)
            if m:
                var = int(m.group(1))
                bfound = 1

            m = re.match(r'-alt=(.*)',arg)
            if m:
                alt = int(m.group(1))
                bfound = 1

            m = re.match(r'-out=(.*)', arg)
            if m:
                out = m.group(1)
                bfound = 1

            m = re.match(r'-alog',arg)
            if m:
                blog = 1
                bfound = 1

            if bfound == 0 and not(arg == argv[0]):
                m = re.search(r'\*',arg)
                if m:
                    filelist += glob(arg)
                else:
                    filelist.append(arg)

    args = {'filelist': filelist,
            'var': var,
            'alt': alt,
            'out': out,
            'log': blog}

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

alt = args['alt']

outpath = args['out']
outpath.replace('/','\\')
if outpath[len(outpath)-1] != '\\':
    outpath += '\\'

maxLon = 360/rtod
minLon = 0
maxLat = 90
minLat = 50

if var == 18 or var == 36:
    cmap = cm.get_cmap('seismic')
elif var == 15:
    cmap = cm.get_cmap('YlOrRd')
elif var == 3 or var == 33:
    cmap = cm.get_cmap('Blues')
else:
    cmap = cm.get_cmap('jet')

vars = [0, 1, 2]

data = read_gitm_one_file(filelist[0], vars)

[nLons, nLats, nAlts] = data[0].shape
Alts = data[2][0, 0, :] / 1000.0
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

if alt < 0:
    altrange = [250, 400]
else:
    altrange = [alt]

fig = plt.figure()
fig.set_size_inches(18, 18)

for alt in altrange:
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

    [vmin, vmax] = get_var_minmax(var, alt)

    if args['log']:
        norm = cm.colors.LogNorm(vmax=vmax, vmin=vmin)
    else:
        norm = cm.colors.Normalize(vmax=vmax, vmin=vmin)

    # time = data['time']
    # title = time.strftime('%b %d, %Y') + '; Alt : ' + "%.2f" % Alt + ' km'
    # plt.title(title)
    # fig.text(0, 0, title)
    # plt.title(title)

    i = 0
    for file in filelist:
        print(file)
        data = read_gitm_one_file(file, [var])

        time = data["time"]
        d2d = np.array(data[var][sLon:eLon, sLat:eLat, iAlt])

        ax = fig.add_subplot(331+i, projection='polar')

        ax.set_ylim(0, 40)
        ax.set_xlim([0, 360])
        # title = time.strftime('%b %d, %Y %H:%M')+'; Alt : '+"%.2f" % Alt + ' km'
        title = time.strftime('%H:%M')
        ax.set_title(title, loc='right')

        ax.set_xticks(np.radians(range(0,360,30)))
        ax.set_xticklabels(map(str,range(0,360,30)))
        ax.set_yticks(range(0,40,10))
        ax.set_yticklabels(map(str,range(90,50,-10)))

        cax = ax.pcolormesh(Lons[sLon:eLon], rmap(Lats[sLat:eLat]), np.transpose(d2d), norm=norm, cmap=cmap)
        ax.grid(True)

        # cbar.set_label(Var, rotation=90)

        data.clear()

        i += 1

    # cbar = fig.colorbar(cax)

    outfile = outpath + Var + time.strftime('_%y%m%d_')+'Alt_'+"%.0f" % Alt + '.png'
    fig.savefig(outfile)
    fig.clear()

print('ok')
