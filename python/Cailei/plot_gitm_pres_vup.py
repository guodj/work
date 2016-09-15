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


def get_g(h):
    m_earth = 5.9722e24
    gc = 6.67384e-11
    r_earth = 6371e3
    return gc * m_earth / ((h+r_earth)**2)

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


def get_p(n, t):
    k = 1.3806488e-23
    return n * k * t


# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


def get_args(argv):

    filelist = []
    blog = 0
    # var = 18
    lon = 120.0
    out = '.\\'

    for arg in argv:

        bfound = 0

        if not bfound:

            # m = re.match(r'-var=(.*)',arg)
            # if m:
            #     var = int(m.group(1))
            #     bfound = 1

            m = re.match(r'-lon=(.*)',arg)
            if m:
                lon = int(m.group(1))
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
            # 'var': var,
            'lon': lon,
            'out': out,
            'log': blog}

    return args

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------

args = get_args(sys.argv)
filelist = args['filelist']

Var = 'gradP÷rg-base&Vup'

lon = args['lon']

outpath = args['out']
outpath.replace('/', '\\')
if outpath[len(outpath)-1] != '\\':
    outpath += '\\'

maxAlt = 500
minAlt = 100
maxLat = 90
minLat = 50

basefile = r'F:\data\Gitm\substorm\bin\3DALL_t000321_100001.bin'
var = list(range(0, 8))
var += [15, 18]


data = read_gitm_one_file(basefile, var)

# [nLons, nLats, nAlts] = data[0].shape
Alts = data[2][0, 0, :] / 1000.0
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

nLats = eLat - sLat
nAlts = eAlt - sAlt

dh = []
hs = []
for ih in range(0, nAlts-2):
    dh.append(Alts[sAlt+ih+2] - Alts[sAlt+ih])
    hs.append((Alts[sAlt + ih + 2] + Alts[sAlt + ih]) / 2)

nn = np.array([])
nn.resize(nLats, nAlts)
for v in range(4, 8):
    nn += np.array(data[v][iLon, sLat:eLat, sAlt:eAlt])

rg = np.array(data[3][iLon, sLat:eLat, sAlt + 1:eAlt - 1])
temp = np.array(data[15][iLon, sLat:eLat, sAlt:eAlt])
gp = np.array([])
gp.resize(nLats, nAlts - 2)

for ih in range(0, nAlts - 2):
    gp[:, ih] = (get_p(nn[:, ih + 2], temp[:, ih + 2]) - get_p(nn[:, ih], temp[:, ih])) / (dh[ih] * 1000)
    rg[:, ih] *= get_g(hs[ih] * 1000)

d2db = -gp / rg

norm_p = cm.colors.Normalize(vmax=0.055, vmin=-0.055)
cmap_p = cm.get_cmap('seismic')
norm_v = cm.colors.Normalize(vmax=52, vmin=-52)
cmap_v = cm.get_cmap('seismic')

fig = plt.figure()
i = 0
for file in filelist:
    data.clear()
    print(file)

    nn = np.array([])
    nn.resize(nLats, nAlts)
    data = read_gitm_one_file(file, var)
    for v in range(4, 8):
        nn += np.array(data[v][iLon, sLat:eLat, sAlt:eAlt])

    rg = np.array(data[3][iLon, sLat:eLat, sAlt + 1:eAlt - 1])
    temp = np.array(data[15][iLon, sLat:eLat, sAlt:eAlt])
    gp = np.array([])
    gp.resize(nLats, nAlts - 2)

    for ih in range(0, nAlts-2):
        gp[:, ih] = (get_p(nn[:, ih + 2], temp[:, ih + 2]) - get_p(nn[:, ih], temp[:, ih])) / (dh[ih] * 1000)
        rg[:, ih] *= get_g(hs[ih] * 1000)

    d2d_p = -gp / rg - d2db

    time = data["time"]
    outfile = outpath + Var + time.strftime('_%y%m%d_%H%M%S.png')

    ax1 = fig.add_subplot(211)

    ax1.set_xlim(minLat, maxLat)
    ax1.set_ylim(minAlt, maxAlt)
    title = time.strftime('%b %d, %Y %H:%M:%S')+'; Lon : '+"%d" % round(Lon)
    ax1.set_title(title)

    cax1 = ax1.pcolormesh(Lats[sLat:eLat], hs, np.transpose(d2d_p), norm=norm_p, cmap=cmap_p)
    ax1.grid(True)

    cbar1 = fig.colorbar(cax1)
    cbar1.set_label('gradP÷rg-base', rotation=90)

    d2d_v = np.array(data[18][iLon, sLat:eLat, sAlt:eAlt])

    ax2 = fig.add_subplot(212)

    ax2.set_xlim(minLat, maxLat)
    ax2.set_ylim(minAlt, maxAlt)
    # title = time.strftime('%b %d, %Y %H:%M:%S')+'; Lon : '+"%d" % round(Lon)
    # ax2.set_title(title)

    cax2 = ax2.pcolormesh(Lats[sLat:eLat], Alts[sAlt:eAlt], np.transpose(d2d_v), norm=norm_v, cmap=cmap_v)
    ax2.grid(True)

    cbar2 = fig.colorbar(cax2)
    cbar2.set_label('V(up)', rotation=90)

    fig.savefig(outfile)
    fig.clear()

    i += 1

print('ok')
