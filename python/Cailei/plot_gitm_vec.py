import matplotlib.pyplot as plt
from pylab import cm
from gitm_routines_new import *
import re
import sys
import math

rtod = 180.0/math.pi
r_earth = 5.e2  # 6371e3

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


def rmap_rad(r):
    r = math.pi/2 - r
    return r


def vn_to_lat(v):
    v /= r_earth


def ve_to_lon(v, lat):
    v /= (r_earth * np.cos(lat))

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


def get_args(argv):

    filelist = []
    blog = 0
    var = 15
    alt = 400
    out = '.\\'

    for arg in argv:

        bfound = 0

        if not bfound:

            m = re.match(r'-var=(.*)', arg)
            if m:
                var = m.group(1)
                bfound = 1

            m = re.match(r'-alt=(.*)', arg)
            if m:
                alt = int(m.group(1))
                bfound = 1

            m = re.match(r'-out=(.*)', arg)
            if m:
                out = m.group(1)
                bfound = 1

            m = re.match(r'-alog', arg)
            if m:
                blog = 1
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
            'out': out,
            'log': blog}

    return args

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


args = get_args(sys.argv)

filelist = args['filelist']
var = args['var']
alt = args['alt']
outpath = args['out']
outpath.replace('/','\\')
if outpath[len(outpath)-1] != '\\':
    outpath += '\\'

maxLon = 360/rtod
minLon = 0
maxLat = 90/rtod
minLat = 50/rtod

vars_to_read = [0, 1, 2]

data = read_gitm_one_file(filelist[0], vars_to_read)

[nlons, nlats, nalts] = data[0].shape
alts_all = data[2][0, 0, :] / 1000.0
lons_all = data[0][:, 0, 0]
lats_all = data[1][0, :, 0]

ilon = 0
while lons_all[ilon] <= maxLon:
    if lons_all[ilon] < minLon:
        ilon_s = ilon
    ilon += 1
ilon_e = ilon + 1

ilat = 0
while lats_all[ilat] <= maxLat:
    if lats_all[ilat] < minLat:
        ilat_s = ilat
    ilat += 1
ilat_e = ilat + 1

if alt < 50:
    ialt = alt
else:
    if alt > alts_all[nalts - 1]:
        ialt = nalts - 3
    else:
        ialt = 2
        while alts_all[ialt] < alt:
            ialt += 1
Alt = alts_all[ialt]

# lons = lons_all[ilon_s:ilon_e]
# lats = lats_all[ilat_s:ilat_e]
lats2d, lons2d = np.meshgrid(lats_all[ilat_s:ilat_e-8:8], lons_all[ilon_s:ilon_e:4])

if var == 'vn':
    vars_to_read = [16, 17]
elif var == 'vi':
    vars_to_read = [36, 37]

# cmap = cm.get_cmap('seismic')
# norm = cm.colors.Normalize(vmax=50, vmin=-50)

fig = plt.figure()

i = 0
for file in filelist:
    print(file)
    data = read_gitm_one_file(file, vars_to_read)
    # data = read_gitm_one_file(file, [18])

    time = data["time"]
    # d2d = np.array(data[18][ilon_s:ilon_e, ilat_s:ilat_e, ialt])
    ve2d = np.array(data[vars_to_read[0]][ilon_s:ilon_e:4, ilat_s:ilat_e-8:8, ialt])
    vn2d = np.array(data[vars_to_read[1]][ilon_s:ilon_e:4, ilat_s:ilat_e-8:8, ialt])

    ve_to_lon(ve2d, lats2d)
    vn_to_lat(vn2d)

    ax = fig.add_subplot(111, projection='polar')

    outfile = outpath + time.strftime('Vi_%y%m%d_%H%M%S.png')

    ax.set_ylim(0, 40/rtod)
    # ax.set_xlim(0, 360)
    title = time.strftime('%b %d, %Y %H:%M:%S')+'; Alt : '+"%.2f" % Alt + ' km'
    ax.set_title(title)

    ax.set_xticks(np.radians(range(0, 360, 30)))
    ax.set_xticklabels(map(str, range(0, 360, 30)))
    ax.set_yticks(np.radians(range(0, 40, 10)))
    ax.set_yticklabels(map(str, range(90, 50, -10)))

    # q = ax.quiver(np.array([math.pi/3, math.pi/6]), np.array([math.pi/6, math.pi/9]), np.array([1, .2]), np.array([1, .8]))
    q = ax.quiver(lons2d, rmap_rad(lats2d), ve2d, vn2d)
    # q = ax.quiver(lons2d, rmap_rad(lats2d), lons2d+0.5, rmap_rad(lats2d)+0.5)
    # qk = ax.quiverkey(q, 0.9, 0.95, 100, 'test', coordinates='figure')

    # cax = ax.pcolormesh(lons2d, rmap_rad(lats2d), d2d, norm=norm, cmap=cmap)
    ax.grid(True)

    # cbar = fig.colorbar(cax)

    fig.savefig(outfile)
    fig.clear()

    data.clear()

    i += 1
