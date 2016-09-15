from gitm_routines_new import *

rtod = 180.0/3.141592

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------

def get_gitm_vars_range(file_to_read, iAlt, vars_to_read=-1):
    data = read_gitm_one_file(file_to_read, vars_to_read)
    minmax = {}
    if (vars_to_read == -1):
        vars_to_read = np.arange(39)
    for iVar in vars_to_read:
        minmax[iVar] = [np.min(data[iVar][:,:,iAlt]), np.max(data[iVar][:,:,iAlt])]

    return minmax

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------

filepath = 'F:\\data\\gitm\\'
filelist = glob(filepath+'3DALL*.bin')

minmax = {}
vars = [0,1,2,3,15,18]
i = 0
alt = 400;
maxLon = 360/rtod
minLon = 0
maxLat = 90
minLat = 50
for file in filelist:
    print(file)
    data = read_gitm_one_file(file, vars)
    if (i == 0):
        [nLons, nLats, nAlts] = data[0].shape
        Alts = data[2][0][0] / 1000.0
        Lons = data[0][:, 0, 0]
        Lats = data[1][0, :, 0] * rtod

        iLon = 0
        while (Lons[iLon] <= maxLon):
            if (Lons[iLon] < minLon):
                sLon = iLon
            iLon = iLon + 1
        eLon = iLon + 1

        iLat = 0
        while (Lats[iLat] <= maxLat):
            if (Lats[iLat] < minLat):
                sLat = iLat
            iLat = iLat + 1
        eLat = iLat + 1

        if (alt < 50):
            iAlt = alt
        else:
            if (alt > Alts[nAlts - 1]):
                iAlt = nAlts - 3
            else:
                iAlt = 2
                while (Alts[iAlt] < alt):
                    iAlt = iAlt + 1
    m = {}
    for iVar in vars:
        d2d = data[iVar][sLon:eLon, sLat:eLat, iAlt]
        m[iVar] = [np.min(d2d), np.max(d2d)]
    #print(m)
    if (i == 0):
        minmax = m;
    else:
        for iVar in vars:
            minmax[iVar] = [min(minmax[iVar][0], m[iVar][0]), max(minmax[iVar][1], m[iVar][1])]
    i += 1

print(minmax)

print('ok')
