# Auroral precipation model using SSUSI data
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import myfunctions as mf
import scipy.signal as signal
import seaborn as sns
sns.set()

# global variables
LLP1 = 200  # Lower limit of pixels in one bin
LLEFTOTH = 200  # lower limit of total energy flux in one hemisphere
LLN = 10  # Lower limit of boundaries found in one hemisphere
BinMLT = 0.25  # MLT range in one bin
LLEFTOT = 10  # lower limit of total energy flux in each bin, for method 1,2
# used by the median filter, can not be neither too large or small
kernel_size = 3  # for method 5, use its own kernel_size
smoothps = 5
datapath = '/home/guod/WD4T/ssusi/'
savepath = '/home/guod/WD4T/work4/'


def image_energy_flux_one_file(
        ax, fn, ns, vmin=None, vmax=None, s=2, alpha=0.8, kernel_size=1):
    '''
    plot ssusi data in a polar coordinate, only for energy flux
    input:
        fn      : file name
        ns      : 'N' or 'S', for NH or SH
        vmin    : min value of var in colorbar
        vmax    : max value of var in colorbar
        s       : scatter points size in points^2
        alpha   : 0 (transparent), 1 (opaque)
    output:
        ax      : axis handle
        hs      : handle of the scatter plot
    Note: for both north and south poles, the parameter 'ns' in set_polar is
          'N'
    '''
    ssusi = Dataset(fn)
    MLat = np.array(ssusi['LATITUDE_GEOMAGNETIC_GRID_MAP'])
    MLT = np.array(ssusi['MLT_GRID_MAP'])
    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE  # value in a no data bin
    if ns == 'N':
        var = 'ENERGY_FLUX_NORTH_MAP'
        ut = ssusi['UT_N'][:]
    if ns == 'S':
        var = 'ENERGY_FLUX_SOUTH_MAP'
        ut = ssusi['UT_S'][:]
    variable = np.array(ssusi[var])
    variable = signal.medfilt2d(variable, kernel_size)
    # Only include data points inside the boundary
    fp = (ut != nodatavalue)
    MLat, MLT, variable = (k[fp] for k in (MLat, MLT, variable))
    # plot variable
    r = 90-MLat
    theta = MLT/12*np.pi
    hs = ax.scatter(theta, r, s=s, c=variable, vmin=vmin,
                    vmax=vmax, alpha=alpha, cmap='Oranges')
    # Set polar coordinates
    mf.set_polar(
            ax, ns='N', boundinglat=50,
            dlat=10, dlon=90, useLT=True)
    ax.set_title(ssusi[var].TITLE)
    ax.text(0.85, 0.05,
            'YEAR: '+str(np.int(ssusi['YEAR'][:])) +
            '\nDOY: '+str(np.int(ssusi['DOY'][:])) +
            '\nORBIT: '+str(np.int(ssusi.STARTING_ORBIT_NUMBER)),
            transform=ax.transAxes,
            fontsize=8)
    # cb = plt.colorbar(hs, ax=ax, pad=0.1)
    # cb.set_label(ssusi[var].UNITS)
    return ax, hs


def find_parameters_one_file_1(fn):
    '''
    Get the peak value, peak latitude, and full width of the aurora.
    The aurora oval should contain 80% of the total energy flux.
    input:
        fn        : ssusi aurora file name
    output:
        pp        : the parameters
    '''
    # 1, Read data
    ssusi = Dataset(fn)
    # --------------------------------------------------
    # --------------adjustable parameters---------------
    PEF = 0.8  # proportion of total energy flux a auroral oval should contain.
    # --------------------------------------------------
    # --------------------------------------------------
    # 2, Read variables that will be used
    # orbitn = ssusi.STARTING_ORBIT_NUMBER # orbit number
    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE  # value in a no data bin
    # magnetic latitude
    MLat = np.array(ssusi['LATITUDE_GEOMAGNETIC_GRID_MAP'])
    # magnetic local time
    MLT = np.array(ssusi['MLT_GRID_MAP'])
    # electron energy flux
    energyfluxn = np.array(ssusi['ENERGY_FLUX_NORTH_MAP'])
    energyfluxs = np.array(ssusi['ENERGY_FLUX_SOUTH_MAP'])
    utn = np.array(ssusi['UT_N'])
    uts = np.array(ssusi['UT_S'])
    # boundary values
    # bv = np.array(ssusi['ELECTRON_ENERGY_FLUX_THRESHOLDS'])
    # Apply median filter to the energy flux
    energyfluxn = signal.medfilt2d(energyfluxn, kernel_size)
    energyfluxs = signal.medfilt2d(energyfluxs, kernel_size)

    # 3, exclude pixels without data
    fpn = utn != nodatavalue
    MLatn, MLTn, energyfluxn, utn = (
            k[fpn] for k in (MLat, MLT, energyfluxn, utn))
    fps = uts != nodatavalue
    MLats, MLTs, energyfluxs, uts = (
            k[fps] for k in (MLat, MLT, energyfluxs, uts))

    # 4, Change UT if part of the previous day is included.
    starttime = ssusi.STARTING_TIME   # format: yyyydddhhmmss
    stoptime = ssusi.STOPPING_TIME
    if starttime[:7] != stoptime[:7]:
        utn[utn > 20] = utn[utn > 20]-24
        uts[uts > 20] = uts[uts > 20]-24
    # Now, utn can be negative, negative valus mean yesterday

    # 5, initialize output
    pp = pd.DataFrame(
            np.ones([192, 5])*np.nan,
            columns=['pbmlat', 'mlatm', 'ebmlat', 'peak', 'datetime'],
            index=pd.MultiIndex.from_product(
                    [['N', 'S'], np.arange(BinMLT/2, 24, BinMLT)]))

    # 6, find parameters
    # 6.1, group data according to hemispheres
    for k00, k0 in enumerate(['N', 'S']):
        # select data in the specified hemisphere
        if k0 == 'N':
            mlat0, mlt0, ef0, ut0 = (MLatn, MLTn, energyfluxn, utn)
            # bvv = bv[0]
        if k0 == 'S':
            mlat0, mlt0, ef0, ut0 = (MLats, MLTs, energyfluxs, uts)
            # bvv = bv[1]

        # 6.2, calculate the cumulative energy flux in the hemisphere,
        # then find the lower and upper limits of the mlat in one bin
        idx = np.argsort(mlat0)
        mlat00, ef00 = mlat0[idx], ef0[idx]
        ef00 = ef00.cumsum()
        if ef00[-1] < LLEFTOTH:
            continue
        sf = 0.9
        fpt = ef00 < (0.5-sf/2)*ef00[-1]
        llmlat = mlat00[fpt][-1]
        fpt = ef00 > (0.5+sf/2)*ef00[-1]
        ulmlat = mlat00[fpt][0]

        # 6.3, group data according to MLT
        for k11, k1 in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
            # Select data in the specified MLT range
            lmlt, rmlt = k1-BinMLT/2, k1+BinMLT/2
            fp = (mlt0 >= lmlt) & (mlt0 <= rmlt)
            if np.sum(fp) <= LLP1:  # lower limit of pixels in one bin
                continue
            mlat1, ef1, ut1 = mlat0[fp], ef0[fp], ut0[fp]
            if mlat1.min() > llmlat or mlat1.max() < ulmlat:
                continue

            # 6.4, sort MLat
            idx = np.argsort(mlat1)
            mlat1, ef1, ut1 = mlat1[idx], ef1[idx], ut1[idx]

            # 6.5, calculate the integral energy flux
            efcum = ef1.cumsum()
            eftot = efcum[-1]
            if eftot < LLEFTOT:
                continue

            # 6.6, find equatorward and poleward boundary
            fpt = efcum < (0.5-PEF/2)*eftot
            ebmlat = mlat1[fpt][-1] if fpt.any() else mlat1[0]
            fpt = efcum > (0.5+PEF/2)*eftot
            pbmlat = mlat1[fpt][0] if fpt.any() else mlat1[-1]

            # 6.7, find max energy flux between ebmlat and pbmlat
            fpt = (mlat1 > ebmlat) & (mlat1 < pbmlat)
            idm = np.argmax(ef1[fpt])
            efm = ef1[fpt][idm]
            mlatm = mlat1[fpt][idm]

            # 6.8, find ut
            utt = np.median(ut1[fpt])

            pp.loc[(k0, k1), :] = [pbmlat, mlatm, ebmlat, efm, utt]

    # 7, remove defective points
    for k00, k0 in enumerate(('N', 'S')):
        # If pbmlat or ebmlat is far apart from that of the neighboring
        # boundary location, it is excluded.
        points = 25  # should be big enough
        mlatt = 5
        pp0 = pp.loc[(k0, slice(None)), ['pbmlat', 'ebmlat']].values.copy()
        ppt = []  # initialize
        for k1 in range(points):
            ppt.append(np.roll(pp0, k1-int((points-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp1 = np.nanmedian(pp1, axis=2)
        fpt = (np.abs(pp1-pp0) > mlatt).any(axis=1)
        pp2 = pp.loc[(k0, slice(None))].values
        pp2[fpt] = np.nan
        pp.loc[(k0, slice(None))] = pp2

    # 8, smooth result
    for k00, k0 in enumerate(('N', 'S')):
        pp0 = pp.loc[(k0, slice(None))].values.copy()
        ppt = []  # initialize
        for k1 in range(smoothps):
            ppt.append(np.roll(pp0, k1-int((smoothps-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp0 = np.median(pp1, axis=2)
        pp.loc[(k0, slice(None))] = pp0

    # 9, If the total number of MLTs < some vlaue, we think there is no aurora
    for k00, k0 in enumerate(('N', 'S')):
        pp0 = pp.loc[(k0, slice(None)), 'pbmlat'].values
        if np.sum(~np.isnan(pp0)) <= LLN:
            pp.loc[(k0, slice(None))] = np.nan

    # 10, calculate datetime
    yyyy = int(ssusi['YEAR'][:])
    ddd = int(ssusi['DOY'][:])
    stime = pd.Timestamp(yyyy, 1, 1, 0, 0, 0)
    dt1 = pd.Timedelta(ddd-1, 'D')
    dt2 = pd.TimedeltaIndex(pp['datetime'].values, 'h')
    pp['datetime'] = stime+dt1+dt2
    return pp


def find_parameters_one_file_2(fn):
    '''
    Get the peak value, peak latitude, and full width of the aurora.
    The aurora oval should contain 80% of the total energy flux.
    input:
        fn        : ssusi aurora file name
    output:
        pp        : the parameters
    '''
    # 1, Read data
    ssusi = Dataset(fn)
    # --------------------------------------------------
    # --------------adjustable parameters---------------
    PEFE = 0.1  # proportion of the total energy flux at boundary
    # --------------------------------------------------
    # --------------------------------------------------
    # 2, Read variables that will be used
    # orbitn = ssusi.STARTING_ORBIT_NUMBER  # orbit number
    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE  # value in a no data bin
    # magnetic latitude
    MLat = np.array(ssusi['LATITUDE_GEOMAGNETIC_GRID_MAP'])
    MLT = np.array(ssusi['MLT_GRID_MAP'])  # magnetic local time
    energyfluxn = np.array(ssusi['ENERGY_FLUX_NORTH_MAP'])
    energyfluxs = np.array(ssusi['ENERGY_FLUX_SOUTH_MAP'])
    utn = np.array(ssusi['UT_N'])
    uts = np.array(ssusi['UT_S'])
    # bv = np.array(ssusi['ELECTRON_ENERGY_FLUX_THRESHOLDS']) # boundary values
    # Apply median filter to the energy flux
    energyfluxn = signal.medfilt2d(energyfluxn, kernel_size)
    energyfluxs = signal.medfilt2d(energyfluxs, kernel_size)

    # 3, exclude pixels without data
    fpn = utn != nodatavalue
    MLatn, MLTn, energyfluxn, utn = (
            k[fpn] for k in (MLat, MLT, energyfluxn, utn))
    fps = uts != nodatavalue
    MLats, MLTs, energyfluxs, uts = (
            k[fps] for k in (MLat, MLT, energyfluxs, uts))

    # 4, Change UT if part of the previous day is included.
    starttime = ssusi.STARTING_TIME   # format: yyyydddhhmmss
    stoptime = ssusi.STOPPING_TIME
    if starttime[:7] != stoptime[:7]:
        utn[utn > 20] = utn[utn > 20]-24
        uts[uts > 20] = uts[uts > 20]-24
    # Now, utn can be negative, negative valus mean yesterday

    # 5, initialize output
    pp = pd.DataFrame(
            np.ones([192, 5])*np.nan,
            columns=['pbmlat', 'mlatm', 'ebmlat', 'peak', 'datetime'],
            index=pd.MultiIndex.from_product(
                    [['N', 'S'], np.arange(BinMLT/2, 24, BinMLT)]))

    # 6, find parameters
    # 6.1, group data according to hemispheres
    for k00, k0 in enumerate(['N', 'S']):
        # select data in the specified hemisphere
        if k0 == 'N':
            mlat0, mlt0, ef0, ut0 = (MLatn, MLTn, energyfluxn, utn)
            # bvv = bv[0]
        if k0 == 'S':
            mlat0, mlt0, ef0, ut0 = (MLats, MLTs, energyfluxs, uts)
            # bvv = bv[1]

        # 6.2, calculate the cumulative energy flux in the hemisphere,
        # then find the lower and upper limits of the mlat in one bin
        idx = np.argsort(mlat0)
        mlat00, ef00 = mlat0[idx], ef0[idx]
        ef00 = ef00.cumsum()
        if ef00[-1] < LLEFTOTH:
            continue
        sf = 0.9
        fpt = ef00 < (0.5-sf/2)*ef00[-1]
        llmlat = mlat00[fpt][-1]
        fpt = ef00 > (0.5+sf/2)*ef00[-1]
        ulmlat = mlat00[fpt][0]

        # 6.3, group data according to MLT
        for k11, k1 in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
            # Select data in the specified MLT range
            lmlt, rmlt = k1-BinMLT/2, k1+BinMLT/2
            fp = (mlt0 >= lmlt) & (mlt0 <= rmlt)
            if np.sum(fp) <= LLP1:
                continue
            mlat1, ef1, ut1 = mlat0[fp], ef0[fp], ut0[fp]
            if mlat1.min() > llmlat or mlat1.max() < ulmlat:
                continue

            # 6.4, sort MLat
            idx = np.argsort(mlat1)
            mlat1, ef1, ut1 = mlat1[idx], ef1[idx], ut1[idx]

            # 6.5, calculate the integral energy flux
            efcum = ef1.cumsum()
            eftot = efcum[-1]
            if eftot < LLEFTOT:
                continue

            # 6.6, find equatorward and poleward boundary
            fpt = efcum < 0.5*eftot
            mlatmed = mlat1[fpt][-1]
            fpt = efcum > PEFE*eftot
            ebmlat = mlat1[fpt][0]
            pbmlat1 = 2*mlatmed-ebmlat
            fpt = efcum < (1-PEFE)*eftot
            pbmlat2 = mlat1[fpt][-1]
            if pbmlat2 < pbmlat1:
                pbmlat = pbmlat2
            else:
                pbmlat = pbmlat1
                eft1 = efcum[mlat1 < pbmlat1][-1]
                eft2 = efcum[mlat1 < pbmlat2][-1]
                sft = 0.5
                eftm = efcum[mlat1 < pbmlat1+sft*(pbmlat2-pbmlat1)][-1]
                if eftm - eft1 > sft*(eft2 - eft1):
                    pbmlat = pbmlat1+sft*(pbmlat2-pbmlat1)

            # 6.7, find max energy flux between ebmlat and pbmlat
            fpt = (mlat1 > ebmlat) & (mlat1 < pbmlat)
            if np.sum(fpt) == 0:
                continue
            idm = np.argmax(ef1[fpt])
            efm = ef1[fpt][idm]
            mlatm = mlat1[fpt][idm]

            # 6.8, find ut
            utt = np.median(ut1[fpt])

            pp.loc[(k0, k1), :] = [pbmlat, mlatm, ebmlat, efm, utt]

    # 7, remove defective points
    for k00, k0 in enumerate(('N', 'S')):
        # If pbmlat or ebmlat is far apart from that of the neighboring
        # boundary location, it is excluded.
        points = 25  # should be big enough
        mlatt = 5
        pp0 = pp.loc[(k0, slice(None)), ['pbmlat', 'ebmlat']].values.copy()
        ppt = []  # initialize
        for k1 in range(points):
            ppt.append(np.roll(pp0, k1-int((points-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp1 = np.nanmedian(pp1, axis=2)
        fpt = (np.abs(pp1-pp0) > mlatt).any(axis=1)
        pp2 = pp.loc[(k0, slice(None))].values
        pp2[fpt] = np.nan
        pp.loc[(k0, slice(None))] = pp2

    # 8, smooth result
    for k00, k0 in enumerate(('N', 'S')):
        pp0 = pp.loc[(k0, slice(None))].values.copy()
        ppt = []  # initialize
        for k1 in range(smoothps):
            ppt.append(np.roll(pp0, k1-int((smoothps-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp0 = np.median(pp1, axis=2)
        pp.loc[(k0, slice(None))] = pp0

    # 9, If the total number of MLTs < some vlaue, we think there is no aurora
    for k00, k0 in enumerate(('N', 'S')):
        pp0 = pp.loc[(k0, slice(None)), 'pbmlat'].values
        if np.sum(~np.isnan(pp0)) <= LLN:
            pp.loc[(k0, slice(None))] = np.nan

    # 10, calculate datetime
    yyyy = int(ssusi['YEAR'][:])
    ddd = int(ssusi['DOY'][:])
    stime = pd.Timestamp(yyyy, 1, 1, 0, 0, 0)
    dt1 = pd.Timedelta(ddd-1, 'D')
    dt2 = pd.TimedeltaIndex(pp['datetime'].values, 'h')
    pp['datetime'] = stime+dt1+dt2
    return pp


def find_parameters_one_file_3(fn):
    '''
    Get the peak value, peak latitude, and full width of the aurora
    The edge of the auroral oval has energy flux <= 20(or 30)% of the maximum
    energy flux. The max energy flux is the same for all the MLT
    input:
        fn        : ssusi aurora file name
    output:
        pp        : the parameters
    '''
    # 1, Read data
    ssusi = Dataset(fn)

    # --------------------------------------------------
    # --------------adjustable parameters---------------
    LLP2 = 20  # minimum pixels in auroral oval
    LLMEF = 1  # Lower limit of the maximum energy flux
    EdgeF = 0.2  # percentage of max energy flux at boundaries
    # Lower limit of pixel proportion with ef > edgeef for the 2nd aurora
    LLA2F = 0.8
    BinW = 0.2  # size of the mlat bin window
    # 'EdgeM' of 'EdgeN' consecutive data points are smaller than
    # boundary values
    EdgeN = 5
    EdgeM = 5
    LLdMLat = 2  # Lower limit of mlat difference of two auroral ovals
    # --------------------------------------------------
    # --------------------------------------------------

    # 2, Read variables that will be used
    # orbitn = ssusi.STARTING_ORBIT_NUMBER  # orbit number
    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE  # value in a no data bin
    # magnetic latitude
    MLat = np.array(ssusi['LATITUDE_GEOMAGNETIC_GRID_MAP'])
    MLT = np.array(ssusi['MLT_GRID_MAP'])  # magnetic local time
    # electron energy flux
    energyfluxn = np.array(ssusi['ENERGY_FLUX_NORTH_MAP'])
    energyfluxs = np.array(ssusi['ENERGY_FLUX_SOUTH_MAP'])
    utn = np.array(ssusi['UT_N'])
    uts = np.array(ssusi['UT_S'])
    # bv = ssusi['ELECTRON_ENERGY_FLUX_THRESHOLDS'][:] # boundary values
    # Apply median filter to the energy flux
    energyfluxn = signal.medfilt2d(energyfluxn, kernel_size)
    energyfluxs = signal.medfilt2d(energyfluxs, kernel_size)

    # 3, exclude pixels without data
    fpn = utn != nodatavalue
    MLatn, MLTn, energyfluxn, utn = (
            k[fpn] for k in (MLat, MLT, energyfluxn, utn))
    fps = uts != nodatavalue
    MLats, MLTs, energyfluxs, uts = (
            k[fps] for k in (MLat, MLT, energyfluxs, uts))

    # 4, Change UT if part of the previous day is included.
    starttime = ssusi.STARTING_TIME   # format: yyyydddhhmmss
    stoptime = ssusi.STOPPING_TIME
    if starttime[:7] != stoptime[:7]:
        utn[utn > 20] = utn[utn > 20]-24
        uts[uts > 20] = uts[uts > 20]-24
    # Now, utn can be negative, negative valus mean yesterday

    # 5, initialize output
    pp = pd.DataFrame(np.ones([192, 9])*np.nan,
                      columns=['pbmlat1', 'mlatm1', 'ebmlat1',
                               'pbmlat2', 'mlatm2', 'ebmlat2',
                               'peak1', 'peak2', 'datetime'],
                      index=pd.MultiIndex.from_product(
                              [['N', 'S'], np.arange(BinMLT/2, 24, BinMLT)]))

    # 6, find parameters
    # 6.1, group data according to hemispheres
    for k00, k0 in enumerate(['N', 'S']):
        # select data in the specified hemisphere
        if k0 == 'N':
            mlat0, mlt0, ef0, ut0 = (MLatn, MLTn, energyfluxn, utn)
            # bvv = bv[0]
        if k0 == 'S':
            mlat0, mlt0, ef0, ut0 = (MLats, MLTs, energyfluxs, uts)
            # bvv = bv[1]

        # 6.2, calculate the cumulative energy flux in the hemisphere,
        # then find the lower and upper limits of the mlat in one bin
        idx = np.argsort(mlat0)
        mlat00, ef00 = mlat0[idx], ef0[idx]
        ef00 = ef00.cumsum()
        if ef00[-1] < LLEFTOTH:
            continue
        sf = 0.9
        fpt = ef00 < (0.5-sf/2)*ef00[-1]
        llmlat = mlat00[fpt][-1]
        fpt = ef00 > (0.5+sf/2)*ef00[-1]
        ulmlat = mlat00[fpt][0]

        # 6.3, Find boundary values for the whole hemisphere.
        maxh = []
        for k11, k1 in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
            # Select data of required MLT range
            lmlt, rmlt = k1-BinMLT/2, k1+BinMLT/2
            fp = (mlt0 >= lmlt) & (mlt0 <= rmlt)
            if np.sum(fp) <= LLP1:
                continue
            maxv = np.max(ef0[fp])
            if maxv > LLMEF:
                maxh.append(maxv)
        if len(maxh) < LLN:
            continue
        else:
            # the edge energy flux of the hemisphere
            edgev = EdgeF*np.median(maxh)

        # 6.4, group data according to MLT
        for k11, k1 in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
            peak1, peak2 = False, False
            # Select data in the specified MLT range
            lmlt, rmlt = k1-BinMLT/2, k1+BinMLT/2
            fp = (mlt0 >= lmlt) & (mlt0 <= rmlt)
            if np.sum(fp) <= LLP1:
                continue
            mlat1, ef1, ut1 = mlat0[fp], ef0[fp], ut0[fp]

            # jump to next loop when only part of the aurora is captured
            if mlat1.min() > llmlat or mlat1.max() < ulmlat:
                continue

            # 6.5, sort MLat
            idx = np.argsort(mlat1)
            mlat1, ef1, ut1 = mlat1[idx], ef1[idx], ut1[idx]

            # 6.6, bin energy flux
            # Is this necessary?
            mlat11 = np.arange(50+BinW/2, 90, BinW)
            ef11 = np.arange(50+BinW/2, 90, BinW)*np.nan
            for k22, k2 in enumerate(mlat11):
                lmlat, rmlat = k2-BinW/2, k2+BinW/2
                fpt = (mlat1 > lmlat) & (mlat1 < rmlat)
                if np.sum(fpt) > 0:
                    ef11[k22] = ef1[fpt].max()  # envelope of energy flux
            fpt = ~np.isnan(ef11)
            mlat11, ef11 = mlat11[fpt], ef11[fpt]

            # 6.7, Find maximum energy flux
            fpmax1 = np.argmax(ef11)
            mlatm1, efm1 = mlat11[fpmax1], ef11[fpmax1]
            if efm1 < LLMEF:
                continue

            # 6.8, Find boundaries
            # 1, energy flux < some value for some consecutive values
            # or 2, difference of Mlat > some value
            # poleward boundary
            efef = np.ones([len(ef11), EdgeN])*np.nan
            for k2 in range(EdgeN):
                efef[:, k2] = np.roll(ef11, -k2)
                efef[-(k2+1):, k2] = 0
            fpt = np.sum(efef < edgev, axis=1) >= EdgeM
            idx2 = np.where(fpt)[0]
            idxpb = idx2[idx2 >= fpmax1].min()
            # equatorward boundary
            efef = np.ones([len(ef11), EdgeN])*np.nan
            for k2 in range(EdgeN):
                efef[:, k2] = np.roll(ef11, k2)
                efef[:k2+1, k2] = 0
            fpt = np.sum(efef < edgev, axis=1) >= EdgeM
            idx2 = np.where(fpt)[0]
            idxeb = idx2[idx2 <= fpmax1].max()
            # equatorward and poleward boundary mlat
            ebmlat1, pbmlat1 = mlat11[idxeb], mlat11[idxpb]

            # 6.9, check the existence of the second auroral
            # If the numbers of (the residual pixels with energy flux > edgev)
            # > some value then, the second maximum exists
            mlat22, ef22 = mlat11.copy(), ef11.copy()
            fpt = (mlat11 >= ebmlat1) & (mlat11 <= pbmlat1)
            mlat22[fpt] = 0
            ef22[fpt] = 0
            if np.sum(ef22 > edgev) > LLA2F*np.sum(fpt):
                # 6.10, Find the second maximum
                fpmax2 = np.argmax(ef22)
                mlatm2, efm2 = mlat22[fpmax2], ef22[fpmax2]

                # 6.11, Find the second boundary
                # poleward boundary
                efef = np.ones([len(ef22), EdgeN])*np.nan
                for k2 in range(EdgeN):
                    efef[:, k2] = np.roll(ef22, -k2)
                    efef[-(k2+1):, k2] = 0
                fpt = np.sum(efef < edgev, axis=1) >= EdgeM
                idx2 = np.where(fpt)[0]
                idxpb = idx2[idx2 >= fpmax2].min()
                # equatorward boundary
                efef = np.ones([len(ef22), EdgeN])*np.nan
                for k2 in range(EdgeN):
                    efef[:, k2] = np.roll(ef22, k2)
                    efef[:k2+1, k2] = 0
                fpt = np.sum(efef < edgev, axis=1) >= EdgeM
                idx2 = np.where(fpt)[0]
                idxeb = idx2[idx2 <= fpmax2].max()
                ebmlat2, pbmlat2 = mlat22[idxeb], mlat22[idxpb]

                # if the two auroral ovals are close, they are combined
                if (ebmlat2-pbmlat1 <= LLdMLat) and \
                   (ebmlat1-pbmlat2 <= LLdMLat):
                    ebmlat1 = np.min([ebmlat1, ebmlat2])
                    pbmlat1 = np.max([pbmlat1, pbmlat2])
                else:
                    fpt = (mlat1 > ebmlat2) & (mlat1 < pbmlat2)
                    if np.sum(fpt) > LLP2:
                        peak2 = True  # The second auroral peak exists
            fpt = (mlat1 > ebmlat1) & (mlat1 < pbmlat1)
            if np.sum(fpt) >= LLP2:
                peak1 = True  # The second auroral peak exists

            if peak1 and peak2:
                ebmlat2, pbmlat2, ebmlat1, pbmlat1 = np.sort(
                        [ebmlat2, pbmlat2, ebmlat1, pbmlat1])
                idx = np.argsort([mlatm1, mlatm2])
                mlatm2, mlatm1 = np.array([mlatm1, mlatm2])[idx]
                efm2, efm1 = np.array([efm1, efm2])[idx]
            elif peak1:
                ebmlat2, pbmlat2, ebmlat1, pbmlat1 = (
                        ebmlat1, np.nan, np.nan, pbmlat1)
                mlatm2, mlatm1 = np.nan, mlatm1
                efm2, efm1 = np.nan, efm1
            elif peak2:
                ebmlat2, pbmlat2, ebmlat1, pbmlat1 = (
                        ebmlat2, np.nan, np.nan, pbmlat2)
                mlatm2, mlatm1 = np.nan, mlatm2
                efm2, efm1 = np.nan, efm2
            else:
                (ebmlat2, pbmlat2, ebmlat1, pbmlat1, mlatm2, mlatm1,
                 efm2, efm1) = np.ones(8)*np.nan
            if peak1 or peak2:
                utt = np.median(ut1[fpt])
            else:
                utt = np.nan
            pp.loc[(k0, k1), :] = [
                    pbmlat1, mlatm1, ebmlat1, pbmlat2, mlatm2, ebmlat2,
                    efm1, efm2, utt]
    # 7, remove defective points
    for k00, k0 in enumerate(('N', 'S')):
        # If pbmlat1 or ebmlat2 is `mlatt` degrees greater or less than the
        # median value of neighboring `points` points, the corresponding
        # points are removed
        points = 25
        mlatt = 5
        pp0 = pp.loc[(k0, slice(None)), ['pbmlat1', 'ebmlat2']].values
        ppt = []  # initialize
        for k1 in range(points):
            ppt.append(np.roll(pp0, k1-int((points-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp1 = np.nanmedian(pp1, axis=2)
        fpt = (np.abs(pp1-pp0) > mlatt).any(axis=1)
        pp2 = pp.loc[(k0, slice(None))].values
        pp2[fpt] = np.nan
        pp.loc[(k0, slice(None))] = pp2

    # 8, smooth result
    for k00, k0 in enumerate(('N', 'S')):
        pp0 = pp.loc[(k0, slice(None))].values
        ppt = []  # initialize
        for k1 in range(smoothps):
            ppt.append(np.roll(pp0, k1-int((smoothps-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp0 = np.median(pp1, axis=2)
        pp.loc[(k0, slice(None))] = pp0

    # 9, If the total number of MLTs < some vlaue, we think there is no auroral
    for k00, k0 in enumerate(('N', 'S')):
        pp0 = pp.loc[(k0, slice(None)), 'pbmlat1'].values
        if np.sum(~np.isnan(pp0)) <= LLN:
            pp.loc[(k0, slice(None))] = np.nan

    # 10, calculate datetime
    yyyy = int(ssusi['YEAR'][:])
    ddd = int(ssusi['DOY'][:])
    stime = pd.Timestamp(yyyy, 1, 1, 0, 0, 0)
    dt1 = pd.Timedelta(ddd-1, 'D')
    dt2 = pd.TimedeltaIndex(pp['datetime'].values, 'h')
    pp['datetime'] = stime+dt1+dt2
    return pp


def find_parameters_one_file_4(fn):
    '''
    Get the peak value, peak latitude, and full width of the aurora
    The edge of the auroral oval has energy flux < 0.5
    input:
        fn        : ssusi aurora file name
    output:
        pp        : the parameters
    '''
    # 1, Read data
    ssusi = Dataset(fn)

    # --------------------------------------------------
    # --------------adjustable parameters---------------
    LLP2 = 20  # minimum pixels in auroral oval
    LLA2F = 0.8
    LLMEF = 1  # Lower limit of the maximum energy flux
    BinW = 0.2  # size of the mlat bin window
    EdgeEF = 0.5  # energy flux at boundary
    # 'EdgeM' of 'EdgeN' consecutive data points are smaller than boundary
    # values
    EdgeN = 5
    EdgeM = 5
    LLdMLat = 2  # Lower limit of mlat difference of two auroral ovals, degree
    # --------------------------------------------------
    # --------------------------------------------------

    # 2, Read variables that will be used
    # orbitn = ssusi.STARTING_ORBIT_NUMBER # orbit number
    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE  # value in a no data bin
    # magnetic latitude
    MLat = np.array(ssusi['LATITUDE_GEOMAGNETIC_GRID_MAP'])
    MLT = np.array(ssusi['MLT_GRID_MAP'])  # magnetic local time
    # electron energy flux
    energyfluxn = np.array(ssusi['ENERGY_FLUX_NORTH_MAP'])
    energyfluxs = np.array(ssusi['ENERGY_FLUX_SOUTH_MAP'])
    utn = np.array(ssusi['UT_N'])
    uts = np.array(ssusi['UT_S'])
    # bv = ssusi['ELECTRON_ENERGY_FLUX_THRESHOLDS'][:]  # boundary values
    # Apply median filter to the energy flux
    energyfluxn = signal.medfilt2d(energyfluxn, kernel_size)
    energyfluxs = signal.medfilt2d(energyfluxs, kernel_size)

    # 3, exclude pixels without data
    fpn = utn != nodatavalue
    MLatn, MLTn, energyfluxn, utn = (
            k[fpn] for k in (MLat, MLT, energyfluxn, utn))
    fps = uts != nodatavalue
    MLats, MLTs, energyfluxs, uts = (
            k[fps] for k in (MLat, MLT, energyfluxs, uts))

    # 4, Change UT if part of the previous day is included.
    starttime = ssusi.STARTING_TIME   # format: yyyydddhhmmss
    stoptime = ssusi.STOPPING_TIME
    if starttime[:7] != stoptime[:7]:
        utn[utn > 20] = utn[utn > 20]-24
        uts[uts > 20] = uts[uts > 20]-24
    # Now, utn can be negative, negative valus mean yesterday

    # 5, initialize output
    pp = pd.DataFrame(np.ones([192, 9])*np.nan,
                      columns=['pbmlat1', 'mlatm1', 'ebmlat1',
                               'pbmlat2', 'mlatm2', 'ebmlat2',
                               'peak1', 'peak2', 'datetime'],
                      index=pd.MultiIndex.from_product(
                              [['N', 'S'], np.arange(BinMLT/2, 24, BinMLT)]))

    # 6, find parameters
    # 6.1, group data according to hemispheres
    for k00, k0 in enumerate(['N', 'S']):
        # select data in the specified hemisphere
        if k0 == 'N':
            mlat0, mlt0, ef0, ut0 = (MLatn, MLTn, energyfluxn, utn)
            # bvv = bv[0]
        if k0 == 'S':
            mlat0, mlt0, ef0, ut0 = (MLats, MLTs, energyfluxs, uts)
            # bvv = bv[1]

        # 6.2, calculate the cumulative energy flux in the hemisphere,
        # then find the lower and upper limits of the mlat in one bin
        idx = np.argsort(mlat0)
        mlat00, ef00 = mlat0[idx], ef0[idx]
        ef00 = ef00.cumsum()
        if ef00[-1] < LLEFTOTH:
            continue
        sf = 0.9
        fpt = ef00 < (0.5-sf/2)*ef00[-1]
        llmlat = mlat00[fpt][-1]
        fpt = ef00 > (0.5+sf/2)*ef00[-1]
        ulmlat = mlat00[fpt][0]

        # 6.3, group data according to MLT
        for k11, k1 in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
            peak1, peak2 = False, False
            # Select data in the specified MLT range
            lmlt, rmlt = k1-BinMLT/2, k1+BinMLT/2
            fp = (mlt0 >= lmlt) & (mlt0 <= rmlt)
            if np.sum(fp) <= LLP1:
                continue
            mlat1, ef1, ut1 = mlat0[fp], ef0[fp], ut0[fp]

            # jump to next loop when only part of the aurora is captured
            if mlat1.min() > llmlat or mlat1.max() < ulmlat:
                continue

            # 6.5, sort MLat
            idx = np.argsort(mlat1)
            mlat1, ef1, ut1 = mlat1[idx], ef1[idx], ut1[idx]

            # 6.6, bin energy flux
            # Is this necessary?
            mlat11 = np.arange(50+BinW/2, 90, BinW)
            ef11 = np.arange(50+BinW/2, 90, BinW)*np.nan
            for k22, k2 in enumerate(mlat11):
                lmlat, rmlat = k2-BinW/2, k2+BinW/2
                fpt = (mlat1 > lmlat) & (mlat1 < rmlat)
                if np.sum(fpt) > 0:
                    ef11[k22] = ef1[fpt].max()  # envelope of energy flux
            fpt = ~np.isnan(ef11)
            mlat11, ef11 = mlat11[fpt], ef11[fpt]

            # 6.7, Find maximum energy flux
            fpmax1 = np.argmax(ef11)
            mlatm1, efm1 = mlat11[fpmax1], ef11[fpmax1]
            if efm1 < LLMEF:
                continue

            # 6.8, Find boundaries
            # 1, energy flux < some value for some consecutive values
            # or 2, difference of Mlat > some value
            # poleward boundary
            efef = np.ones([len(ef11), EdgeN])*np.nan
            for k2 in range(EdgeN):
                efef[:, k2] = np.roll(ef11, -k2)
                efef[-(k2+1):, k2] = 0
            fpt = np.sum(efef < EdgeEF, axis=1) >= EdgeM
            idx2 = np.where(fpt)[0]
            idxpb = idx2[idx2 >= fpmax1].min()
            # equatorward boundary
            efef = np.ones([len(ef11), EdgeN])*np.nan
            for k2 in range(EdgeN):
                efef[:, k2] = np.roll(ef11, k2)
                efef[:k2+1, k2] = 0
            fpt = np.sum(efef < EdgeEF, axis=1) >= EdgeM
            idx2 = np.where(fpt)[0]
            idxeb = idx2[idx2 <= fpmax1].max()
            # equatorward and poleward boundary mlat
            ebmlat1, pbmlat1 = mlat11[idxeb], mlat11[idxpb]

            # 6.9, check the existence of the second auroral
            # If the numbers of (the residual pixels with energy flux > edgev)
            # > some value then, the second maximum exists
            mlat22, ef22 = mlat11.copy(), ef11.copy()
            fpt = (mlat11 >= ebmlat1) & (mlat11 <= pbmlat1)
            mlat22[fpt] = 0
            ef22[fpt] = 0
            if np.sum(ef22 > EdgeEF) > LLA2F*np.sum(fpt):
                # 6.10, Find the second maximum
                fpmax2 = np.argmax(ef22)
                mlatm2, efm2 = mlat22[fpmax2], ef22[fpmax2]

                # 6.11, Find the second boundary
                # poleward boundary
                efef = np.ones([len(ef22), EdgeN])*np.nan
                for k2 in range(EdgeN):
                    efef[:, k2] = np.roll(ef22, -k2)
                    efef[-(k2+1):, k2] = 0
                fpt = np.sum(efef < EdgeEF, axis=1) >= EdgeM
                idx2 = np.where(fpt)[0]
                idxpb = idx2[idx2 >= fpmax2].min()
                # equatorward boundary
                efef = np.ones([len(ef22), EdgeN])*np.nan
                for k2 in range(EdgeN):
                    efef[:, k2] = np.roll(ef22, k2)
                    efef[:k2+1, k2] = 0
                fpt = np.sum(efef < EdgeEF, axis=1) >= EdgeM
                idx2 = np.where(fpt)[0]
                idxeb = idx2[idx2 <= fpmax2].max()
                ebmlat2, pbmlat2 = mlat22[idxeb], mlat22[idxpb]

                # if the two auroral ovals are close, they are combined
                if (ebmlat2-pbmlat1 <= LLdMLat) and \
                   (ebmlat1-pbmlat2 <= LLdMLat):
                    ebmlat1 = np.min([ebmlat1, ebmlat2])
                    pbmlat1 = np.max([pbmlat1, pbmlat2])
                else:
                    fpt = (mlat1 > ebmlat2) & (mlat1 < pbmlat2)
                    if np.sum(fpt) > LLP2:
                        peak2 = True  # The second auroral peak exists
            fpt = (mlat1 > ebmlat1) & (mlat1 < pbmlat1)
            if np.sum(fpt) >= LLP2:
                peak1 = True  # The second auroral peak exists

            if peak1 and peak2:
                ebmlat2, pbmlat2, ebmlat1, pbmlat1 = np.sort(
                        [ebmlat2, pbmlat2, ebmlat1, pbmlat1])
                idx = np.argsort([mlatm1, mlatm2])
                mlatm2, mlatm1 = np.array([mlatm1, mlatm2])[idx]
                efm2, efm1 = np.array([efm1, efm2])[idx]
            elif peak1:
                ebmlat2, pbmlat2, ebmlat1, pbmlat1 = (
                        ebmlat1, np.nan, np.nan, pbmlat1)
                mlatm2, mlatm1 = np.nan, mlatm1
                efm2, efm1 = np.nan, efm1
            elif peak2:
                ebmlat2, pbmlat2, ebmlat1, pbmlat1 = (
                        ebmlat2, np.nan, np.nan, pbmlat2)
                mlatm2, mlatm1 = np.nan, mlatm2
                efm2, efm1 = np.nan, efm2
            else:
                (ebmlat2, pbmlat2, ebmlat1, pbmlat1, mlatm2, mlatm1,
                 efm2, efm1) = np.ones(8)*np.nan
            if peak1 or peak2:
                utt = np.median(ut1[fpt])
            else:
                utt = np.nan
            pp.loc[(k0, k1), :] = [
                    pbmlat1, mlatm1, ebmlat1, pbmlat2, mlatm2, ebmlat2,
                    efm1, efm2, utt]
    # 7, remove defective points
    for k00, k0 in enumerate(('N', 'S')):
        # IF pbmlat1 or ebmlat2 is `mlatt` degrees greater or less than the
        # median value of neighboring `points` points, the corresponding
        # points are removed
        points = 25
        mlatt = 5
        pp0 = pp.loc[(k0, slice(None)), ['pbmlat1', 'ebmlat2']].values
        ppt = []  # initialize
        for k1 in range(points):
            ppt.append(np.roll(pp0, k1-int((points-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp1 = np.nanmedian(pp1, axis=2)
        fpt = (np.abs(pp1-pp0) > mlatt).any(axis=1)
        pp2 = pp.loc[(k0, slice(None))].values
        pp2[fpt] = np.nan
        pp.loc[(k0, slice(None))] = pp2

    # 8, smooth result
    for k00, k0 in enumerate(('N', 'S')):
        pp0 = pp.loc[(k0, slice(None))].values
        ppt = []  # initialize
        for k1 in range(smoothps):
            ppt.append(np.roll(pp0, k1-int((smoothps-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp0 = np.median(pp1, axis=2)
        pp.loc[(k0, slice(None))] = pp0

    # 9, If the total number of MLTs < some vlaue, we think there is no auroral
    for k00, k0 in enumerate(('N', 'S')):
        pp0 = pp.loc[(k0, slice(None)), 'pbmlat1'].values
        if np.sum(~np.isnan(pp0)) <= LLN:
            pp.loc[(k0, slice(None))] = np.nan

    # 10, calculate datetime
    yyyy = int(ssusi['YEAR'][:])
    ddd = int(ssusi['DOY'][:])
    stime = pd.Timestamp(yyyy, 1, 1, 0, 0, 0)
    dt1 = pd.Timedelta(ddd-1, 'D')
    dt2 = pd.TimedeltaIndex(pp['datetime'].values, 'h')
    pp['datetime'] = stime+dt1+dt2
    return pp


def find_parameters_one_file_5(fn):
    '''
    Set flags at the MLats that have 5%, 10%, ..., 95% of the
    total energy flux save MLats and cumulative energy flux
    '''
    # Read data
    ssusi = Dataset(fn)

    # Read variables that will be used
    # orbitn = ssusi.STARTING_ORBIT_NUMBER  # orbit number
    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE  # value in a no data bin
    # magnetic latitude
    MLat = np.array(ssusi['LATITUDE_GEOMAGNETIC_GRID_MAP'])
    MLT = np.array(ssusi['MLT_GRID_MAP'])  # magnetic local time
    # electron energy flux
    energyfluxn = np.array(ssusi['ENERGY_FLUX_NORTH_MAP'])
    energyfluxs = np.array(ssusi['ENERGY_FLUX_SOUTH_MAP'])
    utn = np.array(ssusi['UT_N'])
    uts = np.array(ssusi['UT_S'])
    # bv = np.array(ssusi['ELECTRON_ENERGY_FLUX_THRESHOLDS']) # boundary values

    # Apply median filter to the energy flux to remove impulse noise
    kernel_size_own = 3
    energyfluxn = signal.medfilt2d(energyfluxn, kernel_size_own)
    energyfluxs = signal.medfilt2d(energyfluxs, kernel_size_own)

    # exclude pixels without data
    fpn = (utn != nodatavalue)
    MLatn, MLTn, energyfluxn, utn = (
            k[fpn] for k in (MLat, MLT, energyfluxn, utn))
    fps = uts != nodatavalue
    MLats, MLTs, energyfluxs, uts = (
            k[fps] for k in (MLat, MLT, energyfluxs, uts))

    # Change UT if part of the previous day is included.
    starttime = ssusi.STARTING_TIME   # format: yyyydddhhmmss
    stoptime = ssusi.STOPPING_TIME
    yyyy = int(stoptime[:4])
    ddd = int(stoptime[4:7])
    date = pd.Timestamp(yyyy, 1, 1)+pd.Timedelta(ddd-1, 'D')
    if starttime[:7] != stoptime[:7]:
        utn[utn > 20] = utn[utn > 20]-24
        uts[uts > 20] = uts[uts > 20]-24
    # Now, utn can be negative, negative valus mean yesterday

    # Initialize output
    pp1 = np.ones([np.int(2*24/BinMLT), 19])*np.nan  # MLats
    pp2 = pd.DataFrame(
            index=np.arange(int(2*24/BinMLT)),
            columns=('ns', 'mlt', 'datetime', 'eftot'))

    # group data according to hemispheres
    for k00, k0 in enumerate(['N', 'S']):
        # select data in the specified hemisphere
        if k0 == 'N':
            mlat0, mlt0, ef0, ut0 = (MLatn, MLTn, energyfluxn, utn)
            # bvv = bv[0]
        if k0 == 'S':
            mlat0, mlt0, ef0, ut0 = (MLats, MLTs, energyfluxs, uts)
            # bvv = bv[1]

        # Calculate the cumulative energy flux in the hemisphere,
        # then find the lower and upper limits of mlat
        idx = np.argsort(mlat0)
        mlat00, ef00 = mlat0[idx], ef0[idx]
        ef00 = ef00.cumsum()
        if ef00[-1] < LLEFTOTH:
            continue
        sf = 0.9
        fpt = ef00 < (0.5-sf/2)*ef00[-1]
        llmlat = mlat00[fpt][-1]
        fpt = ef00 > (0.5+sf/2)*ef00[-1]
        ulmlat = mlat00[fpt][0]

        # group data according to MLT
        for k11, k1 in enumerate(np.arange(BinMLT/2, 24, BinMLT)):

            idpp = k00*int(24/BinMLT)+k11
            # Select data in the specified MLT range
            lmlt, rmlt = k1-BinMLT/2, k1+BinMLT/2
            fp = (mlt0 >= lmlt) & (mlt0 <= rmlt)
            if np.sum(fp) <= LLP1:  # lower limit of pixels in one bin
                continue
            mlat1, ef1, ut1 = mlat0[fp], ef0[fp], ut0[fp]
            if mlat1.min() > llmlat or mlat1.max() < ulmlat:
                continue

            # sort MLat
            idx = np.argsort(mlat1)
            mlat1, ef1, ut1 = mlat1[idx], ef1[idx], ut1[idx]

            # calculate the integral energy flux
            efcum = ef1.cumsum()
            eftot = efcum[-1]
            if eftot < LLEFTOT:
                continue

            # Mlats of 5%, 10%... energy flux
            pp1t = np.ones(19)*np.nan
            utt = np.ones(19)*np.nan
            for k22, k2 in enumerate(np.arange(0.05, 0.96, 0.05)):
                fp = (efcum > eftot*k2)
                pp1t[k22] = mlat1[fp][0]
                utt[k22] = ut1[fp][0]
            pp1[idpp, :] = pp1t

            # The UT
            pp2.iloc[idpp, :] = [
                    k0, k1, date+pd.Timedelta(np.median(utt), 'h'), eftot]
    return pp1, pp2


def test_find_parameters_one_file_12(fn, method=1):
    '''
    Input:
        fn        : ssusi file name
    '''
    if method == 1:
        pp = find_parameters_one_file_1(fn)
    if method == 2:
        pp = find_parameters_one_file_2(fn)
    plt.figure(figsize=(8.7, 4.4))
    for k00, k0 in enumerate(['N', 'S']):
        ax = plt.subplot(1, 2, k00+1, polar=True)
        image_energy_flux_one_file(ax, fn, k0, vmin=0, vmax=5, s=1)
        eee = pp.loc[(k0, slice(None))]
        r1 = 90-eee['pbmlat']
        r2 = 90-eee['ebmlat']
        theta = eee.index/12*np.pi
        ax.scatter(theta, r1, c='b', s=5)
        ax.scatter(theta, r2, c='b', s=5)
        plt.title(k0+'H')
    plt.tight_layout()
    return ax


def test_find_parameters_one_file_34(fn, method=3):
    '''
    Input:
        fn        : ssusi file name
    '''
    if method == 3:
        pp = find_parameters_one_file_3(fn)
    if method == 4:
        pp = find_parameters_one_file_4(fn)
    plt.figure(figsize=(8.7, 4.4))
    for k00, k0 in enumerate(['N', 'S']):
        ax = plt.subplot(1, 2, k00+1, polar=True)
        image_energy_flux_one_file(ax, fn, k0, vmin=0, vmax=5, s=1)
        eee = pp.loc[(k0, slice(None))]
        r1 = 90-eee['pbmlat1']
        r2 = 90-eee['ebmlat2']
        r3 = 90-eee['ebmlat1']
        r4 = 90-eee['pbmlat2']
        theta = eee.index/12*np.pi
        ax.scatter(theta, r1, c='b', s=5)
        ax.scatter(theta, r2, c='b', s=5)
        ax.scatter(theta, r3, c='b', s=5)
        ax.scatter(theta, r4, c='b', s=5)
        plt.title(k0+'H')
    plt.tight_layout()
    return ax


def test_find_parameters_one_file_5(fn):
    '''
    Input:
        fn        : ssusi file name
    '''
    pp1, pp2 = find_parameters_one_file_5(fn)
    plt.figure(figsize=(8.7, 4.4))
    for k00, k0 in enumerate(['N', 'S']):
        ax = plt.subplot(1, 2, k00+1, polar=True)
        image_energy_flux_one_file(ax, fn, k0, vmin=0, vmax=5, s=1)
        fp = (pp2.ns == k0).values
        for k11, k1 in enumerate(range(19)):
            r = 90-pp1[fp, k11]
            theta = (pp2.loc[fp, 'mlt']/12*np.pi).values.astype(float)
            ax.scatter(theta, r, c='b', s=4, alpha=0.5)
        plt.title(k0+'H')
    plt.tight_layout()
    return ax


def find_parameters_2013_12(method=1):
    import os
    import fnmatch
    import omni
    fna = os.listdir(datapath)

    # Read IMF By, Bz and AE in 2013
    print('Begin to read solar wind speed, IMF, AE and Dst')
    imfae = omni.get_omni(
            bdate='2013-1-1', edate='2014-1-1',
            variables=['Bym', 'Bzm', 'AE', 'V'], res='1m')
    imfae['Bt'] = np.sqrt(imfae['Bym']**2 + imfae['Bzm']**2)
    imfae['nwcf'] = imfae['V']**(4/3) * imfae['Bt']**(2/3) * \
        ((1-imfae['Bzm']/imfae['Bt'])/2)**(4/3)/100.0
    imfae['EKL'] = 0.5*imfae['V']*(imfae['Bt']-imfae['Bzm'])/1e6
    dst = omni.get_omni(bdate='2013-1-1', edate='2014-1-1',
                        variables=['DST'], res='1h')
    dst = dst.reindex(imfae.index, method='nearest')
    imfae['Dst'] = dst.values
    print('End of reading solar wind speed, IMF, AE and Dst')

    for k0 in range(16, 19):  # satellite
        pp = []
        savefn = 'F{:d}_2013.dat'.format(k0)  # file to save data
        for k1 in range(1, 13):  # month
            print('Begin month {:02d}'.format(k1))
            for k2 in range(1, 32):  # day
                # 00 means full orbit?
                fn = fnmatch.filter(
                        fna,
                        '*F{:d}*2013{:02d}{:02d}*00_DF.NC'.format(k0, k1, k2))
                if not fn:  # empty fn
                    continue
                for k33, k3 in enumerate(fn):
                    print('Processing '
                          'satellite {:s}, orbit {:s}'.
                          format(k3[36:39], k3[-14:-9]))
                    if method == 1:
                        pp0 = find_parameters_one_file_1(datapath+k3)
                    if method == 2:
                        pp0 = find_parameters_one_file_2(datapath+k3)
                    pp0 = pp0.loc[~(pp0['pbmlat'].isnull().values)]
                    if not pp0.empty:
                        pp.append(pp0)

                        fig = plt.figure(figsize=(8.7, 4.4))
                        for k44, k4 in enumerate(['N', 'S']):
                            ax = plt.subplot(1, 2, k44+1, polar=True)
                            image_energy_flux_one_file(
                                    ax, datapath+k3, k4, s=1, vmin=0, vmax=10)
                            ppt = pp0.loc[(k4, slice(None))]
                            r1 = 90-ppt['pbmlat']
                            r2 = 90-ppt['ebmlat']
                            theta = ppt.index/12*np.pi
                            ax.scatter(theta, r1, s=1, c='r', alpha=0.8)
                            ax.scatter(theta, r2, s=1, c='r', alpha=0.8)
                            ax.set_title(k4+'H')
                        plt.tight_layout()
                        plt.savefig(
                                savepath +
                                'method{:d}/F{:d}_2013{:02d}{:02d}_{:s}.png'.
                                format(method, k0, k1, k2, k3[-14:-9]))
                        plt.close(fig)
                    else:
                        print('    No parameters found in this file')
            print('End of month {:02d}'.format(k1))
        pp = pd.concat(pp, axis=0)
        imfaet = imfae.reindex(pp['datetime'], method='ffill')
        pp['Bym'] = imfaet['Bym'].values
        pp['Bzm'] = imfaet['Bzm'].values
        pp['AE'] = imfaet['AE'].values
        pp['Dst'] = imfaet['Dst'].values
        pp['nwcf'] = imfaet['nwcf'].values
        pp['EKL'] = imfaet['EKL'].values
        pp.to_pickle(savepath+'method{:d}/'.format(method)+savefn)
    return


def find_parameters_2013_34(method=3):
    import os
    import fnmatch
    import omni
    fna = os.listdir(datapath)
    # Read IMF By, Bz and AE in 2013
    print('Begin to read solar wind speed, IMF, AE and Dst')
    imfae = omni.get_omni(
            bdate='2013-1-1', edate='2014-1-1',
            variables=['Bym', 'Bzm', 'AE', 'V'], res='1m')
    imfae['Bt'] = np.sqrt(imfae['Bym']**2 + imfae['Bzm']**2)
    imfae['nwcf'] = imfae['V']**(4/3) * \
        imfae['Bt']**(2/3) * \
        ((1-imfae['Bzm']/imfae['Bt'])/2)**(4/3)/100.0
    imfae['EKL'] = 0.5*imfae['V']*(imfae['Bt']-imfae['Bzm'])/1e6
    dst = omni.get_omni(bdate='2013-1-1', edate='2014-1-1',
                        variables=['DST'], res='1h')
    dst = dst.reindex(imfae.index, method='nearest')
    imfae['Dst'] = dst.values
    print('End of reading solar wind speed, IMF, AE and Dst')
    for k0 in range(16, 19):  # satellite
        pp = []
        savefn = 'F{:d}_2013.dat'.format(k0)  # file to save data
        for k1 in range(1, 13):  # month
            print('Begin month {:02d}'.format(k1))
            for k2 in range(1, 32):  # day
                # 00 means full orbit?
                fn = fnmatch.filter(
                        fna,
                        '*F{:d}*2013{:02d}{:02d}*00_DF.NC'.format(k0, k1, k2))
                if not fn:  # empty fn
                    continue
                for k33, k3 in enumerate(fn):
                    print('Processing '
                          'satellite {:s}, orbit {:s}'.
                          format(k3[36:39], k3[-14:-9]))
                    if method == 3:
                        pp0 = find_parameters_one_file_3(datapath+k3)
                    if method == 4:
                        pp0 = find_parameters_one_file_4(datapath+k3)
                    pp0 = pp0.loc[~(pp0['pbmlat1'].isnull().values)]
                    if not pp0.empty:
                        pp.append(pp0)

                        fig = plt.figure(figsize=(8.7, 4.4))
                        for k44, k4 in enumerate(['N', 'S']):
                            ax = plt.subplot(1, 2, k44+1, polar=True)
                            image_energy_flux_one_file(
                                    ax, datapath+k3, k4, s=1, vmin=0, vmax=10)
                            ppt = pp0.loc[(k4, slice(None))]
                            r1 = 90-ppt['pbmlat1']
                            r2 = 90-ppt['ebmlat2']
                            r3 = 90-ppt['ebmlat1']
                            r4 = 90-ppt['pbmlat2']
                            theta = ppt.index/12*np.pi
                            ax.scatter(theta, r1, s=1, c='r', alpha=0.8)
                            ax.scatter(theta, r2, s=1, c='r', alpha=0.8)
                            ax.scatter(theta, r3, s=1, c='r', alpha=0.8)
                            ax.scatter(theta, r4, s=1, c='r', alpha=0.8)
                            ax.set_title(k4+'H')
                        plt.tight_layout()
                        plt.savefig(
                                savepath +
                                'method{:d}/F{:d}_2013{:02d}{:02d}_{:s}.png'.
                                format(method, k0, k1, k2, k3[-14:-9]))
                        plt.close(fig)
                    else:
                        print('    No parameters found in this file')
            print('End of month {:02d}'.format(k1))
        pp = pd.concat(pp, axis=0)
        imfaet = imfae.reindex(pp['datetime'], method='ffill')
        pp['Bym'] = imfaet['Bym'].values
        pp['Bzm'] = imfaet['Bzm'].values
        pp['AE'] = imfaet['AE'].values
        pp['Dst'] = imfaet['Dst'].values
        pp['nwcf'] = imfaet['nwcf'].values
        pp['EKL'] = imfaet['EKL'].values
        pp.to_pickle(savepath+'method{:d}/'.format(method)+savefn)
    return


def parameters_fit(
        ax, pp, method=1, ns='N', hour=19, shour=1/8, vary='width', varx='AE'):
    '''
    linear fit of the parameters
    Input:
        ax          = axis
        pp          = parameter
        method      = 1,2,3,4
        ns          = 'N' or 'S'
        hour        = 0-23
        shour       = 1/8, 3/8, 5/8, 7/8
        vary        = 'width', 'efm', 'mlatm', 'pe', 'ee'
        varx        = 'Bym', 'Bzm', 'AE', 'Dst', 'nwcf', 'EKL'
    '''
    from scipy import stats
    pp0 = pp.loc[ns]
    if hour+shour not in pp0.index:
        return np.nan, np.nan, np.nan
    pp1 = pp.loc[ns].loc[hour+shour]
    vx = pp1[varx]
    if method == 1 or method == 2:
        cc = True
    else:
        cc = False

    # y axis
    if vary == 'efm':
        vy = pp1['peak' if cc else 'peak1']
        # yt = r'Max Energy Flux (ergs/s/$cm^2$)'
    elif vary == 'mlatm':
        vy = pp1['mlatm' if cc else 'mlatm1']
        # yt = 'Peak Mlat'
    elif vary == 'width':
        vy = pp1['pbmlat' if cc else 'pbmlat1'] - \
             pp1['ebmlat' if cc else 'ebmlat2']
        # yt = 'Width (degree)'
    elif vary == 'pe':
        vy = pp1['pbmlat' if cc else 'pbmlat1']
        # yt = 'Poleward MLat (degree)'
    elif vary == 'ee':
        vy = pp1['ebmlat' if cc else 'ebmlat2']
        # yt = 'Equatorward MLat (degree)'

    # median
    if varx == 'AE':
        AEMAX = 400
        xr = np.linspace(0, AEMAX, 21)
        fp = vx < AEMAX
        vx, vy = vx[fp], vy[fp]
        # xt = 'AE'
    elif varx == 'Dst':
        DSTMAX = 0
        xr = np.linspace(-40, DSTMAX, 21)
        fp = vx < DSTMAX
        vx, vy = vx[fp], vy[fp]
        # xt = 'Dst'
    elif varx == 'nwcf':
        xr = np.linspace(0, 100, 21)
        # xt = 'Newell Coupling Function'
    elif varx == 'EKL':
        xr = np.linspace(0, 0.004, 21)
        # xt = 'Kan Lee E'
    elif varx == 'Bym':
        xr = np.linspace(-8, 8, 21)
        # xt = 'GSM By'
    elif varx == 'Bzm':
        xr = np.linspace(-8, 8, 21)
        # xt = 'GSM Bz'

    if np.size(vx) < 100:
        return np.nan, np.nan, np.nan
    plt.sca(ax)
    plt.scatter(vx, vy, s=1, c='b')
    #  yr=xr[:-1] * 1
    #  for k00, k0 in enumerate(xr[:-1]):
    #      fp = (vx>xr[k00]) & (vx<xr[k00+1])
    #      yr[k00] = np.nanmedian(vy.loc[fp])
    #  plt.plot(xr[:-1], yr, 'r')

    # linear fit
    slop, intercept, r_value, p_value, stderr = stats.linregress(vx, vy)
    plt.plot(xr, slop*xr+intercept, 'k')
    plt.text(0.1, 0.7, 'slop: {:.3g}'
             '\nintercept: {:.3g}'
             '\ncorrelation: {:.3f}'.format(slop, intercept, r_value),
             transform=plt.gca().transAxes, fontsize=8)
    # set coordinates
    # plt.xlabel(xt)
    # plt.ylabel(yt)
    plt.xlim(xr.min(), xr.max())
    return slop, intercept, r_value


def parameters_fit_all(method=1):
    fig = plt.figure()
    ax = plt.subplot()
    for k0 in ['F16', 'F17', 'F18']:
        pp = pd.read_pickle(
                savepath+'method{:d}/'.format(method)+k0+'_2013.dat')
        for k1 in ['N', 'S']:
            fit_parameter = pd.DataFrame()
            for k2 in range(24):
                for k3 in [1/8, 3/8, 5/8, 7/8]:
                    if k2+k3 not in pp.loc[k1].index:
                        print('Empty ', k0, ' ', k1, ' ', k2, ' ', k3)
                        continue
                    fig, ax = plt.subplots(
                            5, 2, sharex='col', sharey='row', figsize=(6, 10))
                    for k44, k4 in enumerate(
                            ['ee', 'pe', 'width', 'mlatm', 'efm']):
                        axt = ax[k44]
                        for k55, k5 in enumerate(['AE', 'Dst']):
                            axtt = axt[k55]
                            slop, intercept, cf = parameters_fit(
                                    axtt, pp, ns=k1, hour=k2, shour=k3,
                                    vary=k4, varx=k5, method=method)
                            fff = pd.DataFrame(
                                    [[k4, k5, k2+k3, slop, intercept, cf]],
                                    columns=('y', 'x', 't', 'slop',
                                             'intercept', 'correlation'))
                            fit_parameter = fit_parameter.append(
                                    fff, ignore_index=True)

                    ax[0, 0].set_ylabel('EE')
                    ax[0, 0].set_ylim(60, 80)
                    ax[1, 0].set_ylabel('PE')
                    ax[1, 0].set_ylim(60, 80)
                    ax[2, 0].set_ylabel('Width')
                    ax[2, 0].set_ylim(0, 10)
                    ax[3, 0].set_ylabel('Peak Mlat')
                    ax[3, 0].set_ylim(60, 80)
                    ax[4, 0].set_ylabel('Peak EF')
                    ax[4, 0].set_ylim(0, 20)
                    ax[-1, 0].set_xlabel('AE')
                    ax[-1, 0].set_xlim(0, 400)
                    ax[-1, 1].set_xlabel('Dst')
                    ax[-1, 1].set_xlim(-40, 0)
                    # ax[-1, 2].set_xlabel('NWCF')
                    # ax[-1, 2].set_xlim(0,100)
                    # ax[-1, 3].set_xlabel('Bz')
                    # ax[-1, 3].set_xlim(-8,8)
                    # ax[-1, 4].set_xlabel('By')
                    # ax[-1, 4].set_xlim(-8,8)
                    plt.subplots_adjust(top=0.95, bottom=0.05, wspace=0.1)
                    fig.savefig(
                            savepath+'method{:d}/'.format(method) + 'fit/' +
                            '{:s}_{:s}_{:02d}_{:d}'.format(
                                k0, k1, k2, int((k3*8+1)/2), k4, k5))
                    plt.close(fig)
            pd.to_pickle(fit_parameter, savepath +
                         'method{:d}/fit/fit_parameter_{:s}_{:s}.dat'.
                         format(method, k0, k1))


def compare_methods():
    plt.close('all')
    sat = 'F17'
    ns = 'S'
    x = 'AE'
    y = ('ee', 'pe', 'width', 'efm', 'mlatm')
    yl = ('EE', 'PE', 'Width', 'Peak EF', 'Peak Mlat')
    fig, ax = plt.subplots(5, 1, sharex=True, figsize=(4.7, 8.5))
    for k0 in range(4):
        for k11, k1 in enumerate(y):
            plt.sca(ax[k11])
            p1 = pd.read_pickle(
                    '/home/guod/WD4T/work4/method{:d}/fit/'
                    'fit_parameter_{:s}_{:s}.dat'.format(
                        k0+1, sat, ns))
            fp = (p1.y == k1) & (p1.x == x)
            plt.plot(p1[fp].t, p1[fp].correlation)
            plt.ylabel(yl[k11])
    plt.legend(('80% of total EF', 'method2', 'Edge: 0.2*EF',
                'Edge: 0.5 ergs/s/$cm^2$'), fontsize='x-small')
    ax[0].set_title('Correlation coefficient')
    plt.xlim(0, 24)
    plt.xticks(np.arange(0, 25, 4))
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    plt.close('all')
    fn = (datapath+'PS.APL_V0105S024CE0018_SC.U_DI.A_GP.F16'
          '-SSUSI_PA.APL-EDR-AURORA_DD.20130216_SN.48154-00_DF.NC')
    # fn = (datapath+'PS.APL_V0105S024CE0018_SC.U_DI.A_GP.F16'
    #       '-SSUSI_PA.APL-EDR-AURORA_DD.20130102_SN.47516-00_DF.NC')
    # test_find_parameters_one_file_12(fn,1)
    # test_find_parameters_one_file_12(fn,2)
    # test_find_parameters_one_file_34(fn, 3)
    # test_find_parameters_one_file_34(fn, 4)
    test_find_parameters_one_file_5(fn)
    plt.show()

    # find_parameters_2013_12(method=2)
    # find_parameters_2013_34(method=4)
    # find_parameters_2013_34(method=3)
    # find_parameters_2013_12(method=1)

    # plt.close('all')
    # sat = 'F16'
    # pp = pd.read_pickle(savepath+'method1/'+sat+'_2013.dat')
    # ax = plt.subplot()
    # parameters_fit_12(ax, pp, hour=2, shour=1/8, vary='pe', varx='EKL')
    # plt.show()

    # parameters_fit_all(method=3)
    # parameters_fit_all(method=1)
    # parameters_fit_all(method=2)
    # parameters_fit_all(method=4)

    # compare_methods()

    # plt.close('all')
    # fig, ax = plt.subplots(1, 3, subplot_kw=dict(polar=True), figsize=(13,5))
    # fn = (datapath+'PS.APL_V0105S024CE0018_SC.U_DI.A_GP.F16'
    #       '-SSUSI_PA.APL-EDR-AURORA_DD.20130102_SN.47516-00_DF.NC')
    # image_energy_flux_one_file(
    #     ax[0], fn, ns='S', vmin=0, vmax=5, kernel_size=1)
    # image_energy_flux_one_file(
    #     ax[1], fn, ns='S', vmin=0, vmax=5, kernel_size=3)
    # image_energy_flux_one_file(
    #     ax[2], fn, ns='S', vmin=0, vmax=5, kernel_size=5)
    # plt.show()
    print('----------------------------------------------------------------')
    print('Remaining problems:')
    print('  1, Aurora with complicated structure')
    print('  2, Which IMF values to use?')
    print('  3, Different total energy, different proportion')
    print('----------------------------------------------------------------')
