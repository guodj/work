# Auroral precipation model using SSUSI data
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import myfunctions as mf
from scipy.optimize import curve_fit
import matplotlib.pylab as plt

def image_energy_flux_one_file(ax, fn, ns, vmin=None, vmax=None, s=2, alpha=0.8):
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
    Note: for both north and south poles, the parameter 'ns' in set_polar is 'N'
    '''
    ssusi = Dataset(fn)
    MLat = np.array(ssusi['LATITUDE_GEOMAGNETIC_GRID_MAP'])
    MLT = np.array(ssusi['MLT_GRID_MAP'])
    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE # value in a no data bin
    if ns=='N':
        var='ENERGY_FLUX_NORTH_MAP'
        ut = ssusi['UT_N'][:]
    if ns=='S':
        var='ENERGY_FLUX_SOUTH_MAP'
        ut = ssusi['UT_S'][:]
    variable = np.array(ssusi[var])
    # Only include data points inside the boundary
    fp = (ut != nodatavalue)
    MLat, MLT, variable = (k[fp] for k in (MLat, MLT, variable))
    # plot variable
    r = 90-MLat
    theta = MLT/12*np.pi
    hs = ax.scatter(theta, r, s=s, c=variable, vmin=vmin,
                    vmax=vmax, alpha=alpha)
    # Set polar coordinates
    mf.set_polar(
            ax, ns='N', boundinglat=50,
            dlat=10, dlon=90, useLT=True)
    ax.set_title(ssusi[var].TITLE)
    ax.text(0.8, -0.1,
            'YEAR: '+str(np.int(ssusi['YEAR'][:])) +
            '\nDOY: '+str(np.int(ssusi['DOY'][:])) +
            '\nORBIT: '+str(np.int(ssusi.STARTING_ORBIT_NUMBER)),
            transform=ax.transAxes,
            fontsize=8)
    #cb = plt.colorbar(hs, ax=ax, pad=0.1)
    #cb.set_label(ssusi[var].UNITS)
    return ax, hs


def find_parameters_one_file(fn):
    '''
    Get the peak value, peak latitude, and full width of the aurora
    input:
        fn        : ssusi aurora file name
    output:
        pp        : the parameters
    '''
    # 1, Read data
    ssusi = Dataset(fn)

    #--------------------------------------------------
    #--------------adjustable parameters---------------
    #
    BinMLT = 0.25 # MLT range in one bin
    LLP1 = 200 # Lower limit of pixels in one bin
    LLP2 = 20 # minimum data points in auroral oval
    # Lower limit of residual MLAT bins that has energy > edgev
    # (relative to the number of the first auroral MLAT bins)
    Scale1 = 0.5
    Scale2 = 50  # Lower limit of the maximum energy flux (Scale2*0.2)
    BinW = 1*float(ssusi['PIXELSIZE_GEOMAGNETIC_LATITUDE'][:]) # size of the mlat bin window
    LLN =10 # Lower limit of the maximums and the final boundary numbers
    # Increase these values to get a narrower boundaries
    EdgeF = 0.2 # scale factor of max energy flux at boundaries
    # Increase these values to get a wider boundaries
    EdgeN = 3 # 'EdgeM' of 'EdgeN' consecutive data points are smaller than boundary values
    EdgeM = 3
    # Lower limit of mlat difference of two auroral ovals
    LLdMLat = 5*float(ssusi['PIXELSIZE_GEOMAGNETIC_LATITUDE'][:])
    #--------------------------------------------------
    #--------------------------------------------------

    # 2, Read variables that will be used
    orbitn = ssusi.STARTING_ORBIT_NUMBER # orbit number
    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE # value in a no data bin
    MLat = np.array(ssusi['LATITUDE_GEOMAGNETIC_GRID_MAP']) # magnetic latitude
    MLT = np.array(ssusi['MLT_GRID_MAP']) # magnetic local time
    energyfluxn = np.array(ssusi['ENERGY_FLUX_NORTH_MAP']) # electron energy flux
    energyfluxs = np.array(ssusi['ENERGY_FLUX_SOUTH_MAP'])
    utn = np.array(ssusi['UT_N'])
    uts = np.array(ssusi['UT_S'])
    bv = ssusi['ELECTRON_ENERGY_FLUX_THRESHOLDS'][:] # boundary values

    # 3, exclude pixels without data
    fpn = utn != nodatavalue
    MLatn, MLTn, energyfluxn, utn = (k[fpn] for k in (MLat, MLT, energyfluxn, utn))
    fps = uts != nodatavalue
    MLats, MLTs, energyfluxs, uts = (k[fps] for k in (MLat, MLT, energyfluxs, uts))

    # 4, Change UT if part of the previous day is included.
    starttime = ssusi.STARTING_TIME   # format: yyyydddhhmmss
    stoptime = ssusi.STOPPING_TIME
    if starttime[:7] != stoptime[:7]:
        utn[utn>20] = utn[utn>20]-24  # Now, utn can be negative, negative valus mean yesterday
        uts[uts>20] = uts[uts>20]-24

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
            bvv = bv[0]
        if k0 == 'S':
            mlat0, mlt0, ef0, ut0 = (MLats, MLTs, energyfluxs, uts)
            bvv = bv[1]

        # 6.2, Find boundary values for the whole hemisphere.
        maxh = []
        for k11, k1 in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
            # Select data of required MLT range
            lmlt, rmlt = k1-BinMLT/2, k1+BinMLT/2
            fp = (mlt0>=lmlt) & (mlt0<=rmlt)
            if np.sum(fp)<=LLP1:
                continue
            maxv = np.max(ef0[fp])
            if maxv > bvv*Scale2:
                maxh.append(maxv)
        if len(maxh)<LLN:
            continue
        else:
            edgev = EdgeF*np.median(maxh) # the edge energy flux of the hemisphere

        # 6.3, group data according to MLT
        for k11, k1 in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
            peak1, peak2 = False, False
            # Select data in the specified MLT range
            lmlt, rmlt = k1-BinMLT/2, k1+BinMLT/2
            fp = (mlt0>=lmlt) & (mlt0<=rmlt)
            if np.sum(fp)<=LLP1:
                #print(k0, k1, np.sum(fp))
                continue
            mlat1, ef1, ut1 = mlat0[fp], ef0[fp], ut0[fp]

            # 6.4, sort MLat
            idx = np.argsort(mlat1)
            mlat1, ef1, ut1 = mlat1[idx], ef1[idx], ut1[idx]

            # 6.5, bin energy flux
            mlat11 = np.arange(50+BinW/2, 90, BinW)
            ef11 = np.arange(50+BinW/2, 90, BinW)*np.nan
            for k22, k2 in enumerate(mlat11):
                lmlat, rmlat = k2-BinW/2, k2+BinW/2
                fpt = (mlat1>lmlat) & (mlat1<rmlat)
                if np.sum(fpt)>0:
                    ef11[k22] = ef1[fpt].max()  # envelope of energy flux
            fpt = ~np.isnan(ef11)
            mlat11, ef11 = mlat11[fpt], ef11[fpt]

            # 6.6, Find maximum energy flux
            fpmax1 = np.argmax(ef11)
            mlatm1, efm1 = mlat11[fpmax1], ef11[fpmax1]
            if efm1 < bvv*Scale2:
                continue

            # 6.7, Find boundaries
            # 1, energy flux < some value for some consecutive values
            # or 2, difference of Mlat > some value
            # poleward boundary
            efef = np.ones([len(ef11), EdgeN])*np.nan
            for k2 in range(EdgeN):
                efef[:, k2] = np.roll(ef11, -k2)
                efef[-(k2+1):, k2] = 0
            fpt = np.sum(efef<edgev, axis=1) >= EdgeM
            idx2 = np.where(fpt)[0]
            idxpb = idx2[idx2>=fpmax1].min()
            # equatorward boundary
            efef = np.ones([len(ef11), EdgeN])*np.nan
            for k2 in range(EdgeN):
                efef[:, k2] = np.roll(ef11, k2)
                efef[:k2+1, k2] = 0
            fpt = np.sum(efef<edgev, axis=1) >= EdgeM
            idx2 = np.where(fpt)[0]
            idxeb = idx2[idx2<=fpmax1].max()
            # equatorward and poleward boundary mlat
            ebmlat1, pbmlat1 = mlat11[idxeb], mlat11[idxpb]

            # 6.8, check the existence of the second auroral
            # If the numbers of (the residual pixels with energy flux > edgev) > some value
            # then, the second maximum exists
            mlat22, ef22 = mlat11*1, ef11*1  # *1, so 22 and 11 points to different object
            fpt = (mlat11>=ebmlat1) & (mlat11<=pbmlat1)
            mlat22[fpt] = 0
            ef22[fpt] = 0
            if np.sum(ef22>edgev)>Scale1*np.sum(fpt):
                # 6.9, Find the second maximum
                fpmax2 = np.argmax(ef22)
                mlatm2, efm2 = mlat22[fpmax2], ef22[fpmax2]

                # 6.10, Find the second boundary
                # poleward boundary
                efef = np.ones([len(ef22), EdgeN])*np.nan
                for k2 in range(EdgeN):
                    efef[:, k2] = np.roll(ef22, -k2)
                    efef[-(k2+1):, k2] = 0
                fpt = np.sum(efef<edgev, axis=1) >= EdgeM
                idx2 = np.where(fpt)[0]
                idxpb = idx2[idx2>=fpmax2].min()
                # equatorward boundary
                efef = np.ones([len(ef22), EdgeN])*np.nan
                for k2 in range(EdgeN):
                    efef[:, k2] = np.roll(ef22, k2)
                    efef[:k2+1, k2] = 0
                fpt = np.sum(efef<edgev, axis=1) >= EdgeM
                idx2 = np.where(fpt)[0]
                idxeb = idx2[idx2<=fpmax2].max()
                ebmlat2, pbmlat2 = mlat22[idxeb], mlat22[idxpb]

                # if the two auroral ovals are close, they are combined
                if (ebmlat2-pbmlat1<=LLdMLat) and (ebmlat1-pbmlat2<=LLdMLat):
                    ebmlat1 = np.min([ebmlat1, ebmlat2])
                    pbmlat1 = np.max([pbmlat1, pbmlat2])
                else:
                    fpt = (mlat1>ebmlat2) & (mlat1<pbmlat2)
                    if np.sum(fpt)>LLP2:
                        peak2 = True  # The second auroral peak exists
            fpt = (mlat1>ebmlat1) & (mlat1<pbmlat1)
            if np.sum(fpt)>=LLP2:
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
            pp.loc[(k0, k1),:] = [pbmlat1, mlatm1, ebmlat1, pbmlat2, mlatm2, ebmlat2,
                                  efm1, efm2, utt]
    # 7, remove defective points
    for k00, k0 in enumerate(('N', 'S')):
        # IF pbmlat1 or ebmlat2 is `mlatt` degrees greater or less than the
        # median value of neighboring `points` points, the corresponding
        # points are removed
        points = 21
        mlatt = 10
        pp0 = pp.loc[(k0, slice(None)), ['pbmlat1', 'ebmlat2']].values
        ppt = []  # initialize
        for k1 in range(points):
            ppt.append(np.roll(pp0, k1-int((points-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp1 = np.nanmedian(pp1, axis=2)
        fpt = np.sum(np.abs(pp1-pp0)>mlatt, axis=1)>=1
        pp.loc[fpt] = np.nan

    # 8, smooth result
    for k00, k0 in enumerate(('N', 'S')):
        # size of the smoothing window
        # increase this value if you don't care about small structure
        points = 11
        pp0 = pp.loc[(k0, slice(None))].values
        ppt = []  # initialize
        for k1 in range(points):
            ppt.append(np.roll(pp0, k1-int((points-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp0 = np.nanmedian(pp1, axis=2)
        fp = np.sum(np.isnan(pp1), axis=2)>=(points-1)/2+1
        pp0[fp] = np.nan
        pp.loc[(k0, slice(None))] = pp0

    # 9, If the total number of MLTs < some vlaue, we think there is no auroral
    for k00, k0 in enumerate(('N', 'S')):
        pp0 = pp.loc[(k0, slice(None)), 'pbmlat1'].values
        if np.sum(~np.isnan(pp0))<=LLN:
            pp.loc[(k0,slice(None))] = np.nan

    # 10, calculate datetime
    yyyy = int(ssusi['YEAR'][:])
    ddd = int(ssusi['DOY'][:])
    stime = pd.Timestamp(yyyy, 1, 1, 0, 0, 0)
    dt1 = pd.Timedelta(ddd-1, 'D')
    dt2 = pd.TimedeltaIndex(pp['datetime'].values, 'h')
    pp['datetime'] = stime+dt1+dt2
    return pp


def delete_find_parameters_one_file(fn):
    '''
    Get the peak value, peak latitude, and full width of the aurora
    input:
        fn        : ssusi aurora file name
    output:
        pp        : the parameters
    '''
    # Read data
    ssusi = Dataset(fn)

    #--------------------------------------------------
    #--------------adjustable parameters---------------
    #
    BinMLT = 0.25 # MLT length in one bin
    LLP1 = 200 # minimum data points in one bin
    LLP2 = 20 # minimum data points in auroral oval
    Scale1 = 0.5 # proportion of residual pixels that has energy flux>boundary values
    Scale2 = 50  # a maximum exists when it>Scale2*bvv
    BinW = 1*float(ssusi['PIXELSIZE_GEOMAGNETIC_LATITUDE'][:]) # size of the mlat bin window
    # Increase these values to get a narrower boundaries
    EdgeF = 0.3 # scale factor of max energy flux at boundaries
    # Increase these values to get a wider boundaries
    EdgeN = 3 # 'EdgeM' of 'EdgeN' consecutive data points are smaller than boundary values
    EdgeM = 3
    # minimum delta mlat between two auroral ovals
    LLdMLat = 5*float(ssusi['PIXELSIZE_GEOMAGNETIC_LATITUDE'][:])
    #--------------------------------------------------
    #--------------------------------------------------

    # initialize output
    pp = pd.DataFrame(np.ones([192, 9])*np.nan,
                      columns=['pbmlat1', 'mlatm1', 'ebmlat1',
                               'pbmlat2', 'mlatm2', 'ebmlat2',
                               'peak1', 'peak2', 'datetime'],
                      index=pd.MultiIndex.from_product(
                              [['N', 'S'], np.arange(BinMLT/2, 24, BinMLT)]))

    # Read variables that will be used
    orbitn = ssusi.STARTING_ORBIT_NUMBER # orbit number
    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE # value in a no data bin
    MLat = np.array(ssusi['LATITUDE_GEOMAGNETIC_GRID_MAP']) # magnetic latitude
    MLT = np.array(ssusi['MLT_GRID_MAP']) # magnetic local time
    energyfluxn = np.array(ssusi['ENERGY_FLUX_NORTH_MAP']) # electron energy flux
    energyfluxs = np.array(ssusi['ENERGY_FLUX_SOUTH_MAP'])
    utn = np.array(ssusi['UT_N'])
    uts = np.array(ssusi['UT_S'])
    bv = ssusi['ELECTRON_ENERGY_FLUX_THRESHOLDS'][:]

    # exclude pixels without data
    fpn = utn != nodatavalue
    MLatn, MLTn, energyfluxn, utn = (k[fpn] for k in (MLat, MLT, energyfluxn, utn))
    fps = uts != nodatavalue
    MLats, MLTs, energyfluxs, uts = (k[fps] for k in (MLat, MLT, energyfluxs, uts))

    # The UT time of the first file in one day may be for the previous day
    starttime = ssusi.STARTING_TIME   # format: yyyydddhhmmss
    stoptime = ssusi.STOPPING_TIME
    if starttime[:7] != stoptime[:7]:
        utn[utn>20] = utn[utn>20]-24  # Now, utn can be negative, negative valus mean yesterday
        uts[uts>20] = uts[uts>20]-24

    for k00, k0 in enumerate(['N', 'S']):
        # north or south
        if k0 == 'N':
            mlat0, mlt0, ef0, ut0 = (MLatn, MLTn, energyfluxn, utn)
            bvv = bv[0]
        if k0 == 'S':
            mlat0, mlt0, ef0, ut0 = (MLats, MLTs, energyfluxs, uts)
            bvv = bv[1]

        for k11, k1 in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
            peak1, peak2 = False, False

            # Select data of required MLT range
            lmlt, rmlt = k1-BinMLT/2, k1+BinMLT/2
            fp = (mlt0>=lmlt) & (mlt0<=rmlt)
            if np.sum(fp)<=LLP1:
                #print(k0, k1, np.sum(fp))
                continue
            mlat1, ef1, ut1 = mlat0[fp], ef0[fp], ut0[fp]

            # sort MLat
            idx = np.argsort(mlat1)
            mlat1, ef1, ut1 = mlat1[idx], ef1[idx], ut1[idx]

            # bin energy flux
            mlat11 = np.arange(50+BinW/2, 90, BinW)
            ef11 = np.arange(50+BinW/2, 90, BinW)*np.nan
            for k22, k2 in enumerate(mlat11):
                lmlat, rmlat = k2-BinW/2, k2+BinW/2
                fpt = (mlat1>lmlat) & (mlat1<rmlat)
                if np.sum(fpt)>0:
                    ef11[k22] = ef1[fpt].max()  # envelope of energy flux
            fpt = ~np.isnan(ef11)
            mlat11, ef11 = mlat11[fpt], ef11[fpt]

            # Find maximum energy flux
            fpmax1 = np.argmax(ef11)
            mlatm1, efm1 = mlat11[fpmax1], ef11[fpmax1]
            if efm1 < bvv*Scale2:
                continue

            # Find boundaries, boundaries meet the conditions:
            # 1, energy flux < some value for some consecutive values
            # or 2, difference of Mlat > some value
            efb1 = EdgeF*efm1  # boundary energy flux
            # poleward boundary
            efef = np.ones([len(ef11), EdgeN])*np.nan
            for k2 in range(EdgeN):
                efef[:, k2] = np.roll(ef11, -k2)
                efef[-(k2+1):, k2] = 0
            fpt = np.sum(efef<efb1, axis=1) >= EdgeM
            idx2 = np.where(fpt)[0]
            idxpb = idx2[idx2>=fpmax1].min()
            # equatorward boundary
            efef = np.ones([len(ef11), EdgeN])*np.nan
            for k2 in range(EdgeN):
                efef[:, k2] = np.roll(ef11, k2)
                efef[:k2+1, k2] = 0
            fpt = np.sum(efef<efb1, axis=1) >= EdgeM
            idx2 = np.where(fpt)[0]
            idxeb = idx2[idx2<=fpmax1].max()
            # equatorward and poleward boundary mlat
            ebmlat1, pbmlat1 = mlat11[idxeb], mlat11[idxpb]

            # determine the existence of the second auroral
            mlat22, ef22 = mlat11*1, ef11*1  # *1, so 22 and 11 points to different object
            fpt = (mlat11>=ebmlat1) & (mlat11<=pbmlat1)
            mlat22[fpt] = 0
            ef22[fpt] = 0
            if np.sum(ef22>efb1)>Scale1*np.sum(fpt):
                # Find the second maximum
                fpmax2 = np.argmax(ef22)
                mlatm2, efm2 = mlat22[fpmax2], ef22[fpmax2]
                efb2 = EdgeF*efm2

                # poleward boundary
                efef = np.ones([len(ef22), EdgeN])*np.nan
                for k2 in range(EdgeN):
                    efef[:, k2] = np.roll(ef22, -k2)
                    efef[-(k2+1):, k2] = 0
                fpt = np.sum(efef<efb2, axis=1) >= EdgeM
                idx2 = np.where(fpt)[0]
                idxpb = idx2[idx2>=fpmax2].min()
                # equatorward boundary
                efef = np.ones([len(ef22), EdgeN])*np.nan
                for k2 in range(EdgeN):
                    efef[:, k2] = np.roll(ef22, k2)
                    efef[:k2+1, k2] = 0
                fpt = np.sum(efef<efb2, axis=1) >= EdgeM
                idx2 = np.where(fpt)[0]
                idxeb = idx2[idx2<=fpmax2].max()
                ebmlat2, pbmlat2 = mlat22[idxeb], mlat22[idxpb]

                # if the two auroral ovals are close, they are combined
                if (ebmlat2-pbmlat1<=LLdMLat) and (ebmlat1-pbmlat2<=LLdMLat):
                    ebmlat1 = np.min([ebmlat1, ebmlat2])
                    pbmlat1 = np.max([pbmlat1, pbmlat2])
                else:
                    fpt = (mlat1>ebmlat2) & (mlat1<pbmlat2)
                    if np.sum(fpt)>LLP2:
                        peak2 = True  # The second auroral peak exists
            fpt = (mlat1>ebmlat1) & (mlat1<pbmlat1)
            if np.sum(fpt)>=LLP2:
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
            pp.loc[(k0, k1),:] = [pbmlat1, mlatm1, ebmlat1, pbmlat2, mlatm2, ebmlat2,
                                  efm1, efm2, utt]
    # remove defective points
    for k00, k0 in enumerate(('N', 'S')):
        # IF pbmlat1 or ebmlat2 is `mlatt` degrees greater or less than the
        # median value of neighboring `points` points, the corresponding
        # points are removed
        points = 21
        mlatt = 10
        pp0 = pp.loc[(k0, slice(None)), ['pbmlat1', 'ebmlat2']].values
        ppt = []  # initialize
        for k1 in range(points):
            ppt.append(np.roll(pp0, k1-int((points-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp1 = np.nanmedian(pp1, axis=2)
        fpt = np.sum(np.abs(pp1-pp0)>mlatt, axis=1)>=1
        pp.loc[fpt] = np.nan

    # smooth result
    for k00, k0 in enumerate(('N', 'S')):
        points = 5  # size of the smoothing window
        pp0 = pp.loc[(k0, slice(None))].values
        ppt = []  # initialize
        for k1 in range(points):
            ppt.append(np.roll(pp0, k1-int((points-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp0 = np.nanmedian(pp1, axis=2)
        fp = np.sum(np.isnan(pp1), axis=2)>=(points-1)/2+1
        pp0[fp] = np.nan
        pp.loc[(k0, slice(None))] = pp0
    yyyy = int(ssusi['YEAR'][:])
    ddd = int(ssusi['DOY'][:])
    stime = pd.Timestamp(yyyy, 1, 1, 0, 0, 0)
    dt1 = pd.Timedelta(ddd-1, 'D')
    dt2 = pd.TimedeltaIndex(pp['datetime'].values, 'h')
    pp['datetime'] = stime+dt1+dt2

    return pp


def test_find_parameters_one_file(ax, fn, ns='N'):
    '''
    Input:
        ax        : a polar axis handle
        fn        : ssusi file name
        ns        : 'N' or 'S'
    '''
    ax, hs = image_energy_flux_one_file(ax, fn, ns, vmin=0, vmax=10, s=1)
    pp = find_parameters_one_file(fn)
    eee = pp.loc[(ns, slice(None))][
            ['pbmlat1', 'mlatm1', 'ebmlat1', 'pbmlat2', 'mlatm2', 'ebmlat2']]
    r1 = 90-eee['pbmlat1']
    r2 = 90-eee['ebmlat2']
    r3 = 90-eee['ebmlat1']
    r4 = 90-eee['pbmlat2']
    theta = eee.index/12*np.pi
    ax.scatter(theta, r1, s=1, c='r', alpha=0.8)
    ax.scatter(theta, r2, s=1, c='r', alpha=0.8)
    ax.scatter(theta, r3, s=1, c='r', alpha=0.8)
    ax.scatter(theta, r4, s=1, c='r', alpha=0.8)
    return ax


def find_parameters_2013():
    import os
    import fnmatch
    import omni
    fpath = '/home/guod/big/raid4/lecai/sussi/'
    savepath = '/home/guod/big/raid4/guod/ssusi/'
    fna = os.listdir(fpath)
    # Read IMF By, Bz and AE in 2013
    print('Begin to read solar wind speed, IMF, AE and Dst')
    imfae = omni.get_omni(
            bdate='2013-1-1', edate='2014-1-1',
            variables=['Bym', 'Bzm', 'AE', 'V'], res='1m')
    imfae['Bt'] = np.sqrt(imfae['Bym']**2 + imfae['Bzm']**2)
    imfae['nwcf'] = imfae['V']**(4/3) * \
            imfae['Bt']**(2/3) * \
            ((1-imfae['Bzm']/imfae['Bt'])/2)**(4/3)
    dst = omni.get_omni(bdate='2013-1-1', edate='2014-1-1',
                        variables=['DST'], res='1h')
    dst = dst.reindex(imfae.index, method='nearest')
    imfae['Dst'] = dst.values
    print('End of reading solar wind speed, IMF, AE and Dst')
    for k0 in range(16,19): # satellite
        pp = []
        savefn = 'F{:d}_2013.dat'.format(k0)
        for k1 in range(1,13): # month
            print('Begin month {:02d}'.format(k1))
            for k2 in range(1,32): # day
                # 00 means full orbit?
                fn = fnmatch.filter(
                        fna, '*F{:d}*2013{:02d}{:02d}*00_DF.NC'.format(k0, k1, k2))
                if not fn:
                    continue
                for k33, k3 in enumerate(fn):
                    print('Processing '
                          'satellite {:s}, orbit {:s}'.format(k3[36:39], k3[-14:-9]))
                    pp0 = find_parameters_one_file(fpath+k3)
                    pp0 = pp0.loc[~(pp0['pbmlat1'].isnull().values)]
                    if not pp0.empty:
                        pp.append(pp0)

                        fig = plt.figure()
                        ax = plt.subplot(121, polar=True)
                        # North
                        image_energy_flux_one_file(ax, fpath+k3, 'N', s=1,
                                                   vmin=0, vmax=10)
                        ppt = pp0.loc[('N', slice(None))]
                        r1 = 90-ppt['pbmlat1']
                        r2 = 90-ppt['ebmlat2']
                        r3 = 90-ppt['ebmlat1']
                        r4 = 90-ppt['pbmlat2']
                        theta = ppt.index/12*np.pi
                        ax.scatter(theta, r1, s=1, c='r', alpha=0.8)
                        ax.scatter(theta, r2, s=1, c='r', alpha=0.8)
                        ax.scatter(theta, r3, s=1, c='r', alpha=0.8)
                        ax.scatter(theta, r4, s=1, c='r', alpha=0.8)
                        ax.set_title('North')
                        # South
                        ax = plt.subplot(122, polar=True)
                        image_energy_flux_one_file(ax, fpath+k3, 'S', s=1,
                                                   vmin=0, vmax=10)
                        ppt = pp0.loc[('S', slice(None))]
                        r1 = 90-ppt['pbmlat1']
                        r2 = 90-ppt['ebmlat2']
                        r3 = 90-ppt['ebmlat1']
                        r4 = 90-ppt['pbmlat2']
                        theta = ppt.index/12*np.pi
                        ax.scatter(theta, r1, s=1, c='r', alpha=0.8)
                        ax.scatter(theta, r2, s=1, c='r', alpha=0.8)
                        ax.scatter(theta, r3, s=1, c='r', alpha=0.8)
                        ax.scatter(theta, r4, s=1, c='r', alpha=0.8)
                        ax.set_title('South')
                        plt.savefig(savepath+'F{:d}_2013{:02d}{:02d}_{:s}.png'.format(
                                k0, k1, k2, k3[-14:-9]))
                        plt.close(fig)
                    else:
                        print('    Empty file')
            print('End of month {:02d}'.format(k1))
        pp = pd.concat(pp, axis=0)
        imfaet = imfae.reindex(pp['datetime'], method='nearest')
        pp['Bym'] = imfaet['Bym'].values
        pp['Bzm'] = imfaet['Bzm'].values
        pp['AE'] = imfaet['AE'].values
        pp['Dst'] = imfaet['Dst'].values
        pp['nwcf'] = imfaet['nwcf'].values
        pp.to_pickle(savepath+savefn)
    return


def parameters_imfae():
    import matplotlib.pyplot as plt
    xvar = 'AE'
    ns = 'S'
    iMLT = 18
    sat = 'F17'

    MLT = np.arange(0.25/2, 24, 0.25)
    print(MLT[iMLT])
    fpath = '/home/guod/big/raid4/guod/ssusi/'
    ax = plt.subplot()
    for k0 in range(1, 13):
        pp = pd.read_pickle(fpath+'{:s}_2013{:02d}.dat'.format(sat, k0))
        ppt = pp.loc[(ns, MLT[iMLT]), :]
        #width = ppt['peak1']
        fp = ppt['pbmlat1']-ppt['ebmlat2'] <4
        x = ppt.loc[fp, xvar]
        y = ppt.loc[fp, 'peak1']
        ax.scatter(x, y,s=1,c='b')
    plt.ylim(0,50)
    plt.show()
    return

if __name__=='__main__':
    #  plt.close('all')
    #  ax = plt.subplot(polar=True)
    #  a = test_find_parameters_one_file(ax, '/home/guod/big/raid4/lecai/sussi/'
    #          'PS.APL_V0105S024CE0018_SC.U_DI.A_GP.F16-SSUSI_PA.APL-EDR-'
    #          'AURORA_DD.20130204_SN.47986-00_DF.NC', ns='S')
    #  plt.show()
    find_parameters_2013()

    print('--------------------------------------------------------------------------------')
    print('Remaining problems:')
    print('  1, The smoothing process can not handle all the conditions')
    print('  2, Sometime the swept regions of the satellite cut auroral region'
          ' in the midnight side')
    print('--------------------------------------------------------------------------------')
