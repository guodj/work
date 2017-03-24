# Auroral precipation model using SSUSI data
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import myfunctions as mf
from scipy.optimize import curve_fit

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
    cb = plt.colorbar(hs, ax=ax, pad=0.1)
    cb.set_label(ssusi[var].UNITS)
    return ax, hs


def find_parameters_one_file(fn):
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
    binmlt = 0.25 # MLT length in one bin
    ldp = 200 # minimum data points in one bin
    ldp2 = 20 # minimum data points in auroral oval
    scale2 = 0.5 # proportion of residual pixels that has energy flux>boundary values
    scale3 = 50  # a maximum exists when it>scale3*bvv
    binw = 1*float(ssusi['PIXELSIZE_GEOMAGNETIC_LATITUDE'][:]) # size of the mlat bin window
    # Increase these values to get a narrower boundaries
    edgef = 0.3 # scale factor of max energy flux at boundaries
    # Increase these values to get a wider boundaries
    edgn = 3 # 'edgnm' of 'edgn' consecutive data points are smaller than boundary values
    edgnm = 3
    # minimum delta mlat between two auroral ovals
    edgdmlat2 = 5*float(ssusi['PIXELSIZE_GEOMAGNETIC_LATITUDE'][:])
    #--------------------------------------------------
    #--------------------------------------------------

    # initialize output
    pp = pd.DataFrame(np.ones([192, 9])*np.nan,
                      columns=['pbmlat1', 'mlatm1', 'ebmlat1',
                               'pbmlat2', 'mlatm2', 'ebmlat2',
                               'peak1', 'peak2', 'datetime'],
                      index=pd.MultiIndex.from_product(
                              [['N', 'S'], np.arange(binmlt/2, 24, binmlt)]))

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

        for k11, k1 in enumerate(np.arange(binmlt/2, 24, binmlt)):
            peak1, peak2 = False, False

            # Select data of required MLT range
            lmlt, rmlt = k1-binmlt/2, k1+binmlt/2
            fp = (mlt0>=lmlt) & (mlt0<=rmlt)
            if np.sum(fp)<=ldp:
                #print(k0, k1, np.sum(fp))
                continue
            mlat1, ef1, ut1 = mlat0[fp], ef0[fp], ut0[fp]

            # sort MLat
            idx = np.argsort(mlat1)
            mlat1, ef1, ut1 = mlat1[idx], ef1[idx], ut1[idx]

            # bin energy flux
            mlat11 = np.arange(50+binw/2, 90, binw)
            ef11 = np.arange(50+binw/2, 90, binw)*np.nan
            for k22, k2 in enumerate(mlat11):
                lmlat, rmlat = k2-binw/2, k2+binw/2
                fpt = (mlat1>lmlat) & (mlat1<rmlat)
                if np.sum(fpt)>0:
                    ef11[k22] = ef1[fpt].max()  # envelope of energy flux
            fpt = ~np.isnan(ef11)
            mlat11, ef11 = mlat11[fpt], ef11[fpt]

            # Find maximum energy flux
            fpmax1 = np.argmax(ef11)
            mlatm1, efm1 = mlat11[fpmax1], ef11[fpmax1]
            if efm1 < bvv*scale3:
                continue

            # Find boundaries, boundaries meet the conditions:
            # 1, energy flux < some value for some consecutive values
            # or 2, difference of Mlat > some value
            efb1 = edgef*efm1  # boundary energy flux
            # poleward boundary
            efef = np.ones([len(ef11), edgn])*np.nan
            for k2 in range(edgn):
                efef[:, k2] = np.roll(ef11, -k2)
                efef[-(k2+1):, k2] = 0
            fpt = np.sum(efef<efb1, axis=1) >= edgnm
            idx2 = np.where(fpt)[0]
            idxpb = idx2[idx2>=fpmax1].min()
            # equatorward boundary
            efef = np.ones([len(ef11), edgn])*np.nan
            for k2 in range(edgn):
                efef[:, k2] = np.roll(ef11, k2)
                efef[:k2+1, k2] = 0
            fpt = np.sum(efef<efb1, axis=1) >= edgnm
            idx2 = np.where(fpt)[0]
            idxeb = idx2[idx2<=fpmax1].max()
            # equatorward and poleward boundary mlat
            ebmlat1, pbmlat1 = mlat11[idxeb], mlat11[idxpb]

            # determine the existence of the second auroral
            mlat22, ef22 = mlat11*1, ef11*1  # *1, so 22 and 11 points to different object
            fpt = (mlat11>=ebmlat1) & (mlat11<=pbmlat1)
            mlat22[fpt] = 0
            ef22[fpt] = 0
            if np.sum(ef22>efb1)>scale2*np.sum(fpt):
                # Find the second maximum
                fpmax2 = np.argmax(ef22)
                mlatm2, efm2 = mlat22[fpmax2], ef22[fpmax2]
                efb2 = edgef*efm2

                # poleward boundary
                efef = np.ones([len(ef22), edgn])*np.nan
                for k2 in range(edgn):
                    efef[:, k2] = np.roll(ef22, -k2)
                    efef[-(k2+1):, k2] = 0
                fpt = np.sum(efef<efb2, axis=1) >= edgnm
                idx2 = np.where(fpt)[0]
                idxpb = idx2[idx2>=fpmax2].min()
                # equatorward boundary
                efef = np.ones([len(ef22), edgn])*np.nan
                for k2 in range(edgn):
                    efef[:, k2] = np.roll(ef22, k2)
                    efef[:k2+1, k2] = 0
                fpt = np.sum(efef<efb2, axis=1) >= edgnm
                idx2 = np.where(fpt)[0]
                idxeb = idx2[idx2<=fpmax2].max()
                ebmlat2, pbmlat2 = mlat22[idxeb], mlat22[idxpb]

                # if the two auroral ovals are close, they are combined
                if (ebmlat2-pbmlat1<=edgdmlat2) and (ebmlat1-pbmlat2<=edgdmlat2):
                    ebmlat1 = np.min([ebmlat1, ebmlat2])
                    pbmlat1 = np.max([pbmlat1, pbmlat2])
                else:
                    fpt = (mlat1>ebmlat2) & (mlat1<pbmlat2)
                    if np.sum(fpt)>ldp2:
                        peak2 = True  # The second auroral peak exists
            fpt = (mlat1>ebmlat1) & (mlat1<pbmlat1)
            if np.sum(fpt)>=ldp2:
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
    ax, hs = image_energy_flux_one_file(ax, fn, ns, vmin=0, vmax=10, s=1)
    pp = find_parameters_one_file(fn)
    eee = pp.loc[(ns, slice(None))][
            ['pbmlat1', 'mlatm1', 'ebmlat1', 'pbmlat2', 'mlatm2', 'ebmlat2']]
    r1 = 90-eee['pbmlat1']
    r2 = 90-eee['ebmlat2']
    r3 = 90-eee['ebmlat1']
    r4 = 90-eee['pbmlat2']
    theta = eee.index/12*np.pi
    ax.scatter(theta, r1, s=1, c='r')
    ax.scatter(theta, r2, s=1, c='r')
    ax.scatter(theta, r3, s=1, c='r')
    ax.scatter(theta, r4, s=1, c='r')
    return pp


if __name__=='__main__':
    import matplotlib.pyplot as plt
    plt.close('all')
    fig = plt.figure(figsize=(10.54, 3.8))
    fn = ('/home/guod/big/raid4/lecai/sussi/PS.APL_V0105S024CB0005_SC.U_DI.A_GP.F18'
          '-SSUSI_PA.APL-EDR-AURORA_DD.20130528_SN.18598-00_DF.NC')
    #fn = ('/home/guod/big/raid4/lecai/sussi/PS.APL_V0105S024CE0018_SC.U_DI.A_GP.F16'
    #      '-SSUSI_PA.APL-EDR-AURORA_DD.20130528_SN.49573-00_DF.NC')
    #fn = ('/home/guod/big/raid4/lecai/sussi/PS.APL_V0105S024CB0005_SC.U_DI.A_GP.F18'
    #      '-SSUSI_PA.APL-EDR-AURORA_DD.20130525_SN.18556-00_DF.NC')
    ax=plt.subplot(1,2,1,polar=True)
    a = test_find_parameters_one_file(ax, fn, ns='N')

    ax=plt.subplot(1,2,2,polar=True)
    a = test_find_parameters_one_file(ax, fn, ns='S')
    plt.show()
    print('--------------------------------------------------------------------------------')
    print('Remaining problems:')
    print('  1, The smoothing process can not handle all the conditions')
    print('  2, Sometime the swept regions of the satellite cut auroral region'
          ' in the midnight side')
    print('--------------------------------------------------------------------------------')
