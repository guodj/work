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
    bv = np.array(ssusi['ELECTRON_ENERGY_FLUX_THRESHOLDS']) # energy flux threshold for boundary
    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE # value in a no data bin
    bv1 = bv[0] if ns=='N' else bv[1]
    if ns=='N':
        var='ENERGY_FLUX_NORTH_MAP'
    if ns=='S':
        var='ENERGY_FLUX_SOUTH_MAP'
    variable = np.array(ssusi[var])
    # Only include data points inside the boundary
    fp = (variable != nodatavalue)
    MLat, MLT, variable = (k[fp] for k in (MLat, MLT, variable))
    # plot variable
    r = 90-MLat
    theta = MLT/12*np.pi
    hs = ax.scatter(theta, r, s=s, c=variable, vmin=vmin,
                    vmax=vmax, alpha=alpha)
    # Set polar coordinates
    ax.set_rgrids(np.arange(10, 31, 10), (80, 70, 60))
    ax.set_thetagrids(np.arange(0, 361, 90), (0, 6, 12, 18))
    ax.set_theta_zero_location('S')
    ax.set_rlim(0, 30)
    ax.set_title(ssusi[var].TITLE)
    ax.text(0.8, -0.1,
            'YEAR: '+str(np.int(ssusi['YEAR'][:])) +
            '\nDOY: '+str(np.int(ssusi['DOY'][:])) +
            '\nORBIT: '+str(np.int(ssusi.STARTING_ORBIT_NUMBER)),
            transform=ax.transAxes,
            fontsize=8)
    cb = plt.colorbar(hs, ax=ax, orientation='horizontal')
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

    # adjustable parameters
    binmlt = 0.25 # MLT length in one bin
    ldp = 80 # minimum data points in one bin
    ldp2 = 20 # minimum data points in auroral oval
    scale2 = 0.8 # proportion of residual pixels that has energy flux>boundary values
    binw = 1*float(ssusi['PIXELSIZE_GEOMAGNETIC_LATITUDE'][:]) # size of the mlat bin window
    # Increase these values to get a narrower boundaries
    edgef = 0.3 # scale factor of max energy flux at boundaries
    # Increase these values to get a wider boundaries
    # mlat difference that can be seen as the boundary
    edgdmlat = 3*float(ssusi['PIXELSIZE_GEOMAGNETIC_LATITUDE'][:])
    edgn = 3 # 'edgnm' of 'edgn' consecutive data points are smaller than boundary values
    edgnm = 3
    # minimum delta mlat between two auroral ovals
    edgdmlat2 = 5*float(ssusi['PIXELSIZE_GEOMAGNETIC_LATITUDE'][:])

    # initialize output
    pp = pd.DataFrame(np.ones([192, 8])*np.nan,
                      columns=['pbmlat1', 'mlatm1', 'ebmlat1',
                               'pbmlat2', 'mlatm2', 'ebmlat2',
                               'peak1', 'peak2'],
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

    # exclude pixels without data
    fpn = energyfluxn != nodatavalue
    MLatn, MLTn, energyfluxn, utn = (k[fpn] for k in (MLat, MLT, energyfluxn, utn))
    fps = energyfluxs != nodatavalue
    MLats, MLTs, energyfluxs, uts = (k[fps] for k in (MLat, MLT, energyfluxs, uts))

    for k00, k0 in enumerate(['N', 'S']):
        # north or south
        if k0 == 'N':
            mlat0, mlt0, ef0, ut0 = (MLatn, MLTn, energyfluxn, utn)
        if k0 == 'S':
            mlat0, mlt0, ef0, ut0 = (MLats, MLTs, energyfluxs, uts)

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
                    ef11[k22] = ef1[fpt].max()  # envelope
            fpt = ~np.isnan(ef11)
            mlat11, ef11 = mlat11[fpt], ef11[fpt]

            # Find maximum energy flux
            fpmax1 = np.argmax(ef11)
            mlatm1, efm1 = mlat11[fpmax1], ef11[fpmax1]

            # Find boundaries, boundaries meet the conditions:
            # 1, energy flux < some value for some consecutive values
            # or 2, difference of Mlat > some value
            efb1 = edgef*efm1  # boundary energy flux
            # poleward boundary
            efef = np.ones([len(ef11), edgn])*np.nan
            for k2 in range(edgn):
                efef[:, k2] = np.roll(ef11, -k2)
                if k2 != 0:
                    efef[-k2:] = 0
            fpt0 = np.sum(efef<efb1, axis=1) >= edgnm
            fpt1 = np.append(np.diff(mlat11), 1e10)>= edgdmlat
            fpt = (fpt0) | (fpt1)
            idx2 = np.where(fpt)[0]
            idxpb = idx2[idx2>=fpmax1].min()
            # equatorward boundary
            efef = np.ones([len(ef11), edgn])*np.nan
            for k2 in range(edgn):
                efef[:, k2] = np.roll(ef11, k2)
                if k2 != 0:
                    efef[:k2] = 0
            fpt0 = np.sum(efef<efb1, axis=1) >= edgnm
            fpt1 = np.insert(np.diff(mlat11), 0, 1e10)>= edgdmlat
            fpt = (fpt0) | (fpt1)
            idx2 = np.where(fpt)[0]
            idxeb = idx2[idx2<=fpmax1].max()
            # equatorward and poleward boundary mlat
            ebmlat1, pbmlat1 = mlat11[idxeb], mlat11[idxpb]

            # determine the existence of the second auroral
            fpt = (mlat11<=ebmlat1) | (mlat11>=pbmlat1)
            mlat22, ef22 = mlat11[fpt], ef11[fpt]
            if np.sum(ef22>efb1)>scale2*np.sum((mlat11>=ebmlat1) & (mlat11<=pbmlat1)):
                # Find the second maximum
                fpmax2 = np.argmax(ef22)
                mlatm2, efm2 = mlat22[fpmax2], ef22[fpmax2]
                efb2 = edgef*efm2

                # poleward boundary
                efef = np.ones([len(ef22), edgn])*np.nan
                for k2 in range(edgn):
                    efef[:, k2] = np.roll(ef22, -k2)
                    if k2 != 0:
                        efef[-k2:] = 0
                fpt0 = np.sum(efef<efb2, axis=1) >= edgnm
                fpt1 = np.append(np.diff(mlat22), 1e10)>edgdmlat
                fpt = (fpt0) | (fpt1)
                idx2 = np.where(fpt)[0]
                idxpb = idx2[idx2>=fpmax2].min()
                # equatorward boundary
                efef = np.ones([len(ef22), edgn])*np.nan
                for k2 in range(edgn):
                    efef[:, k2] = np.roll(ef22, k2)
                    if k2 != 0:
                        efef[:k2] = 0
                fpt0 = np.sum(efef<efb2, axis=1) >= edgnm
                fpt1 = np.insert(np.diff(mlat22), 0, 1e10)>edgdmlat
                fpt = (fpt0) | (fpt1)
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
            pp.loc[(k0, k1),:] = [pbmlat1, mlatm1, ebmlat1, pbmlat2, mlatm2, ebmlat2,
                                  efm1, efm2]
    # smooth result
    for k00, k0 in enumerate(('N', 'S')):
        points = 5  # size of the smoothing window
        pp0 = pp.loc[(k0, slice(None))].values
        ppt = []
        for k1 in range(points):
            ppt.append(np.roll(pp0, k1-int((points-1)/2), axis=0))
        pp1 = np.stack(ppt, axis=2)
        pp0 = np.nanmedian(pp1, axis=2)
        fp = np.sum(np.isnan(pp1), axis=2)>=3
        pp0[fp] = np.nan
        pp.loc[(k0, slice(None))] = pp0
    return pp


def save_parameters_to_file(fns, savefn):
    # fns is the file name list
    f = open(savefn, 'w')
    for k00, k0 in enumerate(fns):
        pp = find_parameters_one_file(k0)
def test_find_parameters_one_file(fn, ns='N'):
    plt.close('all')
    plt.figure(figsize=(5.2, 7.8))
    ax = plt.subplot(polar=True)
    ax, hs = image_energy_flux_one_file(ax, fn, ns, vmin=0, vmax=10, s=1)
    pp = find_parameters_one_file(fn)
    eee = pp.loc[(ns, slice(None))][
            ['pbmlat1', 'mlatm1', 'ebmlat1', 'pbmlat2', 'mlatm2', 'ebmlat2']]
    r1 = 90-eee['pbmlat1']
    r2 = 90-eee['ebmlat2']
    r3 = 90-eee['ebmlat1']
    r4 = 90-eee['pbmlat2']
    theta = eee.index/12*np.pi
    ax.scatter(theta, r1, s=3, c='r')
    ax.scatter(theta, r2, s=3, c='r')
    ax.scatter(theta, r3, s=3, c='r')
    ax.scatter(theta, r4, s=3, c='r')
    ax.set_rlim(0, 35)

    #    # adjustable parameters
    #    binmlt = 0.25   # MLT length in one bin
    #    ldp = 80 # least data points in one bin
    #    mlatnum = 20  # number between eemlat and pemlat
    #    # Read data
    #    ssusi = Dataset(fn)
    #    # Read variables that will be used
    #    orbitn = ssusi.STARTING_ORBIT_NUMBER # orbit number
    #    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE # value in a no data bin
    #    MLat = np.array(ssusi['LATITUDE_GEOMAGNETIC_GRID_MAP']) # magnetic latitude
    #    MLT = np.array(ssusi['MLT_GRID_MAP']) # magnetic local time
    #    energyfluxn = np.array(ssusi['ENERGY_FLUX_NORTH_MAP']) # electron energy flux
    #    energyfluxs = np.array(ssusi['ENERGY_FLUX_SOUTH_MAP'])
    #    # exclude points without data
    #    fpn = energyfluxn != nodatavalue
    #    MLatn, MLTn, energyfluxn = (k[fpn] for k in (MLat, MLT, energyfluxn))
    #    fps = energyfluxs != nodatavalue
    #    MLats, MLTs, energyfluxs = (k[fps] for k in (MLat, MLT, energyfluxs))
    #    if ns=='N':
    #        mlat0, mlt0, ef0 = (MLatn, MLTn, energyfluxn)
    #    if ns=='S':
    #        mlat0, mlt0, ef0 = (MLats, MLTs, energyfluxs)
    #    for k11, k1 in enumerate(np.arange(binmlt/2, 24, binmlt)):
    #        lmlt, rmlt = k1-binmlt/2, k1+binmlt/2
    #        fp = (mlt0>=lmlt) & (mlt0<=rmlt)
    #        if np.sum(fp)<=ldp:
    #            #print(k0, k1, np.sum(fp))
    #            continue
    #        mlat1, mlt1, ef1 = (k[fp] for k in (mlat0, mlt0, ef0))
    #        plt.figure()
    #        plt.scatter(mlat1, ef1, s=2, c='k')
    #        pp1 =pp[(pp.NS==ns) & (pp.MLT==k1)][['EeMLat1', 'PeMLat1', 'EF1',
    #                                             'EeMLat2', 'PeMLat2', 'EF2']]
    #        if pp1.empty:
    #            continue
    #        mlatbin = np.arange(pp1.EeMLat1+(pp1.PeMLat1-pp1.EeMLat1)/mlatnum/2,
    #                            pp1.PeMLat1, (pp1.PeMLat1-pp1.EeMLat1)/mlatnum)
    #        plt.plot(mlatbin, pp1.EF1.values[0], 'ro--')
    #        if not np.isnan(pp1.EF2.values[0]).all():
    #            mlatbin = np.arange(pp1.EeMLat2+(pp1.PeMLat2-pp1.EeMLat2)/mlatnum/2,
    #                                pp1.PeMLat2, (pp1.PeMLat2-pp1.EeMLat2)/mlatnum)
    #            plt.plot(mlatbin, pp1.EF2.values[0], 'bo--')
    #        plt.xlim(60, 90)
    #        plt.ylim(0, 40)
    return pp

def delete_test_gaussian_parameter_one_file(fn, ns='N'):
    import matplotlib.pyplot as plt
    rr = 2.35482 # ratio of parameter FWHM and standard deviation (c)
    gaus = lambda x, a1, b1, c1: a1*np.exp(-(x-b1)**2/(2*c1**2))
    binmlt =1 # MLT range in each bin
    fig, ax = plt.subplots(6, 4, sharex=True, sharey=True,figsize=[6.8, 8.8])
    # Read data file and variables
    ssusi = Dataset(fn)
    orbitn = ssusi.STARTING_ORBIT_NUMBER
    nodatavalue = ssusi.NO_DATA_IN_BIN_VALUE
    bv = np.array(ssusi['ELECTRON_ENERGY_FLUX_THRESHOLDS']) # energy flux threshold
    MLat = np.array(ssusi['LATITUDE_GEOMAGNETIC_GRID_MAP'])
    MLT = np.array(ssusi['MLT_GRID_MAP'])
    if False: # Method 1: use ENERGY_FLUX_NORTH_MAP and ENERGY_FLUX_SOUTH_MAP
        energyfluxn = np.array(ssusi['ENERGY_FLUX_NORTH_MAP'])
        energyfluxs = np.array(ssusi['ENERGY_FLUX_SOUTH_MAP'])
        # include data points greater than threshold
        fpn = energyfluxn >= bv[0]
        MLatn, MLTn, energyfluxn = (k[fpn] for k in (MLat, MLT, energyfluxn))
        fps = energyfluxs >= bv[1]
        MLats, MLTs, energyfluxs = (k[fps] for k in (MLat, MLT, energyfluxs))
    else:
        energyfluxn = np.array(ssusi['ELECTRON_FLUX_NORTH_BOUNDARY_MAP'])
        energyfluxs = np.array(ssusi['ELECTRON_FLUX_SOUTH_BOUNDARY_MAP'])
        # include data points > threshold and != 1e10
        fpn = (energyfluxn >= bv[0]) & (energyfluxn != 1e10)
        MLatn, MLTn, energyfluxn = (k[fpn] for k in (MLat, MLT, energyfluxn))
        fps = (energyfluxs >= bv[1]) & (energyfluxs != 1e10)
        MLats, MLTs, energyfluxs = (k[fps] for k in (MLat, MLT, energyfluxs))
    # North or South
    mlat0, mlt0, ef0 = (MLatn, MLTn, energyfluxn)
    if ns == 'S':
        mlat0, mlt0, ef0 = (MLats, MLTs, energyfluxs)
    p = gaussian_fit_one_file(fn)
    for k00, k0 in enumerate(np.arange(binmlt/2, 24, binmlt)):
        # original data
        lmlt, rmlt = k0-binmlt/2, k0+binmlt/2
        fp = (mlt0>=lmlt) & (mlt0<=rmlt)
        mlat1, mlt1, ef1 = (k[fp] for k in (mlat0, mlt0, ef0))
        plt.sca(ax[k00//4, k00%4])
        plt.scatter(mlat1, ef1, s=1)
        plt.text(0.1, 0.7,'MLT: %.1f'%k0, fontsize=10, transform=plt.gca().transAxes)
        if k00//4==5:
            plt.xlabel('MLat')
        if k00%4==0:
            plt.ylabel('Energy flux')
        # fitting values
        p0 = np.array(
                p[(p.MLT==k0) & (p.NS==ns)][['Peak', 'Center_MLat', 'FWHM']]).reshape(-1)
        if np.isnan(p0[0]):
            continue
        lat0 = np.arange(60, 90, 0.1)
        fitef = gaus(lat0, p0[0], p0[1], p0[2]/rr)
        plt.plot(lat0, fitef, 'r')
    plt.ylim(0, 40)
    plt.yticks(np.arange(0,41,10))
    plt.xlim(60, 80)
    plt.xticks(np.arange(60, 81, 5))
    plt.tight_layout(h_pad=0.05, w_pad=0.05)
    plt.show()
    return p



if __name__=='__main__':
    import matplotlib.pyplot as plt
    plt.close('all')
    fig = plt.figure(figsize=(10.54, 10))
    fn = ('/home/guod/big/raid4/lecai/sussi/PS.APL_V0105S024CB0005_SC.U_DI.A_GP.F18'
          '-SSUSI_PA.APL-EDR-AURORA_DD.20130526_SN.18570-02_DF.NC')
    #fn = ('/home/guod/big/raid4/lecai/sussi/PS.APL_V0105S024CB0005_SC.U_DI.A_GP.F18'
    #      '-SSUSI_PA.APL-EDR-AURORA_DD.20130525_SN.18556-00_DF.NC')
    s = Dataset(fn)
    utn = s['UT_N'][:]
    uts = s['UT_S'][:]
    mlat = s['LATITUDE_GEOMAGNETIC_GRID_MAP'][:]
    mlt = s['MLT_GRID_MAP'][:]
    nodatavalue = s.NO_DATA_IN_BIN_VALUE

    fpn = ~(utn==nodatavalue)
    utn, mlatn, mltn = utn[fpn], mlat[fpn], mlt[fpn]
    fps = ~(uts==nodatavalue)
    uts, mlats, mlts = uts[fps], mlat[fps], mlt[fps]

    rn = 90-mlatn
    thetan = mltn/12*np.pi
    ax=plt.subplot(2,2,1,polar=True)
    plt.scatter(thetan, rn, s=1, c=utn)
    ax.set_theta_zero_location('S')
    plt.colorbar()
    ax.text(0, 1.1, 'Max: {:.3f}'.format(utn.max()), transform=ax.transAxes)
    ax.text(0, 1.05, 'Min: {:.3f}'.format(utn.min()), transform=ax.transAxes)

    rs = 90-mlats
    thetas = mlts/12*np.pi
    ax=plt.subplot(2,2,2,polar=True)
    plt.scatter(thetas, rs, s=1, c=uts)
    ax.set_theta_zero_location('S')
    plt.colorbar()
    ax.text(0, 1.1, 'Max: {:.3f}'.format(uts.max()), transform=ax.transAxes)
    ax.text(0, 1.05, 'Min: {:.3f}'.format(uts.min()), transform=ax.transAxes)

    test_find_parameters_one_file
    plt.show()
    print('--------------------------------------------------------------------------------')
    print('Remaining problems:')
    print('  1, The smoothing process can not handle all the conditions')
    print('  2, Sometime the swept regions of the satellite cut auroral region'
          ' in the midnight side')
    print('--------------------------------------------------------------------------------')
