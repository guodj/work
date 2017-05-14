# Auroral precipation model using SSUSI data
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import myfunctions as mf
import scipy.signal as signal
import seaborn as sns
import glob
sns.set('paper', 'whitegrid')

# global variables
LLP1 = 200  # Lower limit of pixels in one bin
LLEFTOTH = 200  # lower limit of total energy flux in one hemisphere
BinMLT = 0.25  # MLT range in one bin
LLEFTOT = 10  # lower limit of total energy flux in each bin, for method 1,2
datapath = '/home/guod/WD4T/ssusi/'
savepath = '/home/guod/WD4T/work4/'
dratio = 0.05 # flag 5%, 10%, ..., 95% of the total energy flux


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


def find_parameters_one_file_5(fn):
    '''
    Set flags at the MLats that have 5%, 10%, ..., 95% of the
    total energy flux save MLats and energy flux
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
    meanenergyn = np.array(ssusi['ELECTRON_MEAN_NORTH_ENERGY_MAP'])
    meanenergys = np.array(ssusi['ELECTRON_MEAN_SOUTH_ENERGY_MAP'])
    utn = np.array(ssusi['UT_N'])
    uts = np.array(ssusi['UT_S'])
    # bv = np.array(ssusi['ELECTRON_ENERGY_FLUX_THRESHOLDS']) # boundary values

    # Apply median filter to the energy flux to remove impulse noise
    kernel_size_own = 3
    energyfluxn = signal.medfilt2d(energyfluxn, kernel_size_own)
    energyfluxs = signal.medfilt2d(energyfluxs, kernel_size_own)
    meanenergyn = signal.medfilt2d(meanenergyn, kernel_size_own)
    meanenergys = signal.medfilt2d(meanenergys, kernel_size_own)

    # exclude pixels without data
    fpn = (utn != nodatavalue)
    MLatn, MLTn, energyfluxn, utn, meanenergyn = (
            k[fpn] for k in (MLat, MLT, energyfluxn, utn, meanenergyn))
    fps = uts != nodatavalue
    MLats, MLTs, energyfluxs, uts, meanenergys = (
            k[fps] for k in (MLat, MLT, energyfluxs, uts, meanenergys))

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
    pp1 = np.ones([np.int(2*24/BinMLT), int(1/dratio-1)])*np.nan  # MLats
    pp2 = pd.DataFrame(
            index=np.arange(int(2*24/BinMLT)),
            columns=('ns', 'mlt', 'datetime', 'eftot'))
    pp3 = np.ones([np.int(2*24/BinMLT), int(1/dratio-1)])*np.nan  # ef
    pp4 = np.ones([np.int(2*24/BinMLT), int(1/dratio-1)])*np.nan  # mean energy

    # group data according to hemispheres
    for k00, k0 in enumerate(['N', 'S']):
        # select data in the specified hemisphere
        if k0 == 'N':
            mlat0, mlt0, ef0, ut0, me0 = (
                    MLatn, MLTn, energyfluxn, utn, meanenergyn)
            # bvv = bv[0]
        if k0 == 'S':
            mlat0, mlt0, ef0, ut0, me0 = (
                    MLats, MLTs, energyfluxs, uts, meanenergys)
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
            # If the total number of pixels in one MLT bin is small, exclude
            if np.sum(fp) <= LLP1:
                continue
            mlat1, ef1, ut1, me1 = mlat0[fp], ef0[fp], ut0[fp], me0[fp]
            # If the Mlat range is too small in one MLT bin, exclude
            if mlat1.min() > llmlat or mlat1.max() < ulmlat:
                continue

            # sort MLat
            idx = np.argsort(mlat1)
            mlat1, ef1, ut1, me1 = mlat1[idx], ef1[idx], ut1[idx], me1[idx]

            # calculate the integral energy flux
            efcum = ef1.cumsum()
            eftot = efcum[-1]
            # If the total energy flux in one MLT bin is small, exclude
            if eftot < LLEFTOT:
                continue

            # Mlats of 5%, 10%... energy flux -->pp1
            # energy fluxes at 5%, 10%, ... --> pp3
            # mean energy --> pp4
            pp1t = np.ones(int(1/dratio-1))*np.nan
            pp3t = np.ones(int(1/dratio-1))*np.nan
            pp4t = np.ones(int(1/dratio-1))*np.nan
            utt = np.ones(int(1/dratio-1))*np.nan
            for k22, k2 in enumerate(np.arange(dratio, 1, dratio)):
                fp = (efcum > eftot*k2)
                pp1t[k22] = mlat1[fp][0]
                pp3t[k22] = ef1[fp][0]
                pp4t[k22] = me1[fp][0]
                utt[k22] = ut1[fp][0]
            pp1[idpp, :] = pp1t
            pp3[idpp, :] = pp3t
            pp4[idpp, :] = pp4t

            # ns, MLT, UT, eftot --> pp2
            pp2.iloc[idpp, :] = [
                    k0, k1, date+pd.Timedelta(np.median(utt), 'h'), eftot]

    # exclude NaN
    fp = (np.isnan(pp1[:, 0]))
    pp1 = pp1[~fp]
    pp3 = pp3[~fp]
    pp4 = pp4[~fp]
    pp2 = pp2.loc[~fp, :]
    pp2 = pp2.astype({'ns': str, 'mlt': float,
                      'datetime': np.datetime64, 'eftot': float})
    return pp1, pp2, pp3, pp4  # Mlats, ..., EFs, MEs


def test_find_parameters_one_file_5(fn):
    '''
    Input:
        fn        : ssusi file name
    '''
    pp1, pp2, pp3, pp4 = find_parameters_one_file_5(fn)
    plt.figure(figsize=(8.7, 4.4))
    for k00, k0 in enumerate(['N', 'S']):
        ax = plt.subplot(1, 2, k00+1, polar=True)
        ax, hs = image_energy_flux_one_file(
                ax, fn, k0, vmin=0, vmax=10, s=1, alpha=1)
        fp = (pp2.ns == k0).values
        for k11, k1 in enumerate(range(int(1/dratio-1))):
            r = 90-pp1[fp, k11]
            theta = pp2.loc[fp, 'mlt']/12*np.pi
            ax.scatter(theta, r, c='b', s=1, alpha=0.5)
        plt.title(k0+'H')
    plt.tight_layout()
    return ax, hs


def find_parameters():
    import os
    import fnmatch
    import omni
    # Read IMF By, Bz and AE in 2013
    print('Begin to read solar wind speed, IMF, AE and Dst')
    imfae = omni.get_omni(
            bdate='2011-1-1', edate='2015-1-1',
            variables=['Bym', 'Bzm', 'AE', 'V'], res='1m')
    imfae['Bt'] = np.sqrt(imfae['Bym']**2 + imfae['Bzm']**2)
    imfae['nwcf'] = imfae['V']**(4/3) * imfae['Bt']**(2/3) * \
        ((1-imfae['Bzm']/imfae['Bt'])/2)**(4/3)/100.0
    imfae['EKL'] = 0.5*imfae['V']*(imfae['Bt']-imfae['Bzm'])/1e6
    dst = omni.get_omni(bdate='2011-1-1', edate='2015-1-1',
                        variables=['DST'], res='1h')
    dst = dst.reindex(imfae.index, method='nearest')
    imfae['Dst'] = dst.values
    print('End of reading solar wind speed, IMF, AE and Dst')
    for k0 in range(16, 19):  # satellite
        pp1, pp2, pp3, pp4 = [], [], [], []
        savefn = 'F{:d}.dat'.format(k0)  # file to save data
        for year in (2011, 2012, 2013, 2014):
            for doy in np.arange(1, 367):
                # 00 means full orbit?
                fn = glob.glob('/home/guod/WD4T/ssusi/ssusi.jhuapl.edu/'
                               'data/f{:d}/apl/edr-aur/{:d}/{:03d}/'
                               '*00_DF.NC'.format(k0, year, doy))
                if not fn:  # empty fn
                    continue
                for k33, k3 in enumerate(fn):
                    print('F{:d}, {:d}-{:03d}, {:s}'.format(
                        k0, year, doy, k3[-14:-9]))
                    pp1t, pp2t, pp3t, pp4t = find_parameters_one_file_5(k3)
                    if np.size(pp1t) != 0:
                        pp1.append(pp1t)
                        pp2.append(pp2t)
                        pp3.append(pp3t)
                        pp4.append(pp4t)

                        # fig = plt.figure(figsize=(8.7, 4.4))
                        # for k44, k4 in enumerate(['N', 'S']):
                        #     ax = plt.subplot(1, 2, k44+1, polar=True)
                        #     image_energy_flux_one_file(
                        #             ax, datapath+k3, k4, s=1, vmin=0, vmax=5,
                        #             alpha=1)
                        #     fp = (pp2t.ns == k4).values
                        #     for k55, k5 in enumerate(range(int(1/dratio-1))):
                        #         r = 90-pp1t[fp, k55]
                        #         theta = pp2t.loc[fp, 'mlt']/12*np.pi
                        #         ax.scatter(theta, r, c='b', s=1, alpha=0.5)
                        #     ax.set_title(k4+'H')
                        # plt.tight_layout()
                        # plt.savefig(
                        #         savepath +
                        #         'method5/F{:d}_2013{:02d}{:02d}_{:s}.png'.
                        #         format(k0, k1, k2, k3[-14:-9]))
                        # plt.close(fig)
                    else:
                        print('    No parameters found in this file')
        pp1 = np.concatenate(pp1, axis=0)
        pp2 = pd.concat(pp2, axis=0, ignore_index=True)
        pp3 = np.concatenate(pp3, axis=0)
        pp4 = np.concatenate(pp4, axis=0)
        imfaet = imfae.reindex(pp2['datetime'], method='nearest')
        pp2['Bym'] = imfaet['Bym'].values
        pp2['Bzm'] = imfaet['Bzm'].values
        pp2['AE'] = imfaet['AE'].values
        pp2['Dst'] = imfaet['Dst'].values
        pp2['nwcf'] = imfaet['nwcf'].values
        pp2['EKL'] = imfaet['EKL'].values
        pd.to_pickle((pp1, pp2, pp3, pp4), savepath+'method5/'+savefn)
    return


def aurora_reconstruction_statistics(
        ax, sat='F16', ns='S', AErange=[0, 50], cmap='hot_r', whichp='ef'):
    """
    Reconstruct Aurora to a map
    """
    s=4
    marker='o'
    MLats, Mixps, EFs, AveEs = pd.read_pickle(
            savepath+'method5/'+'{:s}.dat'.format(sat))
    fp = ((Mixps.ns == ns) & (Mixps.AE>AErange[0]) &
          (Mixps.AE<AErange[1])).values
    mlats1, mixps1, efs1, avees1 = (
            MLats[fp, :], Mixps.loc[fp, :], EFs[fp, :], AveEs[fp, :])
    mlatsout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    psout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    if whichp == 'ef':
        vmin, vmax = 0, 6
        usep = efs1
    else:
        vmin, vmax = 0, 5
        usep = avees1
    for k00, k0 in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
        fp = (mixps1.mlt == k0).values
        if np.sum(fp) < 200:
            continue
        mlatsout[:, k00] = np.median(mlats1[fp, :], axis=0)
        psout[:, k00] = np.median(usep[fp, :], axis=0)
        ax.scatter(
                (np.ones(mlatsout[:, k00].shape)*k0)/12*np.pi,
                90-mlatsout[:, k00],
                vmin=vmin, vmax=vmax, c=psout[:, k00],
                cmap=cmap, s=s, marker=marker)
    mlatsout_ext = 0.5*(mlatsout+np.roll(mlatsout, 1, axis=1))
    psout_ext = 0.5*(psout+np.roll(psout, 1, axis=1))
    for k00, k0 in enumerate(np.arange(0, 24, BinMLT)):
        hs = ax.scatter(
                (np.ones(mlatsout_ext[:, k00].shape)*k0)/12*np.pi,
                90-mlatsout_ext[:, k00],
                vmin=vmin, vmax=vmax, c=psout_ext[:, k00],
                cmap=cmap, s=s, marker=marker)
    return ax, hs


def statistics():
    plt.close('all')
    AErange = np.ones([8, 2])
    AErange[:, 0] = np.arange(0, 400, 50)
    AErange[:, 1] = AErange[:, 0]+50

    # 1, sample number
    #    ns = 'S'
    #    sat = 'F16'
    #    p1, p2, p3, p4 = pd.read_pickle(
    #            '/home/guod/WD4T/work4/method5/{:s}.dat'.format(sat))
    #    fig, ax = plt.subplots(2, 4, figsize=(11, 6), sharex=True, sharey=True)
    #    x = np.arange(BinMLT/2, 24, BinMLT)
    #    y = np.ones(x.shape)*np.nan
    #    for k00, k0 in enumerate(range(2)):
    #        for k11, k1 in enumerate(range(4)):
    #            plt.sca(ax[k00, k11])
    #            AEranget = AErange[k00*4+k11, :]
    #            for k22, k2 in enumerate(x):
    #                fp1 = ((p2.ns==ns) & (p2.AE>AEranget[0]) &
    #                       (p2.AE<AEranget[1]) & (p2.mlt==k2)).values
    #                y[k22] = np.sum(fp1)
    #            plt.bar(x, y, align='center', width=BinMLT)
    #            plt.title('AE: [{:.0f}, {:.0f}]'.format(AEranget[0], AEranget[1]))
    #            if k00 ==1:
    #                plt.xlabel('MLT')
    #            if k1 ==0:
    #                plt.ylabel('Num')
    #    plt.ylim(0, 2500)
    #    plt.xlim(0, 24)
    #    plt.xticks(np.arange(0, 25, 4))
    #    plt.tight_layout()

    # Can we combine N and S?
    #    s = 3
    #    p1, p2, p3, p4 = pd.read_pickle(
    #            '/home/guod/WD4T/work4/method5/F17.dat')
    #    fig, ax = plt.subplots(2, 4, subplot_kw={'polar':True}, figsize=(9.7, 5.4))
    #    for k00, k0 in enumerate(range(2)):
    #        for k11, k1 in enumerate(range(4)):
    #            AEranget = AErange[k00*4+k11, :]
    #            plt.sca(ax[k00, k11])
    #            for k22, ns in enumerate(['N', 'S']):
    #                cc = 'r' if ns=='N' else 'b'
    #                for k33 , mlt in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
    #                    fp = ((p2.ns==ns) & (p2.AE>AEranget[0]) &
    #                          (p2.AE<AEranget[1]) & (p2.mlt==mlt)).values
    #                    if np.sum(fp) < 200:
    #                        continue
    #                    # plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 9]),
    #                    #             c=cc, s=s)
    #                    plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 0]),
    #                                c=cc, s=s)
    #                    plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 18]),
    #                                c=cc, s=s)
    #            mf.set_polar(plt.gca(), boundinglat=60)

    # Can we combine F16, F17 and F18?
    #    s = 3
    #    AEranget = [100, 150]
    #    fig, ax = plt.subplots(1, 2, subplot_kw={'polar':True})
    #    for k00, k0 in enumerate(['F16', 'F17', 'F18']):
    #        p1, p2, p3, p4 = pd.read_pickle(
    #                '/home/guod/WD4T/work4/method5/{:s}.dat'.format(k0))
    #        cc = ['k', 'b', 'r']
    #        for k11, ns in enumerate(['N', 'S']):
    #            plt.sca(ax[k11])
    #            for k22 , mlt in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
    #                fp = ((p2.ns==ns) & (p2.AE>AEranget[0]) &
    #                      (p2.AE<AEranget[1]) & (p2.mlt==mlt)).values
    #                if np.sum(fp) < 200:
    #                    continue
    #                #plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 9]),
    #                #            c=cc[k00], s=s)
    #                plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 0]),
    #                            c=cc[k00], s=s)
    #                plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 18]),
    #                            c=cc[k00], s=s)
    #            mf.set_polar(plt.gca(), boundinglat=60)

    # Distribution of MLats and EFs
    #    s=3
    #    ns = 'S'
    #    sat = 'F16'
    #    idratio = 9  # 0-->5%, 9-->50%, 18-->95%
    #    p1, p2, p3, p4 = pd.read_pickle(
    #            '/home/guod/WD4T/work4/method5/{:s}.dat'.format(sat))
    #    # MLat
    #    fig, ax = plt.subplots(2, 4, figsize=(11, 6), sharex=True, sharey=True)
    #    x = np.arange(BinMLT/2, 24, BinMLT)
    #    for k00, k0 in enumerate(range(2)):
    #        for k11, k1 in enumerate(range(4)):
    #            plt.sca(ax[k00, k11])
    #            AEranget = AErange[k00*4+k11, :]
    #            for k22, k2 in enumerate(x):
    #                fp1 = ((p2.ns==ns) & (p2.AE>AEranget[0]) &
    #                       (p2.AE<AEranget[1]) & (p2.mlt==k2)).values
    #                if np.sum(fp1)<200:
    #                    continue
    #                # plt.boxplot(p1[fp1, idratio], notch=True,
    #                #             positions=[k2], widths=0.6)
    #                plt.scatter(k2, np.median(p1[fp1, idratio]),
    #                            c=sns.color_palette()[0], s=s)
    #                plt.scatter(k2, np.percentile(p1[fp1, idratio], 25),
    #                            c='gray', s=s)
    #                plt.scatter(k2, np.percentile(p1[fp1, idratio], 75),
    #                            c='gray', s=s)
    #            plt.title('AE: [{:.0f}, {:.0f}]'.format(AEranget[0], AEranget[1]))
    #            if k0 ==1:
    #                plt.xlabel('MLT')
    #            if k1 ==0:
    #                plt.ylabel('MLat')
    #    plt.ylim(60, 90)
    #    plt.xlim(0, 24)
    #    plt.xticks(np.arange(0, 25, 4))
    #    plt.tight_layout()
    #    # EF
    #    fig, ax = plt.subplots(2, 4, figsize=(11, 6), sharex=True, sharey=True)
    #    x = np.arange(BinMLT/2, 24, BinMLT)
    #    for k00, k0 in enumerate(range(2)):
    #        for k11, k1 in enumerate(range(4)):
    #            plt.sca(ax[k00, k11])
    #            AEranget = AErange[k00*4+k11, :]
    #            for k22, k2 in enumerate(x):
    #                fp1 = ((p2.ns==ns) & (p2.AE>AEranget[0]) &
    #                       (p2.AE<AEranget[1]) & (p2.mlt==k2)).values
    #                if np.sum(fp1)<200:
    #                    continue
    #                # plt.boxplot(p3[fp1, idratio], notch=True,
    #                #             positions=[k2], widths=0.6)
    #                plt.scatter(k2, np.median(p3[fp1, idratio]),
    #                            c=sns.color_palette()[0], s=s)
    #                plt.scatter(k2, np.percentile(p3[fp1, idratio], 25),
    #                            c='gray', s=s)
    #                plt.scatter(k2, np.percentile(p3[fp1, idratio], 75),
    #                            c='gray', s=s)
    #            plt.title('AE: [{:.0f}, {:.0f}]'.format(AEranget[0], AEranget[1]))
    #            if k0 ==1:
    #                plt.xlabel('MLT')
    #            if k1 ==0:
    #                plt.ylabel('Energy Flux')
    #    plt.ylim(0, 12)
    #    plt.xlim(0, 24)
    #    plt.xticks(np.arange(0, 25, 4))
    #    plt.tight_layout()

    plt.show()


def epstein(x, a, b, c, d):
    return a*np.exp((x-b)/c)/((1+np.exp((x-b)/d))**2)


def epstein_fourier_fit(sat='F16', ns='S', AErange=(0, 50), n=4):
    from scipy.optimize import curve_fit
    MLats, Mixps, EFs, AveEs = pd.read_pickle(
            savepath+'method5/'+'{:s}.dat'.format(sat))
    fp = (Mixps.ns == ns) & (Mixps.AE > AErange[0]) & (Mixps.AE < AErange[1])
    fp = fp.values
    mlats1, mixps1, efs1, avees1 =(
            MLats[fp, :], Mixps.loc[fp, :], EFs[fp, :], AveEs[fp, :])
    mlatsout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    efsout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    aveesout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    for imlt, mlt in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
        fp = (mixps1.mlt==mlt).values
        if np.sum(fp) < 500:
            continue
        mlatsout[:, imlt] = np.median(mlats1[fp, :], axis=0)
        efsout[:, imlt] = np.median(efs1[fp, :], axis=0)
        aveesout[:, imlt] = np.median(avees1[fp, :], axis=0)
    p4ef = np.ones((4, int(24/BinMLT)))*np.nan
    p4avee = np.ones((4, int(24/BinMLT)))*np.nan
    for imlt, mlt in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
        if np.isnan(mlatsout[:, imlt]).all():
            continue
        p4ef[:, imlt], pcov = curve_fit(
                epstein, mlatsout[:, imlt], efsout[:, imlt],
                bounds=([0,60, 0, 0], [100, 90, 20, 20]), max_nfev=30000)
        p4avee[:, imlt], pcov = curve_fit(
                epstein, mlatsout[:, imlt], aveesout[:, imlt],
                bounds=([0,60, 0, 0], [100, 90, 20, 20]), max_nfev=30000)
    cossin = np.ones((int(24/BinMLT), (n+1)*2))*np.nan
    for k0 in np.arange(n+1):
        cossin[:, k0*2] = np.cos(k0*np.pi*np.arange(BinMLT/2, 24, BinMLT)/12)
        cossin[:, k0*2+1] = np.sin(k0*np.pi*np.arange(BinMLT/2, 24, BinMLT)/12)
    fp = ~np.isnan(p4ef[0, :])
    fgp4ef = np.linalg.lstsq(cossin[fp, :], p4ef[:, fp].T)[0]
    fp = ~np.isnan(p4avee[0, :])
    fgp4avee = np.linalg.lstsq(cossin[fp, :], p4avee[:, fp].T)[0]
    if True: # for test
        plt.close('all')
        which = 3
        mlt = np.arange(BinMLT/2, 24, BinMLT)
        plt.plot(mlt, p4ef[which, :])
        ya = np.dot(cossin, fgp4ef)[:,which]
        plt.plot(mlt, ya)
        plt.show()

    return fgp4ef, fgp4avee


def fourier_fit(AErange=(0, 50), n=4):
    # F16
    MLats16, Mixps16, EFs16, AveEs16 = pd.read_pickle(
            savepath+'method5/F16.dat')
    fp16 = ((Mixps16.AE > AErange[0]) & (Mixps16.AE < AErange[1])).values
    mlats16_1, mixps16_1, efs16_1, avees16_1 =(
            MLats16[fp16, :], Mixps16.loc[fp16, :],
            EFs16[fp16, :], AveEs16[fp16, :])
    # F17
    MLats17, Mixps17, EFs17, AveEs17 = pd.read_pickle(
            savepath+'method5/F17.dat')
    fp17 = ((Mixps17.AE > AErange[0]) & (Mixps17.AE < AErange[1])).values
    mlats17_1, mixps17_1, efs17_1, avees17_1 =(
            MLats17[fp17, :], Mixps17.loc[fp17, :],
            EFs17[fp17, :], AveEs17[fp17, :])
    # combine
    mlats1 = np.concatenate((mlats16_1, mlats17_1), axis=0)
    mixps1 = pd.concat((mixps16_1, mixps17_1), axis=0, ignore_index=True)
    efs1 = np.concatenate((efs16_1, efs17_1), axis=0)
    avees1 = np.concatenate((avees16_1, avees17_1), axis=0)

    mlatsout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    efsout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    aveesout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    for imlt, mlt in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
        fp = (mixps1.mlt==mlt).values
        if np.sum(fp) < 200:
            continue
        mlatsout[:, imlt] = np.median(mlats1[fp, :], axis=0)
        efsout[:, imlt] = np.median(efs1[fp, :], axis=0)
        aveesout[:, imlt] = np.median(avees1[fp, :], axis=0)
    cossin = np.ones((int(24/BinMLT), 2*(n+1)))*np.nan
    for k0 in np.arange(n+1):
        cossin[:, k0*2] = np.cos(k0*np.pi*np.arange(BinMLT/2, 24, BinMLT)/12)
        cossin[:, k0*2+1] = np.sin(k0*np.pi*np.arange(BinMLT/2, 24, BinMLT)/12)
    fp = ~np.isnan(mlatsout[0, :])
    fg4mlat = np.linalg.lstsq(cossin[fp, :], mlatsout[:, fp].T)[0]
    #mlat50p = mlatsout[9, :]  # 50% of total energy
    #fg4mlat50p = np.linalg.lstsq(cossin[fp, :], mlat50p[fp])[0]
    #dmlatsout = mlatsout - mlat50p.reshape(1, -1)
    #fg4dmlat = np.linalg.lstsq(cossin[fp, :], dmlatsout[:, fp].T)[0]

    #fg4mlat = fg4dmlat + fg4mlat50p.reshape(-1, 1)
    fg4ef = np.linalg.lstsq(cossin[fp, :], efsout[:, fp].T)[0]
    fg4avee = np.linalg.lstsq(cossin[fp, :], aveesout[:, fp].T)[0]
    if True:   # for test
        plt.close('all')
        mlt = np.arange(0, 24, 0.01)
        cossint = np.ones((len(mlt), 2*(n+1)))*np.nan
        for k0 in np.arange(n+1):
            cossint[:, k0*2] = np.cos(k0*np.pi*mlt/12)
            cossint[:, k0*2+1] = np.sin(k0*np.pi*mlt/12)
        mlatst = np.dot(cossint, fg4mlat)
        efst = np.dot(cossint, fg4ef)
        aveest = np.dot(cossint, fg4avee)

        fig = plt.figure()
        ax = plt.subplot()
        for k0 in np.arange(10, 15):
            imlt, = np.where(mlt==k0)
            plt.plot(mlatst[imlt, :].reshape(-1), efst[imlt, :].reshape(-1))
        plt.xlim(75, 80)
        plt.ylim(0, 5)
        fig = plt.figure()
        ax = plt.subplot()
        plt.plot(mlt, mlatst)
        fig = plt.figure()
        ax = plt.subplot(polar=True)
        for k0 in np.arange(int(1/dratio-1)):
            plt.scatter(mlt/12*np.pi,
                        90-mlatst[:, k0], s=3, c = aveest[:, k0],
                        vmin=0, vmax=6, cmap='hot_r')
        mf.set_polar(ax, boundinglat=60)
        plt.show()
    return


def aurora_reconstruction_fit(
        sat='F16', ns='S', AErange=(0, 50), n=4, whichp='ef'):
    fgp4ef, fgp4avee = epstein_fourier_fit(sat, ns, AErange, n)
    fg = fgp4ef if whichp=='ef' else fgp4avee
    mlt = np.arange(0, 24.1, 1)
    cossin = np.ones((len(mlt), (n+1)*2))*np.nan
    for k0 in np.arange(n+1):
        cossin[:, k0*2] = np.cos(k0*np.pi*mlt/12)
        cossin[:, k0*2+1] = np.sin(k0*np.pi*mlt/12)
    mabcd = np.dot(cossin, fg)
    mlat = np.arange(60, 90)
    yout = np.ones((len(mlat), len(mlt)))*np.nan
    for imlt in np.arange(len(mlt)):
        yout[:, imlt] = epstein(mlat, mabcd[imlt, 0], mabcd[imlt, 1],
                                mabcd[imlt, 2], mabcd[imlt, 3])
    mlt, mlat = np.meshgrid(mlt, mlat)
    plt.close('all')
    ax = plt.subplot(polar=True)
    plt.contourf(mlt/12*np.pi, 90-mlat, yout, cmap='hot_r',
                 levels=np.arange(100), zorder=0, extend='both')
    mf.set_polar(ax, boundinglat=60)
    plt.show()


def linear_fit():
    from scipy import stats

    sat = 'F18'

    MLT = np.arange(BinMLT/2, 24, BinMLT)
    # NS, MLT, 5%-95%, (mlat, ef, average energy)
    slope = np.ones((2, len(MLT), 19, 3))*np.nan
    intercept = np.ones((2, len(MLT), 19, 3))*np.nan
    rvalue = np.ones((2, len(MLT), 19, 3))*np.nan

    mlats, p2, efs, aes = pd.read_pickle(
            '/home/guod/WD4T/work4/method5/{:s}.dat'.format(sat))
    for k00, k0 in enumerate(['N', 'S']):
        for k11, k1 in enumerate(MLT):
            fp = ((p2.mlt == k1) & (p2.ns == k0)).values
            if np.sum(fp) < 150:
                continue
            mlats0, efs0, aes0 = (k[fp, :] for k in [mlats, efs, aes])
            AE = (p2.loc[fp, 'AE']).values
            for k22, k2 in enumerate(range(19)):
                mlats1, efs1, aes1 = (k[:, k22] for k in [mlats0, efs0, aes0])
                for k33, k3 in enumerate([mlats1, efs1, aes1]):
                    slopet, interceptt, rvaluet, pvalue, stderr = \
                            stats.linregress(AE, k3)
                    if rvaluet>0.4:
                        print(slopet, interceptt, rvaluet)
                    slope[k00, k11, k22, k33] = slopet
                    intercept[k00, k11, k22, k33] = interceptt
                    rvalue[k00, k11, k22, k33] = rvaluet
    pd.to_pickle([slope, intercept, rvalue],
                 savepath+'/method5/{:s}_fit.dat'.format(sat))
    return


def test_linear_fit():
    plt.close('all')
    sat='F16'
    MLT = np.arange(BinMLT/2, 24, BinMLT)
    slope, intercept, rvalue = pd.read_pickle(
            savepath+'/method5/{:s}_fit.dat'.format(sat))
    fig, ax = plt.subplots(3, 3, sharex=True)
    yl = ['Slope', 'Intercept', 'r-value']
    for k00, k0 in enumerate([slope, intercept, rvalue]):
        for k11, k1 in enumerate(['MLat', 'Energy Flux', 'Average Energy']):
            plt.sca(ax[k00, k11])
            plt.plot(MLT, k0[1, :, 2, k11])
            if k00 == 0:
                plt.title(k1)
            if k00 == 2:
                plt.xlabel('MLT')
            if k11 == 0:
                plt.ylabel(yl[k00])
    plt.xlim(0, 24)
    plt.xticks(range(0, 25, 4))
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    # plt.close('all')
    # fn = (datapath+'PS.APL_V0105S024CE0018_SC.U_DI.A_GP.F16'
    #       '-SSUSI_PA.APL-EDR-AURORA_DD.20130102_SN.47516-00_DF.NC')
    # test_find_parameters_one_file_5(fn)
    # plt.show()

    # find_parameters()

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

    # Aurora reconstruction
    # plt.close('all')
    # fig, ax = plt.subplots(2, 4, subplot_kw={'polar': True}, figsize=(14, 8))
    # AErange = np.ones([8, 2])
    # AErange[:, 0] = np.arange(0, 400, 50)
    # AErange[:, 1] = AErange[:, 0]+50
    # whichp = 'avees' if False else 'ef'
    # for k0 in range(2):
    #     for k1 in range(4):
    #         axt = ax[k0, k1]
    #         a, b = aurora_reconstruction_statistics(
    #                 axt, AErange=AErange[k0*4+k1, :], whichp=whichp)
    #         mf.set_polar(axt, boundinglat=60)
    #         axt.set_title('AE: [{:.0f}, {:.0f}]'.\
    #                 format(AErange[k0*4+k1, 0], AErange[k0*4+k1, 1]), y=1.06)
    # cax = plt.axes([0.35, 0.06, 0.3, 0.03])
    # plt.colorbar(b, cax=cax, orientation='horizontal')
    # if whichp == 'ef':
    #     plt.xlabel('Energy flux (ergs/s/cm$^2$)')
    # else:
    #     plt.xlabel('Mean Energy (keV)')
    # plt.show()

    # plt.close('all')
    # fn = (datapath+'PS.APL_V0105S024CE0018_SC.U_DI.A_GP.F16'
    #       '-SSUSI_PA.APL-EDR-AURORA_DD.20130118_SN.47748-00_DF.NC')
    # ax = plt.subplot(polar=True)
    # #ax, hs = image_energy_flux_one_file(
    # #    ax, fn, 'S', vmin=0, vmax=10, s=1, alpha=1, kernel_size=1)
    # ax, hs = test_find_parameters_one_file_5(fn)
    # plt.colorbar(hs)
    # plt.show()

    # linear_fit()

    # test_linear_fit()

    # statistics()

    # a, b = epstein_fourier_fit()
    # aurora_reconstruction_fit()

    fourier_fit()

    print('----------------------------------------------------------------')
    print('Remaining problems:')
    print('  1, Aurora with complicated structure')
    print('  2, Different total energy, different proportion?')
    print('  3, Remove dayglow?')
    print('  4, Use median or mean?')
    print('----------------------------------------------------------------')
