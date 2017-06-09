# Auroral precipation model using SSUSI data
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import myfunctions as mf
import scipy.signal as signal
import seaborn as sns
import glob
import os
from scipy import stats
sns.set('paper', 'whitegrid')

# global variables
LLP1 = 200  # Lower limit of pixels in one bin
LLEFTOTH = 200  # lower limit of total energy flux in one hemisphere
BinMLT = 0.25  # MLT range in one bin
LLEFTOT = 10  # lower limit of total energy flux in each bin, for method 1,2
datapath = '/home/guod/WD4T/ssusi/'
savepath = '/home/guod/WD4T/work4/'
dratio = 0.05 # flag 5%, 10%, ..., 95% of the total energy flux
LLSMP = 200 # lower limit of samples
n = 6 # number of fourier fitting terms
RE = 6371000 # unit: m
AErange = np.ones([8, 2])
AErange[:, 0] = np.arange(0, 400, 50)
AErange[:, 1] = AErange[:, 0]+50
PIXELSIZE_GEOMAGNETIC_LATITUDE=0.22077596187591553 # unit: degree
PIXELSIZE_GEOMAGNETIC_LONGITUDE=0.22077596187591553


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
                    vmax=vmax, alpha=alpha, cmap='nipy_spectral')
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


def find_parameters_one_file_5(fn, test=False):
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
    # kernel_size_own = 3
    # energyfluxn = signal.medfilt2d(energyfluxn, kernel_size_own)
    # energyfluxs = signal.medfilt2d(energyfluxs, kernel_size_own)
    # meanenergyn = signal.medfilt2d(meanenergyn, kernel_size_own)
    # meanenergys = signal.medfilt2d(meanenergys, kernel_size_own)

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
    # pp1: MLats; pp2: other parameters; pp3: ef; pp4: mean energy
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
            llmlat, ulmlat = 65, 75
        else:
            sf = 0.90
            fpt = ef00 >= (0.5-sf/2)*ef00[-1]
            llmlat = mlat00[fpt][0]
            fpt = ef00 <= (0.5+sf/2)*ef00[-1]
            ulmlat = mlat00[fpt][-1]

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

            # Mlats of 5%, 10%... energy flux -->pp1
            # energy fluxes at 5%, 10%, ... --> pp3
            # mean energy --> pp4
            pp1t = np.ones(int(1/dratio-1))*np.nan
            pp3t = np.ones(int(1/dratio-1))*np.nan
            pp4t = np.ones(int(1/dratio-1))*np.nan
            utt = np.ones(int(1/dratio-1))*np.nan
            for k22, k2 in enumerate(np.arange(dratio, 1, dratio)):
                fp = (efcum >= eftot*k2)
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
    if test: # for test
        plt.figure(figsize=(8.7, 4.4))
        for k00, k0 in enumerate(['N', 'S']):
            ax = plt.subplot(1, 2, k00+1, polar=True)
            ax, hs = image_energy_flux_one_file(
                    ax, fn, k0, vmin=0, vmax=10, s=1, alpha=1)
            fp = (pp2.ns == k0).values
            for k11, k1 in enumerate(range(int(1/dratio-1))):
                r = 90-pp1[fp, k11]
                theta = pp2.loc[fp, 'mlt']/12*np.pi
                ax.scatter(theta, r, c='r', s=10, alpha=0.5)
            plt.title(k0+'H')
        plt.tight_layout()
    return pp1, pp2, pp3, pp4  # Mlats, ..., EFs, MEs


def find_parameters(saveplot=False):
    import os
    import fnmatch
    import omni
    # Read IMF By, Bz and AE
    print('Begin to read solar wind speed, IMF, AE and Dst')
    imfae = omni.get_omni(
            bdate='2011-1-1', edate='2015-1-1',
            variables=['Bym', 'Bzm', 'AE', 'V'], res='1m')
    imfae['Bt'] = np.sqrt(imfae['Bym']**2 + imfae['Bzm']**2)
    imfae['nwcf'] = (imfae['V']**(4/3)) * (imfae['Bt']**(2/3)) * \
        (((1-imfae['Bzm']/imfae['Bt'])/2)**(4/3))/100.0
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

                        if saveplot:
                            fig = plt.figure(figsize=(8.7, 4.4))
                            for k44, k4 in enumerate(['N', 'S']):
                                ax = plt.subplot(1, 2, k44+1, polar=True)
                                image_energy_flux_one_file(
                                        ax, datapath+k3, k4, s=1, vmin=0,
                                        vmax=5, alpha=1)
                                fp = (pp2t.ns == k4).values
                                for k55, k5 in enumerate(range(int(1/dratio-1))):
                                    r = 90-pp1t[fp, k55]
                                    theta = pp2t.loc[fp, 'mlt']/12*np.pi
                                    ax.scatter(theta, r, c='b', s=1, alpha=0.5)
                                ax.set_title(k4+'H')
                            plt.tight_layout()
                            plt.savefig(
                                    savepath +
                                    'method5/F{:d}_2013{:02d}{:02d}_{:s}.png'.
                                    format(k0, k1, k2, k3[-14:-9]))
                            plt.close(fig)
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
        ax, AErange=(0, 50), cmap='viridis', whichp='ef', vmin=0, vmax=5):
    """
    Reconstruct Aurora to a map
    """
    s=5
    marker='o'
    # F16
    MLats, Mixps, EFs, AveEs = pd.read_pickle(savepath+'method5/F16.dat')
    fpmlat = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1]) &
          (Mixps.eftot!=0)).values
    fpef = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1])).values
    mlats16_1, mixps16_1, mixps16_2, efs16_2, avees16_1 = (
            MLats[fpmlat, :], Mixps.loc[fpmlat, :],
            Mixps.loc[fpef, :], EFs[fpef, :], AveEs[fpmlat, :])
    # F17
    MLats, Mixps, EFs, AveEs = pd.read_pickle(savepath+'method5/F17.dat')
    fpmlat = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1]) &
          (Mixps.eftot!=0)).values
    fpef = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1])).values
    mlats17_1, mixps17_1, mixps17_2, efs17_2, avees17_1 = (
            MLats[fpmlat, :], Mixps.loc[fpmlat, :],
            Mixps.loc[fpef, :], EFs[fpef, :], AveEs[fpmlat, :])
    # combine
    mlats1 = np.concatenate((mlats16_1, mlats17_1), axis=0)
    mixps1 = pd.concat((mixps16_1, mixps17_1), axis=0, ignore_index=True)
    mixps2 = pd.concat((mixps16_2, mixps17_2), axis=0, ignore_index=True)
    efs2 =np.concatenate((efs16_2, efs17_2), axis=0)
    avees1 =np.concatenate((avees16_1, avees17_1), axis=0)

    mlatsout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    psout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    if whichp == 'ef':
        usep = efs2
    else:
        usep = avees1
    for k00, k0 in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
        if (k0>=10) & (k0<=14):
            continue
        fp1 = (mixps1.mlt == k0).values
        fp2 = (mixps2.mlt == k0).values
        if np.sum(fp1) < LLSMP:
            continue
        mlatsout[:, k00] = np.median(mlats1[fp1, :], axis=0)
        if whichp == 'ef':
            psout[:, k00] = np.median(usep[fp2, :], axis=0)
        else:
            psout[:, k00] = np.median(usep[fp1, :], axis=0)
        hs = ax.scatter(
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


def aurora_reconstruction_statistics_all(whichp='ef', cmap='nipy_spectral',
                                         vmin=0, vmax=5):
    plt.close('all')
    fig, ax = plt.subplots(2, 4, subplot_kw={'polar': True}, figsize=(14, 8))
    AErange = np.ones([8, 2])
    AErange[:, 0] = np.arange(0, 400, 50)
    AErange[:, 1] = AErange[:, 0]+50
    for k0 in range(2):
        for k1 in range(4):
            axt = ax[k0, k1]
            a, b = aurora_reconstruction_statistics(
                    axt, AErange=AErange[k0*4+k1, :],
                    whichp=whichp, cmap=cmap, vmin=vmin, vmax=vmax)
            mf.set_polar(axt, boundinglat=60)
            axt.set_title('AE: [{:.0f}, {:.0f}]'.\
                    format(AErange[k0*4+k1, 0], AErange[k0*4+k1, 1]), y=1.06)
    cax = plt.axes([0.35, 0.06, 0.3, 0.03])
    plt.colorbar(b, cax=cax, orientation='horizontal')
    if whichp == 'ef':
        plt.xlabel('Energy flux (ergs/s/cm$^2$)')
    else:
        plt.xlabel('Mean Energy (keV)')
    plt.show()


def statistics():
    plt.close('all')

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
    #                       (p2.AE<AEranget[1]) & (p2.mlt==k2) &
    #                       (p2.eftot>LLEFTOT)).values
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
    # s = 3
    # p1, p2, p3, p4 = pd.read_pickle(
    #         '/home/guod/WD4T/work4/method5/F17.dat')
    # fig, ax = plt.subplots(2, 4, subplot_kw={'polar':True}, figsize=(9.7, 5.4))
    # for k00, k0 in enumerate(range(2)):
    #     for k11, k1 in enumerate(range(4)):
    #         AEranget = AErange[k00*4+k11, :]
    #         plt.sca(ax[k00, k11])
    #         for k22, ns in enumerate(['N', 'S']):
    #             cc = 'r' if ns=='N' else 'b'
    #             for k33 , mlt in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
    #                 fp = ((p2.ns==ns) & (p2.AE>AEranget[0]) &
    #                       (p2.AE<AEranget[1]) & (p2.mlt==mlt) &
    #                       (p2.eftot!=0)).values
    #                 if np.sum(fp) < LLSMP:
    #                     continue
    #                 # plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 9]),
    #                 #             c=cc, s=s)
    #                 plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 0]),
    #                             c=cc, s=s)
    #                 plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 18]),
    #                             c=cc, s=s)
    #         mf.set_polar(plt.gca(), boundinglat=60)

    # Can we combine F16, F17 and F18?
    # s = 3
    # AEranget = [100, 150]
    # fig, ax = plt.subplots(1, 2, subplot_kw={'polar':True})
    # for k00, k0 in enumerate(['F16', 'F17', 'F18']):
    #     p1, p2, p3, p4 = pd.read_pickle(
    #             '/home/guod/WD4T/work4/method5/{:s}.dat'.format(k0))
    #     cc = ['k', 'b', 'r']
    #     for k11, ns in enumerate(['N', 'S']):
    #         plt.sca(ax[k11])
    #         for k22 , mlt in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
    #             fp = ((p2.ns==ns) & (p2.AE>AEranget[0]) &
    #                   (p2.AE<AEranget[1]) & (p2.mlt==mlt) &
    #                   (p2.eftot>LLEFTOT)).values
    #             if np.sum(fp) < LLSMP:
    #                 continue
    #             #plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 9]),
    #             #            c=cc[k00], s=s)
    #             plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 0]),
    #                         c=cc[k00], s=s)
    #             plt.scatter(mlt/12*np.pi, 90-np.median(p1[fp, 18]),
    #                         c=cc[k00], s=s)
    #         mf.set_polar(plt.gca(), boundinglat=60)

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
    #                       (p2.AE<AEranget[1]) & (p2.mlt==k2) &
    #                       (p2.eftot>LLEFTOT)).values
    #                if np.sum(fp1)<LLSMP:
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
    #                       (p2.AE<AEranget[1]) & (p2.mlt==k2) &
    #                       (p2.eftot>LLEFTOT)).values
    #                if np.sum(fp1)<LLSMP:
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

    # What are discarded
    AErange = (200, 250)
    # F16
    MLats, Mixps, EFs, AveEs = pd.read_pickle(savepath+'method5/F16.dat')
    fp = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1])).values
    mlats16_1, mixps16_1, efs16_1, avees16_1 = (
            MLats[fp, :], Mixps.loc[fp, :], EFs[fp, :], AveEs[fp, :])
    # F17
    MLats, Mixps, EFs, AveEs = pd.read_pickle(savepath+'method5/F17.dat')
    fp = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1])).values
    mlats17_1, mixps17_1, efs17_1, avees17_1 = (
            MLats[fp, :], Mixps.loc[fp, :], EFs[fp, :], AveEs[fp, :])
    # combine
    mlats1 = np.concatenate((mlats16_1, mlats17_1), axis=0)
    mixps1 = pd.concat((mixps16_1, mixps17_1), axis=0, ignore_index=True)
    efs1 =np.concatenate((efs16_1, efs17_1), axis=0)
    avees1 =np.concatenate((avees16_1, avees17_1), axis=0)

    fig, ax = plt.subplots(4, 6, sharex=True, sharey=True)
    for imlt, mlt in enumerate(np.arange(BinMLT/2, 24, 1)):
        plt.sca(ax[imlt//6, imlt%6])
        fp = (mixps1.mlt==mlt)
        xx = mixps1.loc[fp, 'eftot'].values
        print(mlt, np.sum(xx==0)/xx.size)
        for p in np.arange(0, 200, LLEFTOT):
            xxt = np.sum((xx>=p) & (xx<p+LLEFTOT))/xx.size
            plt.bar(p+LLEFTOT/2, xxt, width=LLEFTOT)
        plt.xlim(0, 100)
        plt.ylim(0, 1)
        plt.title('MLT={:6.3f}'.format(mlt))
    [ax[k, 0].set_ylabel('P') for k in range(4)]
    [ax[-1, k].set_xlabel('Sum of EF') for k in range(6)]

    plt.show()


def epstein(x, a, b, c, d):
    return a*np.exp((x-b)/c)/((1+np.exp((x-b)/d))**2)


def fourier_fit(AErange=(100, 150), test=False):
    # F16
    MLats, Mixps, EFs, AveEs = pd.read_pickle(savepath+'method5/F16.dat')
    fpmlat = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1]) &
          (Mixps.eftot!=0)).values
    fpef = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1])).values
    mlats16_1, mixps16_1, mixps16_2, efs16_2, avees16_1 = (
            MLats[fpmlat, :], Mixps.loc[fpmlat, :],
            Mixps.loc[fpef, :], EFs[fpef, :], AveEs[fpmlat, :])
    # F17
    MLats, Mixps, EFs, AveEs = pd.read_pickle(savepath+'method5/F17.dat')
    fpmlat = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1]) &
          (Mixps.eftot!=0)).values
    fpef = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1])).values
    mlats17_1, mixps17_1, mixps17_2, efs17_2, avees17_1 = (
            MLats[fpmlat, :], Mixps.loc[fpmlat, :],
            Mixps.loc[fpef, :], EFs[fpef, :], AveEs[fpmlat, :])
    # combine
    mlats1 = np.concatenate((mlats16_1, mlats17_1), axis=0)
    mixps1 = pd.concat((mixps16_1, mixps17_1), axis=0, ignore_index=True)
    mixps2 = pd.concat((mixps16_2, mixps17_2), axis=0, ignore_index=True)
    efs2 =np.concatenate((efs16_2, efs17_2), axis=0)
    avees1 =np.concatenate((avees16_1, avees17_1), axis=0)

    mlatsout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    efsout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    aveesout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    for imlt, mlt in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
        if (mlt>10) & (mlt<14):
            continue
        fp1 = (mixps1.mlt==mlt).values
        fp2 = (mixps2.mlt==mlt).values
        if np.sum(fp1) < LLSMP:
            continue
        mlatsout[:, imlt] = np.median(mlats1[fp1, :], axis=0)
        efsout[:, imlt] = np.median(efs2[fp2, :], axis=0)
        aveesout[:, imlt] = np.median(avees1[fp1, :], axis=0)

    # complement data gap
    cossin = np.ones((int(24/BinMLT), 2*(1+1)))*np.nan
    for k0 in np.arange(1+1):
        cossin[:, k0*2] = np.cos(k0*np.pi*np.arange(BinMLT/2, 24, BinMLT)/12)
        cossin[:, k0*2+1] = np.sin(k0*np.pi*np.arange(BinMLT/2, 24, BinMLT)/12)
    fp = ~np.isnan(mlatsout[0, :])
    fg4mlat = np.linalg.lstsq(cossin[fp, :], mlatsout[:, fp].T)[0]
    fg4ef = np.linalg.lstsq(cossin[fp, :], efsout[:, fp].T)[0]
    fg4avee = np.linalg.lstsq(cossin[fp, :], aveesout[:, fp].T)[0]

    mlatsout1 = np.dot(cossin, fg4mlat).T
    efsout1 = np.dot(cossin, fg4ef).T
    aveesout1 = np.dot(cossin, fg4avee).T

    mlatsout[:, ~fp] = mlatsout1[:, ~fp]
    efsout[:, ~fp] = efsout1[:, ~fp]
    aveesout[:, ~fp] = aveesout1[:, ~fp]

    # Fit
    cossin = np.ones((int(24/BinMLT), 2*(n+1)))*np.nan
    for k0 in np.arange(n+1):
        cossin[:, k0*2] = np.cos(k0*np.pi*np.arange(BinMLT/2, 24, BinMLT)/12)
        cossin[:, k0*2+1] = np.sin(k0*np.pi*np.arange(BinMLT/2, 24, BinMLT)/12)
    fg4mlat = np.linalg.lstsq(cossin, mlatsout.T)[0]
    fg4ef = np.linalg.lstsq(cossin, efsout.T)[0]
    fg4avee = np.linalg.lstsq(cossin, aveesout.T)[0]

    if test:   # for test
        plt.close('all')
        mlt = np.arange(0, 24, 0.05)
        cossint = np.ones((len(mlt), 2*(n+1)))*np.nan
        for k0 in np.arange(n+1):
            cossint[:, k0*2] = np.cos(k0*np.pi*mlt/12)
            cossint[:, k0*2+1] = np.sin(k0*np.pi*mlt/12)
        mlatst = np.dot(cossint, fg4mlat)
        efst = np.dot(cossint, fg4ef)
        aveest = np.dot(cossint, fg4avee)

        fig = plt.figure()
        ax = plt.subplot()
        for k0 in np.arange(9, 16):
            imlt, = np.where(mlt==k0)
            plt.plot(mlatst[imlt, :].reshape(-1), efst[imlt, :].reshape(-1))
        plt.xlim(70, 80)
        plt.ylim(0, 5)
        fig = plt.figure()
        ax = plt.subplot()
        plt.plot(mlt, mlatst)
        fig = plt.figure()
        ax = plt.subplot(polar=True)
        for k0 in np.arange(int(1/dratio-1)):
            plt.scatter(mlt/12*np.pi,
                        90-mlatst[:, k0], s=3, c = efst[:, k0],
                        vmin=0, vmax=6, cmap='viridis')
        mf.set_polar(ax, boundinglat=60)
        plt.show()
    return fg4mlat, fg4ef, fg4avee


def fourier_fit_save_parameters():
    fg = []
    for k0 in range(8):
        print(n, ', ', AErange[k0, :])
        fg4mlat, fg4ef, fg4avee = fourier_fit(AErange=AErange[k0, :])
        fg.append(np.stack([fg4mlat, fg4ef, fg4avee], axis=2))
    fg = np.stack(fg, axis=3)
    pd.to_pickle(fg, '/home/guod/WD4T/work4/method5/fg.dat')
    return


def aurora_reconstruction_fit(
        ax, AErange=(0, 50), whichp='ef', vmin=0, vmax=8, cmap='viridis', s=3):
    iAE = int(AErange[0]/50)
    fg = pd.read_pickle('/home/guod/WD4T/work4/method5/fg.dat')
    fg = fg[:, :, :, iAE]
    mlt = np.arange(0, 24, 0.1)
    cossin = np.ones((len(mlt), (n+1)*2))*np.nan
    for k0 in np.arange(n+1):
        cossin[:, k0*2] = np.cos(k0*np.pi*mlt/12)
        cossin[:, k0*2+1] = np.sin(k0*np.pi*mlt/12)
    mlat = np.dot(cossin, fg[:, :, 0])
    ef = np.dot(cossin, fg[:, :, 1])
    avee = np.dot(cossin, fg[:, :, 2])
    yy = ef if whichp == 'ef' else avee
    for imlt, mltv in enumerate(mlt):
        hc =ax.scatter(np.ones(mlat.shape[1])*mltv/12*np.pi, 90-mlat[imlt, :],
                       c=yy[imlt, :], s=s, vmin=vmin, vmax=vmax, cmap=cmap)
    return ax, hc

def aurora_reconstruction_fit_all(whichp='avee', cmap='nipy_spectral'):
    plt.close('all')
    fig, ax = plt.subplots(
            2, 4, subplot_kw={'polar': True}, figsize=(10.4, 6.5))
    AErange = np.ones([8, 2])
    AErange[:, 0] = np.arange(0, 400, 50)
    AErange[:, 1] = AErange[:, 0]+50
    for k0 in range(2):
        for k1 in range(4):
            axt = ax[k0, k1]
            axt, hc = aurora_reconstruction_fit(
                    axt, AErange=AErange[k0*4+k1, :], whichp=whichp,
                    vmin=0, vmax=5, cmap=cmap)
            axt.set_title('AE: [{:.0f}, {:.0f}]'.\
                    format(AErange[k0*4+k1, 0], AErange[k0*4+k1, 1]), y=1.06)
            mf.set_polar(axt, boundinglat=60)
    cax = plt.axes([0.35, 0.06, 0.3, 0.03])
    plt.colorbar(hc, cax=cax, orientation='horizontal')
    if whichp == 'ef':
        plt.xlabel('Energy flux (ergs/s/cm$^2$)')
    else:
        plt.xlabel('Mean Energy (keV)')
    plt.show()
    return


def diff_fit_statistics(
        ax, AErange=(0, 50), cmap='viridis', whichp='ef', vmin=0, vmax=5, s=5):
    # statistics
    # F16
    MLats, Mixps, EFs, AveEs = pd.read_pickle(savepath+'method5/F16.dat')
    fp = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1])).values
    mlats16_1, mixps16_1, efs16_1, avees16_1 = (
            MLats[fp, :], Mixps.loc[fp, :], EFs[fp, :], AveEs[fp, :])
    # F17
    MLats, Mixps, EFs, AveEs = pd.read_pickle(savepath+'method5/F17.dat')
    fp = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1])).values
    mlats17_1, mixps17_1, efs17_1, avees17_1 = (
            MLats[fp, :], Mixps.loc[fp, :], EFs[fp, :], AveEs[fp, :])
    # combine
    mlats1 = np.concatenate((mlats16_1, mlats17_1), axis=0)
    mixps1 = pd.concat((mixps16_1, mixps17_1), axis=0, ignore_index=True)
    efs1 =np.concatenate((efs16_1, efs17_1), axis=0)
    avees1 =np.concatenate((avees16_1, avees17_1), axis=0)

    mlatsout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    psout = np.ones((int(1/dratio-1), int(24/BinMLT)))*np.nan
    if whichp == 'ef':
        usep = efs1
    else:
        usep = avees1
    mlt = np.arange(BinMLT/2, 24, BinMLT)
    for k00, k0 in enumerate(mlt):
        fp = (mixps1.mlt == k0).values
        if np.sum(fp) < LLSMP:
            continue
        mlatsout[:, k00] = np.median(mlats1[fp, :], axis=0)
        psout[:, k00] = np.median(usep[fp, :], axis=0)
    fp = ~np.isnan(mlatsout[0, :])
    mlt = mlt[fp]
    mlatsout = mlatsout[:, fp].T
    psout = psout[:, fp].T

    # fit
    iAE = int(AErange[0]/50)
    fg = pd.read_pickle('/home/guod/WD4T/work4/method5/fg.dat')
    fg = fg[:, :, :, iAE]
    cossin = np.ones((len(mlt), (n+1)*2))*np.nan
    for k0 in np.arange(n+1):
        cossin[:, k0*2] = np.cos(k0*np.pi*mlt/12)
        cossin[:, k0*2+1] = np.sin(k0*np.pi*mlt/12)
    mlatfit = np.dot(cossin, fg[:, :, 0])
    effit = np.dot(cossin, fg[:, :, 1])
    aveefit = np.dot(cossin, fg[:, :, 2])
    psfit = effit if whichp == 'ef' else aveefit

    for imlt, mltv in enumerate(mlt):
        hc =ax.scatter(
                np.ones(mlatfit.shape[1])*mltv/12*np.pi, 90-mlatfit[imlt, :],
                c=(100*(psfit-psout)/psout)[imlt, :],
                s=s, vmin=vmin, vmax=vmax, cmap=cmap)
    return ax, hc


def hp(test=True):
    dmlt = 0.01
    mlt = np.arange(0, 24, dmlt)
    cossin = np.ones((len(mlt), (n+1)*2))*np.nan
    for k0 in np.arange(n+1):
        cossin[:, k0*2] = np.cos(k0*np.pi*mlt/12)
        cossin[:, k0*2+1] = np.sin(k0*np.pi*mlt/12)
    fg = pd.read_pickle('/home/guod/WD4T/work4/method5/fg.dat')
    HP = np.arange(8)*np.nan
    for iAE in np.arange(AErange.shape[0]):
        fgt = fg[:, :, :, iAE]
        mlatfit = np.dot(cossin, fgt[:, :, 0])
        effit = np.dot(cossin, fgt[:, :, 1])
        aveefit = np.dot(cossin, fgt[:, :, 2])
        dmlat = np.diff(mlatfit, axis=1)
        if np.sum(dmlat<=0)>=1:
            print('Something wront')
            return
        dmlat = np.append(dmlat, dmlat[:, -1].reshape(-1, 1), axis=1)
        HP[iAE] = np.sum(
                (effit*1e-3)*(RE**2)*np.cos(mlatfit*np.pi/180)*
                (dmlat*np.pi/180)*(dmlt*np.pi/12)*1e-9)  # unit: kw
    if test:
        plt.close('all')
        AE = np.mean(AErange, axis=1)
        plt.plot(AE, HP, '-o')
        #plt.plot(AE, np.sqrt(AE)+2)
        plt.show()
        plt.xlim(0, 400)
        plt.ylim(0, 25)
        plt.xlabel('AE')
        plt.ylabel('Hemisphere Power (GW)')
    return HP



def linear_fit_AE(test=True):
    dmlt = 1
    mlt = np.arange(0, 24, dmlt)
    cossin = np.ones((len(mlt), (n+1)*2))*np.nan
    for k0 in np.arange(n+1):
        cossin[:, k0*2] = np.cos(k0*np.pi*mlt/12)
        cossin[:, k0*2+1] = np.sin(k0*np.pi*mlt/12)
    fg = pd.read_pickle('/home/guod/WD4T/work4/method5/fg.dat')
    mlatsout = np.ones((len(mlt), int(1/dratio-1), AErange.shape[0]))*np.nan
    efsout = np.ones((len(mlt), int(1/dratio-1), AErange.shape[0]))*np.nan
    aveesout = np.ones((len(mlt), int(1/dratio-1), AErange.shape[0]))*np.nan
    for iAE in np.arange(AErange.shape[0]):
        fgt = fg[:, :, :, iAE]
        mlatsout[:, :, iAE] = np.dot(cossin, fgt[:, :, 0])
        efsout[:, :, iAE] = np.dot(cossin, fgt[:, :, 1])
        aveesout[:, :, iAE] = np.dot(cossin, fgt[:, :, 2])
    if test: # Test
        plt.close('all')
        fig, ax = plt.subplots(4, 6, sharex=True, sharey=True)
        for k in range(24):
            plt.sca(ax[k//6, k%6])
            plt.plot(AErange.mean(axis=1),
                     mlatsout[k, 1::8,  :].T, 'o-', lw=1.5)
            plt.xlim(0, 400)
            plt.ylim(60, 80)
            plt.title('MLT: {:02d}'.format(k))
        # ax[1, -1].legend([2, 6, 10, 14, 18, 22],
        #                  bbox_to_anchor=(1, 1.05), frameon=True)
        [ax[-1, k].set_xlabel('AE') for k in range(6)]
        [ax[k, 0].set_ylabel('MLat') for k in range(4)]
        plt.show()
    return


def imf_by_effect(
        AErange=(00, 50), cmap='viridis', whichp='ef', vmin=0, vmax=3, byb=0):
    # F16
    MLats, Mixps, EFs, AveEs = pd.read_pickle(savepath+'method5/F16.dat')
    fp = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1]) &
          (Mixps.eftot!=0)).values
    mlats16_1, mixps16_1, efs16_1, avees16_1 = (
            MLats[fp, :], Mixps.loc[fp, :], EFs[fp, :], AveEs[fp, :])
    # F17
    MLats, Mixps, EFs, AveEs = pd.read_pickle(savepath+'method5/F17.dat')
    fp = ((Mixps.AE>AErange[0]) & (Mixps.AE<AErange[1]) &
          (Mixps.eftot!=0)).values
    mlats17_1, mixps17_1, efs17_1, avees17_1 = (
            MLats[fp, :], Mixps.loc[fp, :], EFs[fp, :], AveEs[fp, :])
    # combine
    mlats1 = np.concatenate((mlats16_1, mlats17_1), axis=0)
    mixps1 = pd.concat((mixps16_1, mixps17_1), axis=0, ignore_index=True)
    efs1 =np.concatenate((efs16_1, efs17_1), axis=0)
    avees1 =np.concatenate((avees16_1, avees17_1), axis=0)

    plt.close('all')
    fig, ax = plt.subplots(2, 2, subplot_kw={'polar':True}, figsize=[6.5, 6])
    for ins, ns in enumerate(['N', 'S']):
        fp1 = (mixps1.ns==ns)
        for iby, by in enumerate(['positive', 'negative']):
            plt.sca(ax[ins, iby])
            if by=='positive':
                fp2 = (mixps1.Bym>byb)
            else:
                fp2 = (mixps1.Bym<-byb)
            fp = ((fp1) & (fp2)).values
            mlats2 = mlats1[fp, :]
            mixps2 = mixps1.loc[fp, :]
            efs2 = efs1[fp, :]
            avees2 = avees1[fp, :]
            for imlt, mlt in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
                if (mlt>10) & (mlt<14):
                    continue
                fp = (mixps2.mlt==mlt).values
                if np.sum(fp)< LLSMP:
                    continue
                mlatsave = np.median(mlats2[fp, :], axis=0)
                efsave = np.median(efs2[fp, :], axis=0)
                aveesave = np.median(avees2[fp, :], axis=0)
                mltave = np.ones(efsave.shape)*mlt
                if whichp == 'ef':
                    hc = plt.scatter(mltave/12*np.pi, 90-mlatsave, c=efsave,
                                     cmap=cmap, vmin=vmin, vmax=vmax, s=6)
                else:
                    hc = plt.scatter(mltave/12*np.pi, 90-mlatsave, c=aveesave,
                                     cmap=cmap, vmin=vmin, vmax=vmax, s=6)
            if by=='positive':
                plt.title('{:s}, By > {:d} nt'.format(ns, byb))
            else:
                plt.title('{:s}, By < -{:d} nt'.format(ns, byb))
            mf.set_polar(plt.gca(), 'N', boundinglat=60)
    cax = plt.axes([0.35, 0.06, 0.3, 0.03])
    plt.colorbar(hc, cax=cax, orientation='horizontal')
    if whichp == 'ef':
        plt.xlabel('Energy flux (ergs/s/cm$^2$)')
    else:
        plt.xlabel('Mean Energy (keV)')
    plt.show()
    return


def hp_statistics():
    # F16
    MLats16, Mixps16, EFs16, AveEs16 = pd.read_pickle(
            savepath+'method5/F16.dat')
    # F17
    MLats17, Mixps17, EFs17, AveEs17 = pd.read_pickle(
            savepath+'method5/F17.dat')
    # combine
    MLats, Mixps, EFs, AveEs = (
            np.concatenate((MLats16, MLats17), axis=0),
            pd.concat((Mixps16, Mixps17), axis=0, ignore_index=True),
            np.concatenate((EFs16, EFs17), axis=0),
            np.concatenate((AveEs16, AveEs17), axis=0))

    hp = np.zeros(
            (AErange.shape[0], np.size(np.arange(BinMLT/2, 24, BinMLT)))
            )*np.nan
    for ii in np.arange(AErange.shape[0]):
        for imlt, mlt in enumerate(np.arange(BinMLT/2, 24, BinMLT)):
            if (mlt>6) & (mlt<18):
                continue
            fp = ((Mixps.AE>AErange[ii, 0]) & (Mixps.AE>AErange[ii, 0]) &
                  (Mixps.mlt==mlt)).values
            hp[ii, imlt] = (
                    (np.mean(Mixps.loc[fp, 'eftot'])*1e-3)*
                    (RE*PIXELSIZE_GEOMAGNETIC_LATITUDE/180*np.pi*
                     RE*PIXELSIZE_GEOMAGNETIC_LATITUDE/180*np.pi)*1e-9)
    plt.close('all')
    plt.plot(np.mean(AErange, axis=1), np.nansum(hp, axis=1), '-o')
    plt.plot(np.arange(1,400),
             1.3+0.048*np.arange(1,400)+0.241*np.sqrt(np.arange(1,400)))
    plt.xlim(0, 400)
    plt.ylim(0,30)
    plt.xlabel('AE')
    plt.ylabel('HP')
    plt.show()

if __name__ == '__main__':
    # plt.close('all')
    # sat = 'f16'
    # year = 2011
    # doy = 360
    # path = '/home/guod/WD4T/ssusi/ssusi.jhuapl.edu/data/{:s}/'\
    #        'apl/edr-aur/{:d}/{:03d}/'.format(sat, year, doy)
    # fns = os.listdir(path)
    # for fn in fns:
    #     find_parameters_one_file_5(path+fn, test=True)
    # plt.show()

    # find_parameters()

    # aurora_reconstruction_statistics_all('avee', cmap='nipy_spectral')

    # statistics()

    # fourier_fit(test=True)

    # fourier_fit_save_parameters()

    # aurora_reconstruction_fit_all('avee')

    # difference between statistics and fourier fit
    # plt.close('all')
    # whichp = 'avee' if True else 'ef'
    # fig, ax = plt.subplots(
    #         2, 4, subplot_kw={'polar': True}, figsize=(10.4, 6.5))
    # AErange = np.ones([8, 2])
    # AErange[:, 0] = np.arange(0, 400, 50)
    # AErange[:, 1] = AErange[:, 0]+50
    # for k0 in range(2):
    #     for k1 in range(4):
    #         axt = ax[k0, k1]
    #         axt, hc = diff_fit_statistics(
    #             axt, AErange=AErange[k0*4+k1, :], cmap='seismic',
    #             whichp=whichp, vmin=-10, vmax=10, s=3)
    #         axt.set_title('AE: [{:.0f}, {:.0f}]'.\
    #                 format(AErange[k0*4+k1, 0], AErange[k0*4+k1, 1]), y=1.06)
    #         mf.set_polar(axt, boundinglat=60)
    # cax = plt.axes([0.35, 0.06, 0.3, 0.03])
    # plt.colorbar(hc, cax=cax, orientation='horizontal')
    # plt.show()

    # hp()

    # linear_fit_AE()

    # imf_by_effect(cmap='nipy_spectral', vmax=3)

    hp_statistics()

    print('----------------------------------------------------------------')
    print('Remaining problems:')
    print('  1, Aurora with complicated structure')
    print('  2, Different total energy, different proportion?')
    print('  3, Remove dayglow?')
    print('  4, Use median or mean?')
    print('----------------------------------------------------------------')
