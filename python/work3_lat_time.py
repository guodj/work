#Global imports
import numpy as np
import matplotlib.pyplot as plt
import gitm
import glob
import pandas as pd
import calc_rusanov as cr


def plot_rho(read=True):
    if read:
        rho = []
        for fn1, fn2 in zip(fns1,fns2):
            g1, g2 = gitm.GitmBin(fn1), gitm.GitmBin(fn2)
            rho1, rho2 = g1['Rho'], g2['Rho']
            rhod = 100*(rho2 - rho1)/rho1
            rhot1 = []
            for alt in alts:
                ialt = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]-alt*1000))+2
                rhot2 = []
                for lt in lts:
                    ilt = np.argmin(np.abs(g1['LT'][2:-2, 0, 0]-lt))+2
                    rhot2.append(
                        np.array(rhod[ilt, 2:-2, ialt]).reshape(1,1,1,-1))
                rhot2 = np.concatenate(rhot2, axis=2)
                rhot1.append(rhot2)
            rhot1 = np.concatenate(rhot1, axis=1)
            rho.append(rhot1)
        rho = np.concatenate(rho, axis=0)
        lat = g1['dLat'][0, 2:-2, 0]
        np.savez('/home/guod/tmp/rho', rho=rho, lat=lat)
    ldf = np.load('/home/guod/tmp/rho.npz')
    rho, lat = ldf['rho'], ldf['lat']
    plt.close()
    plt.contourf(
        timeidx, lat, rho[:,3, 1, :].T, levels=np.linspace(-5,5,21),
        cmap='seismic')
    plt.ylim(-90, -40)
    plt.show()
    return


def plot_ewind(read=True):
    if read:
        ewind= []
        for fn1, fn2 in zip(fns1,fns2):
            g1, g2 = gitm.GitmBin(fn1), gitm.GitmBin(fn2)
            ewind1, ewind2 = g1['V!Dn!N (east)'], g2['V!Dn!N (east)']
            ewindd = ewind2 - ewind1
            ewindt1 = []
            for alt in alts:
                ialt = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]-alt*1000))+2
                ewindt2 = []
                for lt in lts:
                    ilt = np.argmin(np.abs(g1['LT'][2:-2, 0, 0]-lt))+2
                    ewindt2.append(
                        np.array(ewindd[ilt, 2:-2, ialt]).reshape(1,1,1,-1))
                ewindt2 = np.concatenate(ewindt2, axis=2)
                ewindt1.append(ewindt2)
            ewindt1 = np.concatenate(ewindt1, axis=1)
            ewind.append(ewindt1)
        ewind = np.concatenate(ewind, axis=0)
        lat = g1['dLat'][0, 2:-2, 0]
        np.savez('/home/guod/tmp/ewind', ewind=ewind, lat=lat)
    ldf = np.load('/home/guod/tmp/ewind.npz')
    ewind, lat = ldf['ewind'], ldf['lat']
    plt.close()
    plt.contourf(
        timeidx, lat, ewind[:,3, 1, :].T, levels=np.linspace(-100,100,21),
        cmap='seismic')
    plt.ylim(-90, -40)
    plt.show()
    return


def plot_divrhov_rho(read=True):
    if read:
        divrhov = []
        for fn1, fn2 in zip(fns1,fns2):
            g1, g2 = gitm.GitmBin(fn1), gitm.GitmBin(fn2)
            # density change (shrink)
            lon1 = np.array(g1['Longitude'])
            lat1 = np.array(g1['Latitude'])
            alt1 = np.array(g1['Altitude'])
            Re = 6371*1000 # Earth radius, unit: m
            RR = Re+alt1
            omega = 2*np.pi/(24*3600)
            rho1 = np.array(g1['Rho'])
            nwind1 = np.array(g1['V!Dn!N (north)'])
            ewind1 = np.array(g1['V!Dn!N (east)']) + omega*RR*np.cos(lat1)
            uwind1 = np.array(g1['V!Dn!N (up)'])
            div_rhov1 = \
                (cr.calc_div_hozt(lon1, lat1, alt1, rho1*nwind1, rho1*ewind1)\
                +cr.calc_div_vert(alt1, rho1*uwind1))/rho1

            # density change (no shrink)
            lon2 = np.array(g2['Longitude'])
            lat2 = np.array(g2['Latitude'])
            alt2 = np.array(g2['Altitude'])
            Re = 6371*1000 # Earth radius, unit: m
            RR = Re+alt2
            omega = 2*np.pi/(24*3600)
            rho2 = np.array(g2['Rho'])
            nwind2 = np.array(g2['V!Dn!N (north)'])
            ewind2 = np.array(g2['V!Dn!N (east)']) + omega*RR*np.cos(lat2)
            uwind2 = np.array(g2['V!Dn!N (up)'])
            div_rhov2 = \
                (cr.calc_div_hozt(lon2, lat2, alt2, rho2*nwind2, rho2*ewind2)\
                +cr.calc_div_vert(alt2, rho2*uwind2))/rho2
            div_rhov = div_rhov1 - div_rhov2

            divrhov1 = []
            for alt in alts:
                ialt = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]-alt*1000))+2
                divrhov2 = []
                for lt in lts:
                    ilt = np.argmin(np.abs(g1['LT'][2:-2, 0, 0]-lt))+2
                    divrhov2.append(
                        np.array(div_rhov[ilt, 2:-2, ialt]).reshape(1,1,1,-1))
                divrhov2 = np.concatenate(divrhov2, axis=2)
                divrhov1.append(divrhov2)
            divrhov1 = np.concatenate(divrhov1, axis=1)
            divrhov.append(divrhov1)
        divrhov = np.concatenate(divrhov, axis=0)
        lat = g1['dLat'][0, 2:-2, 0]
        np.savez('/home/guod/tmp/divrhov_rho', divrhov=divrhov, lat=lat)
    ldf = np.load('/home/guod/tmp/divrhov_rho.npz')
    divrhov, lat = ldf['divrhov'], ldf['lat']
    plt.close()
    plt.contourf(
        timeidx, lat, divrhov[:,3, 1, :].T, levels=np.linspace(-1,1,21)*1e-4,
        cmap='seismic')
    plt.ylim(-90, -40)
    plt.show()
    return


def plot_vert_divv(read=True):
    if read:
        divv = []
        for fn1, fn2 in zip(fns1,fns2):
            g1, g2 = gitm.GitmBin(fn1), gitm.GitmBin(fn2)
            velr = np.array(g1['V!Dn!N (up)'])
            div_v1 = cr.calc_div_vert(g1['Altitude'], velr)

            velr = np.array(g2['V!Dn!N (up)'])
            div_v2 = cr.calc_div_vert(g2['Altitude'], velr)
            div_v = div_v1 - div_v2

            divv1 = []
            for alt in alts:
                ialt = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]-alt*1000))+2
                divv2 = []
                for lt in lts:
                    ilt = np.argmin(np.abs(g1['LT'][2:-2, 0, 0]-lt))+2
                    divv2.append(
                        np.array(div_v[ilt, 2:-2, ialt]).reshape(1,1,1,-1))
                divv2 = np.concatenate(divv2, axis=2)
                divv1.append(divv2)
            divv1 = np.concatenate(divv1, axis=1)
            divv.append(divv1)
        divv = np.concatenate(divv, axis=0)
        lat = g1['dLat'][0, 2:-2, 0]
        np.savez('/home/guod/tmp/vert_divv', divv=divv, lat=lat)
    ldf = np.load('/home/guod/tmp/vert_divv.npz')
    divv, lat = ldf['divv'], ldf['lat']
    plt.close()
    plt.contourf(
        timeidx, lat, divv[:,3, 1, :].T, levels=np.linspace(-1,1,21)*1e-4,
        cmap='seismic')
    plt.ylim(-90, -40)
    plt.show()
    return


def plot_hozt_divv(read=True):
    if read:
        divv = []
        for fn1, fn2 in zip(fns1,fns2):
            g1, g2 = gitm.GitmBin(fn1), gitm.GitmBin(fn2)
            lat1 = np.array(g1['Latitude'])
            alt1 = np.array(g1['Altitude'])
            Re = 6371*1000 # Earth radius, unit: m
            RR = Re+alt1
            omega = 2*np.pi/(24*3600)
            nwind1 = np.array(g1['V!Dn!N (north)'])
            ewind1 = np.array(g1['V!Dn!N (east)']) + omega*RR*np.cos(lat1)
            div_v1 = cr.calc_div_hozt(
                g1['Longitude'], g1['Latitude'], g1['Altitude'],
                nwind1, ewind1)

            lat2 = np.array(g2['Latitude'])
            alt2 = np.array(g2['Altitude'])
            Re = 6371*1000 # Earth radius, unit: m
            RR = Re+alt2
            omega = 2*np.pi/(24*3600)
            nwind2 = np.array(g2['V!Dn!N (north)'])
            ewind2 = np.array(g2['V!Dn!N (east)']) + omega*RR*np.cos(lat2)
            div_v2 = cr.calc_div_hozt(
                g2['Longitude'], g2['Latitude'], g2['Altitude'],
                nwind2, ewind2)
            div_v = div_v1 - div_v2

            divv1 = []
            for alt in alts:
                ialt = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]-alt*1000))+2
                divv2 = []
                for lt in lts:
                    ilt = np.argmin(np.abs(g1['LT'][2:-2, 0, 0]-lt))+2
                    divv2.append(
                        np.array(div_v[ilt, 2:-2, ialt]).reshape(1,1,1,-1))
                divv2 = np.concatenate(divv2, axis=2)
                divv1.append(divv2)
            divv1 = np.concatenate(divv1, axis=1)
            divv.append(divv1)
        divv = np.concatenate(divv, axis=0)
        lat = g1['dLat'][0, 2:-2, 0]
        np.savez('/home/guod/tmp/hozt_divv', divv=divv, lat=lat)
    ldf = np.load('/home/guod/tmp/hozt_divv.npz')
    divv, lat = ldf['divv'], ldf['lat']
    plt.close()
    plt.contourf(
        timeidx, lat, divv[:,3, 1, :].T, levels=np.linspace(-1,1,21)*1e-4,
        cmap='seismic')
    plt.ylim(-90, -40)
    plt.show()
    return


def plot_vert_vgradrho_rho(read=True):
    if read:
        gradrho = []
        for fn1, fn2 in zip(fns1,fns2):
            g1, g2 = gitm.GitmBin(fn1), gitm.GitmBin(fn2)
            rho1 = np.array(g1['Rho'])
            vgradrho1 = \
                g1['V!Dn!N (up)'] \
                * cr.calc_rusanov_alts_ausm(g1['Altitude'],rho1)/g1['Rho']

            rho2 = np.array(g2['Rho'])
            vgradrho2 = \
                g2['V!Dn!N (up)']\
                * cr.calc_rusanov_alts_ausm(g2['Altitude'],rho2)/g2['Rho']
            vgradrhod = vgradrho1 - vgradrho2

            gradrho1 = []
            for alt in alts:
                ialt = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]-alt*1000))+2
                gradrho2 = []
                for lt in lts:
                    ilt = np.argmin(np.abs(g1['LT'][2:-2, 0, 0]-lt))+2
                    gradrho2.append(
                        np.array(vgradrhod[ilt, 2:-2, ialt]).reshape(1,1,1,-1))
                gradrho2 = np.concatenate(gradrho2, axis=2)
                gradrho1.append(gradrho2)
            gradrho1 = np.concatenate(gradrho1, axis=1)
            gradrho.append(gradrho1)
        gradrho = np.concatenate(gradrho, axis=0)
        lat = g1['dLat'][0, 2:-2, 0]
        np.savez('/home/guod/tmp/vert_vgradrho', gradrho=gradrho, lat=lat)
    ldf = np.load('/home/guod/tmp/vert_vgradrho.npz')
    gradrho, lat = ldf['gradrho'], ldf['lat']
    plt.close()
    plt.contourf(
        timeidx, lat, gradrho[:,3, 1, :].T, levels=np.linspace(-1,1,21)*1e-4,
        cmap='seismic')
    plt.ylim(-90, -40)
    plt.show()
    return


def plot_hozt_vgradrho_rho(read=True):
    if read:
        gradrho = []
        for fn1, fn2 in zip(fns1,fns2):
            g1, g2 = gitm.GitmBin(fn1), gitm.GitmBin(fn2)

            lon1 = np.array(g1['Longitude'])
            lat1 = np.array(g1['Latitude'])
            alt1 = np.array(g1['Altitude'])
            Re = 6371*1000 # Earth radius, unit: m
            RR = Re+alt1
            omega = 2*np.pi/(24*3600)
            rho1 = np.array(g1['Rho'])
            nwind1 = np.array(g1['V!Dn!N (north)'])
            ewind1 = np.array(g1['V!Dn!N (east)']) + omega*RR*np.cos(lat1)
            vgradrho1 = \
                (nwind1*cr.calc_rusanov_lats(lat1,alt1,rho1)\
                +ewind1*cr.calc_rusanov_lons(lon1,lat1,alt1,rho1))/rho1

            lon2 = np.array(g2['Longitude'])
            lat2 = np.array(g2['Latitude'])
            alt2 = np.array(g2['Altitude'])
            Re = 6371*1000 # Earth radius, unit: m
            RR = Re+alt2
            omega = 2*np.pi/(24*3600)
            rho2 = np.array(g2['Rho'])
            nwind2 = np.array(g2['V!Dn!N (north)'])
            ewind2 = np.array(g2['V!Dn!N (east)']) + omega*RR*np.cos(lat2)
            vgradrho2 = \
                (nwind2*cr.calc_rusanov_lats(lat2,alt2,rho2)\
                 +ewind2*cr.calc_rusanov_lons(lon2,lat2,alt2,rho2))/rho2
            vgradrhod = vgradrho1 - vgradrho2

            gradrho1 = []
            for alt in alts:
                ialt = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]-alt*1000))+2
                gradrho2 = []
                for lt in lts:
                    ilt = np.argmin(np.abs(g1['LT'][2:-2, 0, 0]-lt))+2
                    gradrho2.append(
                        np.array(vgradrhod[ilt, 2:-2, ialt]).reshape(1,1,1,-1))
                gradrho2 = np.concatenate(gradrho2, axis=2)
                gradrho1.append(gradrho2)
            gradrho1 = np.concatenate(gradrho1, axis=1)
            gradrho.append(gradrho1)
        gradrho = np.concatenate(gradrho, axis=0)
        lat = g1['dLat'][0, 2:-2, 0]
        np.savez('/home/guod/tmp/hozt_vgradrho', gradrho=gradrho, lat=lat)
    ldf = np.load('/home/guod/tmp/hozt_vgradrho.npz')
    gradrho, lat = ldf['gradrho'], ldf['lat']
    plt.close()
    plt.contourf(
        timeidx, lat, gradrho[:,3, 1, :].T, levels=np.linspace(-1,1,21)*1e-4,
        cmap='seismic')
    plt.ylim(-90, -40)
    plt.show()
    return


def plot_all():
    from matplotlib.ticker import MultipleLocator
    ldf = np.load('/home/guod/tmp/rho.npz')
    rho, lat = ldf['rho'], ldf['lat']
    ldf = np.load('/home/guod/tmp/ewind.npz')
    ewind, lat = ldf['ewind'], ldf['lat']
    ldf = np.load('/home/guod/tmp/divrhov_rho.npz')
    divrhov, lat = ldf['divrhov'], ldf['lat']
    ldf = np.load('/home/guod/tmp/vert_divv.npz')
    vdivv, lat = ldf['divv'], ldf['lat']
    ldf = np.load('/home/guod/tmp/hozt_divv.npz')
    hdivv, lat = ldf['divv'], ldf['lat']
    ldf = np.load('/home/guod/tmp/vert_vgradrho.npz')
    vgradrho, lat = ldf['gradrho'], ldf['lat']
    ldf = np.load('/home/guod/tmp/hozt_vgradrho.npz')
    hgradrho, lat = ldf['gradrho'], ldf['lat']
    plt.close()
    xtime = (timeidx - stime)/pd.Timedelta('1min')

    fig, ax = plt.subplots(7,1,sharex=True,sharey=True,figsize=[6.07, 9.35])
    plt.subplots_adjust(top=0.95, bottom=0.05)

    plt.sca(ax[0])
    plt.contourf(
        xtime, lat, rho[:,iialt, iilt, :].T,
        levels=np.linspace(-10,10,21), cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.set_label(r'$\rho$')
    plt.ylabel('Latitude')

    plt.sca(ax[1])
    plt.contourf(
        xtime, lat, ewind[:,iialt, iilt, :].T,
        levels=np.linspace(-200,200,21), cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.set_label(r'u (east)')
    plt.ylabel('Latitude')

    plt.sca(ax[2])
    plt.contourf(
        xtime, lat, divrhov[:,iialt, iilt, :].T,
        levels=np.linspace(-5,5,21)*1e-5, cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.formatter.set_powerlimits((0,0))
    hcb.update_ticks()
    hcb.set_label(r'$\frac{\nabla\cdot(\rho u)}{\rho}$')
    plt.ylabel('Latitude')

    plt.sca(ax[3])
    plt.contourf(
        xtime, lat, vdivv[:,iialt, iilt, :].T,
        levels=np.linspace(-5,5,21)*1e-5, cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.formatter.set_powerlimits((0,0))
    hcb.update_ticks()
    hcb.set_label(r'$\nabla\cdot w$')
    plt.ylabel('Latitude')

    plt.sca(ax[4])
    plt.contourf(
        xtime, lat, hdivv[:,iialt, iilt, :].T,
        levels=np.linspace(-5,5,21)*1e-5, cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.formatter.set_powerlimits((0,0))
    hcb.update_ticks()
    hcb.set_label(r'$\nabla\cdot u$')
    plt.ylabel('Latitude')

    plt.sca(ax[5])
    plt.contourf(
        xtime, lat, vgradrho[:,iialt, iilt, :].T,
        levels=np.linspace(-5,5,21)*1e-5, cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.formatter.set_powerlimits((0,0))
    hcb.update_ticks()
    hcb.set_label(r'$w\cdot\frac{\nabla\rho}{\rho}$')
    plt.ylabel('Latitude')

    plt.sca(ax[6])
    plt.contourf(
        xtime, lat, hgradrho[:,iialt, iilt, :].T,
        levels=np.linspace(-5,5,21)*1e-5, cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.formatter.set_powerlimits((0,0))
    hcb.update_ticks()
    hcb.set_label(r'$u\cdot\frac{\nabla\rho}{\rho}$')
    plt.ylabel('Latitude')
    plt.xlabel('Time (minute)')

    ax[6].xaxis.set_major_locator(MultipleLocator(20))
    ax[6].xaxis.set_minor_locator(MultipleLocator(5))
    ax[6].yaxis.set_major_locator(MultipleLocator(10))
    ax[6].yaxis.set_minor_locator(MultipleLocator(5))


    plt.ylim(-90, -40)
    plt.show()

    return


if __name__=='__main__':
    import gc
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 01:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='1min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/'\
         + 'run_shrink_iondrift_4_c1_high_resolution/data/'
    fns1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fp2 = '/home/guod/simulation_output/momentum_analysis/'\
         + 'run_no_shrink_iondrift_4_1_high_resolution/data/'
    fns2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
            for k in timeidx]

    # save path
    spath = '/home/guod/Documents/work/fig/density_cell/' \
            'why_no_low_density_cell_at_high_latitude/time_lat/lt07/'

    alts = [130, 200, 300, 400, 500, 600]
    lts = [6, 7, 8, 9, 10, 11, 12]
    iialt = 2
    iilt = 1
    #plot_rho(read=False)
    #plot_ewind(read=False)
    #plot_divrhov_rho(read=False)
    #plot_vert_divv(read=False)
    #plot_hozt_divv(read=False)
    #plot_vert_vgradrho_rho(read=False)
    #plot_hozt_vgradrho_rho(read=False)
    plot_all()
    gc.collect()
