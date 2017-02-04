#!/home/guod/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt

def func1():
    x, y = np.arange(4), [0, 1, 1, 1]
    plt.plot(x, y)
    plt.xlim(0, 3)
    plt.ylim(0, 1.5)
    plt.xticks([])
    plt.yticks([])
    plt.xlabel('Time')
    plt.ylabel('Density difference')
    plt.show()


def func2():
    def f1():
        plt.plot(np.arange(10))
    def f2():
        plt.plot(np.arange(10)[::-1])
    f1()
    f2()
    plt.show()


def gitm_lat_lt_nspolar():
    import gitm
    import gitm_3D_global_plots as gpt
    import myfunctions as mf
    from scipy.interpolate import griddata
    from scipy import interpolate

    def gitm_lat_lt_grid(gitmdata, height):
        """
        Derive rho in Latitude-LT grid
        input:
            gitmdata: GitmBin
            height: which height (km) do you want?
            centerlat: 90 or -90
            boundinglat: latitude at the edge
        output:
            theta0: localtime(0-24) -> radian(0-360)
            r0: latitude(0-90)
            rho0: gitm density
        """
        h = gitmdata['Altitude'][0,0,:]
        ih = np.argmin(abs(h-height*1000))
        lat_data = np.array(gitmdata['Latitude'][:, :, ih])/np.pi*180 # degree
        lt_data = np.array(gitmdata['LT'][:, :, ih]) # hour
        rho_data = np.array(gitmdata['Rho'][:, :, ih])

        # Set index of 0 LT at 0
        ind = np.argmin(lt_data[:,0])
        lat0 = np.roll(lat_data, -ind, axis=0)
        lt0 = np.roll(lt_data, -ind, axis=0)
        rho0 = np.roll(rho_data, -ind, axis=0)
        return lat0, lt0, rho0

    def lt2theta(lat, lt, rho, centerlat=90, boundinglat=0):
        # Exclude latitudes greater than 90 or smaller than -90
        fp = (lat[0, :]>=-90) &(lat[0, :]<90)
        lat0 = lat[:, fp]
        lt0 = lt[:, fp]
        rho0 = rho[:, fp]

        csign = np.sign(centerlat)
        r0 = 90 - csign*lat0
        theta0 = lt0/12*np.pi #radian

        fp = r0[0, :]<=90
        r0 = r0[:, fp]
        theta0 = theta0[:, fp]
        rho0 = rho0[:, fp]
        return theta0, r0, rho0

    def gitm_mlat_mlt_grid(gitmdata, height, centerlat=90, boundinglat=0):
        """
        Derive rho in Latitude-LT grid
        input:
            gitmdata: GitmBin
            height: which height (km) do you want?
            centerlat: 90 or -90
            boundinglat: latitude at the edge
        """
        from apexpy import Apex
        h = gitmdata['Altitude'][0,0,:]
        ih = np.argmin(abs(h-height*1000))
        lat_data = np.array(gitmdata['Latitude'][:, :, ih])/np.pi*180 # degree
        lon_data = np.array(gitmdata['Longitude'][:, :, ih])/np.pi*180
        lt_data = np.array(gitmdata['LT'][:, :, ih]) # hour
        rho_data = np.array(gitmdata['Rho'][:, :, ih])

        lat0 = lat_data.reshape(-1)
        lt0 = lt_data.reshape(-1)
        rho0 = rho_data.reshape(-1)
        lon0 = lon_data.reshape(-1)
        fp = (lat0<=90) & (lat0>=-90)
        lat0, lt0, rho0, lon0 = lat0[fp], lt0[fp], rho0[fp], lon0[fp]
        mlat0, mlt0 = Apex().convert(lat=lat0, lon=lon0, source='geo',
                                     dest='mlt', datetime=gitmdata['time'],
                                     height=h[ih]/1000)
        #mlat0, tmp1 = Apex().convert(lat=lat0, lon=lon0, source='geo',
        #                             dest='qd', datetime=gitmdata['time'],
        #                             height=h[ih]/1000)

        csign = np.sign(centerlat)
        mlat0 = csign*mlat0
        fp = (mlat0>=csign*boundinglat) & (mlat0<=90)
        mlat0, mlt0, rho0 = mlat0[fp], mlt0[fp], rho0[fp]

        r0 = 90 - mlat0
        theta0 = mlt0/12*np.pi-np.pi/2 # -pi/2 to make 0 LT at bottom

        x0 = r0*np.cos(theta0)
        y0 = r0*np.sin(theta0)

        edger = 90-csign*boundinglat
        x1, y1 = np.meshgrid(np.arange(-edger, edger, 0.5),
                             np.arange(-edger, edger, 0.5))
        z1 = griddata((x0, y0), rho0, (x1, y1), method='linear')
        return x1, y1, z1

    ff = '/home/guod/tmp/3DALL_t1002'
    fname = ['26_210000.bin', '26_223000.bin', '27_000000.bin',
             '27_013000.bin', '27_030000.bin']
    gitm0 = gitm.GitmBin(ff+'27_000000.bin',
                         varlist=['Latitude', 'Longitude', 'LT', 'Rho'])
    #    theta0, r0, rho0 = gitm_lat_lt_grid(gitm0, 400, centerlat=90)
    #    ax = plt.subplot(111, polar=True)
    #    plt.contourf(theta0, r0, rho0, 20)
    fig, ax = plt.subplots(5, 2, figsize=(4.6,11.5), subplot_kw=dict(polar=True))
    for k11 in range(5):
        gitm1 = gitm.GitmBin(ff+fname[k11],
                             varlist=['Latitude', 'Longitude', 'LT', 'Rho'])
        for k00, k0 in enumerate([90, -90]):
            plt.sca(ax[k11, k00])
            #mf.add_polar_coordinate(plt.gca(), centerlat=k0,
            #                        rxl=0.25, rtl=0.1,
            #                        lonlabel='LT')
            x0, y0, z0 = gitm_lat_lt_grid(gitm0, 400)
            x1, y1, z1 = gitm_lat_lt_grid(gitm1, 400)
            x2, y2, z2 = lt2theta(x0, y0, 100*(z1-z0)/z0, centerlat=k0)
            #plt.axis('equal')
            #plt.axis('off')
            conh = plt.contourf(x2, y2, z2,
                                np.linspace(-25, 25, 11),
                                cmap='seismic_r', extend='both')
            #plt.xlim(-92,92)
            #plt.ylim(-92,92)
    plt.subplots_adjust(left=0.03, right=0.97, top=0.95, bottom=0.10)
    cax = plt.axes([0.2, 0.05, 0.6, 0.01])
    cbh = plt.colorbar(conh, cax=cax, orientation='horizontal')
    cbh.set_ticks(np.arange(-25, 26, 5))
    return gitm0


def func4():
    """
    Test (lat, lt) to (mlat, mlt) conversions.

    1, Use a = Apex(date=...); "date" determines which IGRF coefficients are
       used in conversions. Uses current date as default.
    2, height is needed for better conversion.
       champ and grace files use qd coordinates.
    3, mlt in qd and apex coordinates are the same.
    """
    import champ_grace as cg
    from apexpy import Apex as Apex
    import matplotlib.pyplot as plt
    a = cg.ChampDensity('2005-1-1', '2005-1-2')
    b = Apex(date=2005)
    mlatt, mlt = b.convert(
            lat=a.lat, lon=a.long, source='geo', dest='mlt',
            datetime=a.index, height=a.height)
    mlat, mlongt = b.convert(
            lat=a.lat, lon=a.long, source='geo', dest='qd',
            height=a.height)
    mlat2 = np.array(a.Mlat)
    mlt2 = np.array(a.MLT)
    plt.plot(mlat-mlat2)
    plt.plot(abs(mlt-mlt2) % 24, '.')
    plt.show()
    return


if __name__ == '__main__':
    plt.close('all')
    import gc
    gc.collect()
    a = func4()
    plt.show()
