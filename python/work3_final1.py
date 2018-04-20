#Global imports
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from aacgmv2 import convert
import gitm_new as gitm
import gitm_3D_const_alt as g3ca
import gitm_create_coordinate as gcc
import cartopy.crs as ccrs
import calc_rusanov as cr
import pandas as pd
import glob
import datetime as dt
from mpl_toolkits.mplot3d import Axes3D

def figure1():
    plt.close('all')
    fig = plt.figure(figsize=(5.41, 8.54))
    ax = fig.add_subplot(111, projection='3d')

    levels = [np.linspace(6.1e-10, 17.8e-10,21), np.linspace(1.1e-10, 2.7e-10, 21),
              np.linspace(1.3e-11, 3.9e-11,21), np.linspace(2.6e-12, 9.503e-12,21),
              np.linspace(4.3e-13, 2.71e-12,21), np.linspace(7.4e-14, 8.2e-13,21)]
    olat = -40
    olatr = np.abs(olat)/180*pi

    glats, glons = convert(-90,0,0,date=dt.date(2003,3,22), a2g=True)
    glats, glons = glats/180*pi, glons/180*pi

    artalt = [10, 200, 380, 540, 680, 820]
    realalt = [150, 200, 300, 400, 500, 600]
    for ik, alt in enumerate([150, 200, 300, 400, 500, 600]):
        lon0, lat0, rho1 = g3ca.contour_data('Rho', g1, alt=alt)
        lat0, lon0 = lat0/180*pi, lon0/180*pi
        ilat = lat0[:, 0]<=olat/180*pi
        lat0, lon0, rho1 = lat0[ilat,:], lon0[ilat,:], rho1[ilat,:]
        r = lat0+pi/2
        theta = lon0
        x = r*np.cos(theta)
        x = -x # for sorthern hemisphere
        y = r*np.sin(theta)
        v = rho1
        print (np.min(v), np.max(v))
        ax.contourf(
            x, y, v, zdir='z', offset=artalt[ik], levels=levels[ik],
            cmap='jet', zorder=100)
        ax.scatter(
            -(glats+pi/2)*np.cos(glons), (glats+pi/2)*np.sin(glons),
            artalt[ik], s=50, c='k', zorder=200)
        for rtick in np.arange(10,50,10)/180*pi:
            thetatick = np.arange(0, 2*pi, 0.01)
            xtick = rtick*np.cos(thetatick)
            ytick = rtick*np.sin(thetatick)
            hp = ax.plot(
                xtick,ytick, artalt[ik], color='0.7',linestyle=':',
                zorder = 300, linewidth=0.8)
    plt.subplots_adjust(top=1,bottom=0)
    ax.set_xlim(-olatr, olatr)
    ax.set_ylim(-olatr, olatr)
    ax.set_zlim(0,800)
    ax.axis('off')

    for ik, alt in enumerate(artalt):
        ialt = np.argmin(np.abs(g1['Altitude'][0,0,:]/1000-realalt[ik]))
        ax.plot([-olatr, -olatr+0.1], [-olatr, -olatr], alt, color='k')
        ax.text(-olatr-0.1, -olatr, alt,
                '%.0f km'%(g1['Altitude'][0,0,ialt]/1000),
                verticalalignment='center')
    ax.plot([-olatr, -olatr], [-olatr, -olatr], [artalt[0], artalt[-1]], color='k')

    ax.text(-olatr-0.3, 0, artalt[-1], '06', horizontalalignment='center',
            zdir='x', zorder=1000)
    ax.text(0, -olatr-0.3, artalt[-1], '00', horizontalalignment='center',
            zdir='x', zorder=1000)
    ax.text(olatr+0.3, 0, artalt[-1], '18', horizontalalignment='center',
            zdir='x', zorder=1000)
    ax.text(0, olatr+0.3, artalt[-1], '12', horizontalalignment='center',
            zdir='x', zorder=1000)
    plt.show()
    return

def figure2_3(run=1):
    plt.close('all')
    plt.figure(figsize=(8.28, 9.25))
    if run==1:
        g = g1
    else:
        g = g2
    rholevels = [np.linspace(8, 16, 21)*1e-10, np.linspace(1,4,21)*1e-11,
                 np.linspace(0.5, 8,21)*1e-13]
    rhoticks = [np.arange(8, 17, 2)*1e-10, np.arange(1,5,1)*1e-11,
                np.arange(1, 9, 2)*1e-13]
    templevels = [np.linspace(600, 800, 21), np.linspace(900,1500,21),
                 np.linspace(900, 1500,21)]
    tempticks = [np.arange(600, 900, 100), np.arange(900,1600,200),
                np.arange(900, 1600, 200)]
    labels =[('(a)','(b)','(c)'), ('(d)','(e)','(f)'), ('(g)','(h)','(i)')]
    qlat, qlon = convert(-90, 0, 0, date=dt.date(2003,3,22), a2g=True)
    for ialt, alt in enumerate([150, 300, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, 2:-2]/1000-alt))+2
        alt_str = '%6.0f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
            3, 3, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
            coastlines=False, lonticklabel=[1,1,1,1])

        # density
        lon0, lat0, zdata0 = g3ca.contour_data('Rho', g, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(
            lon0, lat0, zdata0, rholevels[ialt], transform=ccrs.PlateCarree(),
            cmap='jet', extend='both')
        gcapos = ax.get_position()
        axcb = [
            gcapos.x0, gcapos.y0+0.028, gcapos.width, gcapos.height/20]
        axcb = plt.axes(axcb)
        hcb = plt.colorbar(
            hc, cax=axcb, extendrect=True, ticks=rhoticks[ialt],
            orientation='horizontal')

        # wind
        lon0, lat0, ewind, nwind = g3ca.vector_data(g, 'neutral', alt=alt)
        lon0, lat0, ewind, nwind = g3ca.convert_vector(
            lon0, lat0, ewind, nwind, plot_type='polar', projection=projection)
        hq = ax.quiver(
            lon0, lat0, ewind, nwind, scale=1000, scale_units='inches',
            regrid_shape=20)
        ax.quiverkey(hq, 0.97, 0.92, 500, '500 m/s')

        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        ax.scatter(
            0, -90, color='w', s=30, transform=ccrs.PlateCarree(), zorder=1000)

        ax.set_title('%s km' % alt_str, y=1.08)

        if ialt == 0:
            plt.text(
                -0.2, 0.5, r'$\rho$ and wind', rotation=90,
                transform=ax.transAxes, verticalalignment='center')
        plt.text(0, 0.9, labels[0][ialt],transform=ax.transAxes,color='r')

        #----------------------------------------------------------------------
        ax, projection = gcc.create_map(
            3, 3, ialt+4, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
            coastlines=False, lonticklabel=[1,1,1,1])
        # temperature
        lon0, lat0, zdata0 = g3ca.contour_data('Temperature', g, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(
            lon0, lat0, zdata0, templevels[ialt], transform=ccrs.PlateCarree(),
            cmap='jet', extend='both')
        gcapos = ax.get_position()
        axcb = [
            gcapos.x0, gcapos.y0-0.017, gcapos.width, gcapos.height/20]
        axcb = plt.axes(axcb)
        hcb = plt.colorbar(
            hc, cax=axcb, extendrect=True, ticks=tempticks[ialt],
            orientation='horizontal')

        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        ax.scatter(
            0, -90, color='w', s=30, transform=ccrs.PlateCarree(), zorder=1000)

        if ialt == 0:
            plt.text(
                -0.2, 0.5, 'Temperature (K)', rotation=90,
                transform=ax.transAxes, verticalalignment='center')
        plt.text(0, 0.9, labels[1][ialt],transform=ax.transAxes,color='r')
        #----------------------------------------------------------------------
        ax, projection = gcc.create_map(
            3, 3, ialt+7, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
            coastlines=False, lonticklabel=[1,1,1,1])
        # ion drift
        lon0, lat0, ewind, nwind = g3ca.vector_data(g, 'ion', alt=alt)
        lon0, lat0, ewind, nwind = g3ca.convert_vector(
            lon0, lat0, ewind, nwind, plot_type='polar', projection=projection)
        if run == 2:
            ewind, nwind = 0.1*ewind, 0.1*nwind
        hq = ax.quiver(
            lon0, lat0, ewind, nwind, scale=1000, scale_units='inches',
            regrid_shape=20)
        ax.quiverkey(hq, 0.97, 0.92, 500, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        ax.scatter(
            0, -90, color='w', s=30, transform=ccrs.PlateCarree(), zorder=1000)

        if ialt == 0:
            plt.text(
                -0.2, 0.5, 'Ion Drift', rotation=90,
                transform=ax.transAxes, verticalalignment='center')
        plt.text(0, 0.9, labels[2][ialt],transform=ax.transAxes,color='r')
    plt.subplots_adjust(wspace=0.25, hspace=0.25, top=0.95, bottom=0.05)

    if run==1:
        plt.savefig(savepath+'Figure2.eps')
        plt.savefig(savepath+'Figure2.jpeg')
    else:
        plt.savefig(savepath+'Figure3.eps')
        plt.savefig(savepath+'Figure3.jpeg')
    return


def figure4():
    plt.close('all')
    plt.figure(figsize=(8.28, 9.25))
    labels =[('(a)','(b)','(c)'), ('(d)','(e)','(f)'), ('(g)','(h)','(i)')]
    qlat, qlon = convert(-90, 0, 0, date=dt.date(2003,3,22), a2g=True)
    ylim = [[0.6e-9,1.6e-9], [1.5e-11, 4e-11], [1e-13, 9e-13]]
    for ialt, alt in enumerate([150, 300, 600]):
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]/1000-alt))+2
        alt_str = '%6.0f' % (g1['Altitude'][0, 0, alt_ind]/1000)
        rho1, rho2 = g1['Rho'], g2['Rho']
        rhod = 100*(rho1-rho2)/rho2
        lats, lons  = g1['Latitude'], g1['Longitude']

        #----------------------------------------------------------------------
        ax, projection = gcc.create_map(
            3, 3, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g1, 'polar',  useLT=True),
            coastlines=False, lonticklabel=[1,1,1,1])

        # density difference
        lon1, lat1, zdata1 = g3ca.contour_data('Rho', g1, alt=alt)
        lon2, lat2, zdata2 = g3ca.contour_data('Rho', g2, alt=alt)
        fp = (lat1[:,0]>slat) & (lat1[:,0]<nlat)
        lon0, lat0 = lon1[fp, :], lat1[fp,:]
        zdata0 = (100*(zdata1-zdata2)/zdata2)[fp, :]
        print(zdata0.min())
        hc = ax.contourf(
            lon0, lat0, zdata0, np.linspace(-30,30,21),
            transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        # colorbar
        gcapos = ax.get_position()
        axcb = [
            gcapos.x0, gcapos.y0+0.035, gcapos.width, gcapos.height/20]
        axcb = plt.axes(axcb)
        hcb = plt.colorbar(
            hc, cax=axcb, extendrect=True, ticks=np.arange(-30, 31, 15),
            orientation='horizontal')

        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind1-ewind2, nwind1-nwind2, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.97, 0.92, 500, '500 m/s')

        # magnetic pole
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        ax.scatter(
            0, -90, color='w', s=30, transform=ccrs.PlateCarree(), zorder=1000)

        latt1 = [-90, nlat]
        lont1 = [15*(12-ut), 15*(12-ut)]
        ax.plot(lont1, latt1, transform = ccrs.PlateCarree(), c='g')
        latt1 = [-90, nlat]
        lont1 = [15*(0-ut), 15*(0-ut)]
        ax.plot(lont1, latt1, transform = ccrs.PlateCarree(), c='g')
        ax.set_title('%s km' % alt_str, y=1.07)
        if ialt==0:
            ax.text(-0.2, 0.5, r'$\delta \rho$ (%)',
                    transform=ax.transAxes,rotation=90, fontsize=14,
                    verticalalignment='center', horizontalalignment='center')
        plt.text(
            0.05, 0.9, labels[0][ialt],transform=ax.transAxes,color='r')

        #----------------------------------------------------------------------
        # density in run1 and run2
        plt.subplot(3,3,ialt+7)
        ilat = (lats[0,:,0]/pi*180<nlat) &(lats[0,:,0]/pi*180>slat)
        lt = 12
        ilt = np.argmin(np.abs(g1['LT'][2:-2,0,0]-lt))+2
        plt.plot(
            180+lats[ilt,ilat,alt_ind]/pi*180,rho1[ilt,ilat,alt_ind], 'r',
            label = 'Run 1')
        plt.plot(
            180+lats[ilt,ilat,alt_ind]/pi*180,rho2[ilt,ilat,alt_ind], 'b',
            label = 'Run 2')
        plt.legend(loc='upper center')
        lat_1, rho_11, rho_12, rho_1d = \
            180+lats[ilt,2,alt_ind]/pi*180, rho1[ilt,2, alt_ind], \
            rho2[ilt,2, alt_ind], rhod[ilt,2,alt_ind]

        lt = 0
        ilt = np.argmin(np.abs(g1['LT'][2:-2,0,0]-lt))+2
        plt.plot(-lats[ilt,ilat,alt_ind]/pi*180,rho1[ilt,ilat,alt_ind], 'r')
        plt.plot(-lats[ilt,ilat,alt_ind]/pi*180,rho2[ilt,ilat,alt_ind], 'b')
        lat_2, rho_21, rho_22, rho_2d = \
            -lats[ilt,2,alt_ind]/pi*180, rho1[ilt,2, alt_ind], \
            rho2[ilt,2, alt_ind], rhod[ilt,2,alt_ind]

        plt.plot([lat_1,lat_2], [rho_11,rho_21], 'r')
        plt.plot([lat_1,lat_2], [rho_12,rho_22], 'b')
        plt.grid('on')

        plt.xlim([-nlat,180+nlat])
        plt.ylim(ylim[ialt])
        plt.xticks(np.arange(-nlat,181+nlat,30),
            np.concatenate(
                [np.arange(nlat,-90,-30), np.arange(-90,nlat+1,30)]))

        if ialt==0:
            plt.ylabel(r'$\rho$ (kg/m$^3$)', fontsize=14)
        plt.xlabel('Latitude')
        plt.text(
            0.05, 0.9, labels[2][ialt],transform=plt.gca().transAxes,color='r')

        #----------------------------------------------------------------------
        # density difference
        plt.subplot(3,3,ialt+4)
        lt = 12
        ilt = np.argmin(np.abs(g1['LT'][2:-2,0,0]-lt))+2
        plt.plot(
            180+lats[ilt,ilat,alt_ind]/pi*180,rhod[ilt,ilat,alt_ind], 'k')

        lt = 0
        ilt = np.argmin(np.abs(g1['LT'][2:-2,0,0]-lt))+2
        plt.plot(
            -lats[ilt,ilat,alt_ind]/pi*180,rhod[ilt,ilat,alt_ind], 'k')
        plt.plot([lat_1,lat_2], [rho_1d,rho_2d], 'k')

        plt.grid('on')

        plt.xlim([-nlat,180+nlat])
        plt.xticks(np.arange(-nlat,181+nlat,30),
            np.concatenate(
                [np.arange(nlat,-90,-30), np.arange(-90,nlat+1,30)]))
        plt.ylim(-35, 5)
        plt.yticks(np.arange(-30, 0.1, 10))
        if ialt==0:
            plt.ylabel(r'$\delta\rho$ (%)', fontsize=14)
        plt.xlabel('Latitude')
        plt.text(
            0.05, 0.9, labels[1][ialt],transform=plt.gca().transAxes,color='r')

    plt.subplots_adjust(wspace=0.25, hspace=0.25, top=0.95, bottom=0.05)
    plt.savefig(savepath+'Figure4.eps')
    plt.savefig(savepath+'Figure4.jpeg')
    return


def figure5_7(alt=150):
    times = pd.to_datetime([
        '2003-03-22 00:00:00', '2003-03-22 00:30:00', '2003-03-22 01:00:00',
        '2003-03-22 02:00:00', '2003-03-22 03:00:00', '2003-03-22 06:00:00'])
    plt.close('all')
    plt.figure(figsize=(9.5, 9.25))
    xtitle = [
        r'$-\nabla\cdot w$', r'$-\nabla\cdot u$',
        r'$-\frac{w\cdot\nabla\rho}{\rho}$',
        r'$-\frac{u\cdot\nabla\rho}{\rho}$',
        r'$\frac{\partial \rho}{\rho\partial t}$',
        r'$\delta\rho$ (%)']
    ylabel = times.strftime('%H:%M')
    qlat, qlon = convert(-90, 0, 0, date=dt.date(2003,3,22), a2g=True)
    for itime, time in enumerate(times):
        strtime = time.strftime('%y%m%d_%H%M')
        fn1 = glob.glob(
            '/home/guod/simulation_output/momentum_analysis/'\
            'run_no_shrink_iondrift_4_1/data/3DALL_t%s*.bin' % strtime)
        fn2 = glob.glob(
            '/home/guod/simulation_output/momentum_analysis/'\
            'run_shrink_iondrift_4_c1/data/3DALL_t%s*.bin' % strtime)
        g1 = gitm.read(fn1[0])
        g2 = gitm.read(fn2[0])

        lon1 = g1['Longitude']
        lat1 = g1['Latitude']
        alt1 = g1['Altitude']
        Re = 6371*1000 # Earth radius, unit: m
        RR = Re+alt1
        omega = 2*pi/(24*3600)
        rho1 = np.array(g1['Rho'])
        nwind1 = g1['V!Dn!N (north)']
        ewind1 = g1['V!Dn!N (east)'] + omega*RR*np.cos(lat1)
        uwind1 = g1['V!Dn!N (up)']
        div_rhov1 = \
            (cr.calc_div_hozt(lon1, lat1, alt1, rho1*nwind1, rho1*ewind1)\
            +cr.calc_div_vert(alt1, rho1*uwind1))/rho1
        vert_divv1 = cr.calc_div_vert(alt1, uwind1)
        hozt_divv1 = cr.calc_div_hozt(lon1, lat1, alt1, nwind1, ewind1)
        vert_vgradrho1 = uwind1*cr.calc_rusanov_alts_ausm(alt1,rho1)/rho1
        hozt_vgradrho1 = \
            (nwind1*cr.calc_rusanov_lats(lat1,alt1,rho1)\
            +ewind1*cr.calc_rusanov_lons(lon1,lat1,alt1,rho1))/rho1

        dc1 = [
            vert_divv1, hozt_divv1, vert_vgradrho1, hozt_vgradrho1,
            div_rhov1, rho1]

        lon2 = g2['Longitude']
        lat2 = g2['Latitude']
        alt2 = g2['Altitude']
        Re = 6371*1000 # Earth radius, unit: m
        RR = Re+alt2
        omega = 2*pi/(24*3600)
        rho2 = g2['Rho']
        nwind2 = g2['V!Dn!N (north)']
        ewind2 = g2['V!Dn!N (east)'] + omega*RR*np.cos(lat2)
        uwind2 = g2['V!Dn!N (up)']
        div_rhov2 = \
            (cr.calc_div_hozt(lon2, lat2, alt2, rho2*nwind2, rho2*ewind2)\
            +cr.calc_div_vert(alt2, rho2*uwind2))/rho2
        vert_divv2 = cr.calc_div_vert(alt2, uwind2)
        hozt_divv2 = cr.calc_div_hozt(lon2, lat2, alt2, nwind2, ewind2)
        vert_vgradrho2 = uwind2*cr.calc_rusanov_alts_ausm(alt2,rho2)/rho2
        hozt_vgradrho2 = \
            (nwind2*cr.calc_rusanov_lats(lat2,alt2,rho2)\
            +ewind2*cr.calc_rusanov_lons(lon2,lat2,alt2,rho2))/rho2

        dc2 = [
            vert_divv2, hozt_divv2, vert_vgradrho2, hozt_vgradrho2,
            div_rhov2, rho2]

        alt_ind = np.argmin(np.abs(alt1[0, 0, :]/1000-alt))
        alt_str = '%6.0f' % (alt1[0, 0, alt_ind]/1000)

        for idc in range(6):
            lonticklabel = [0, 0, 0, 0]
            if idc == 0:
                lonticklabel[3] = lonticklabel[-1]+1
            if idc == 5:
                lonticklabel[1] = lonticklabel[-1]+1
            if itime == 0:
                lonticklabel[0] = lonticklabel[-2]+1
            if itime == 5:
                lonticklabel[2] = lonticklabel[-2]+1
            ax, projection = gcc.create_map(
                6, 6, itime*6+idc+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1, 'polar',  useLT=True),
                coastlines=False, lonticklabel=lonticklabel)
            g1['temp'] = -(dc1[idc] - dc2[idc])
            if idc==5:
                g1['temp'] = 100*(dc1[idc] - dc2[idc])/dc2[idc]
            lon0, lat0, zdata0 = g3ca.contour_data('temp', g1, alt=alt)
            fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
            lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
            levels = np.linspace(-8,8,21)*1e-5
            if idc==5:
                levels = np.linspace(-30,30,21)
            hc = ax.contourf(
                lon0, lat0, zdata0, levels,
                transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
            # wind difference
            alt_ind = np.argmin(np.abs(alt1[0, 0, :]/1000-alt))
            alt_str = '%6.0f' % (alt1[0, 0, alt_ind]/1000)
            lon1, lat1, ewind1, nwind1 = g3ca.vector_data(
                g1, 'neutral', alt=alt)
            lon2, lat2, ewind2, nwind2 = g3ca.vector_data(
                g2, 'neutral', alt=alt)
            lon0, lat0 = lon1, lat1
            lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                    lon0, lat0, ewind1-ewind2, nwind1-nwind2,
                    plot_type='polar', projection=projection)
            hq = ax.quiver(
                    lon0, lat0, ewind0, nwind0, scale=1500,
                    scale_units='inches', regrid_shape=20, headwidth=5)
            ax.scatter(
                qlon, qlat, s=10, color='k', transform=ccrs.PlateCarree())
            ax.scatter(
                0, -90, color='w', s=10, transform=ccrs.PlateCarree(),
                zorder=1000)
            if idc==0:
                plt.text(
                    -0.3, 0.5, ylabel[itime], rotation=90,
                    transform=ax.transAxes, verticalalignment='center')
            if itime==0:
                plt.title(xtitle[idc], y=1.15, fontsize=14)
                texts = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
                plt.text(0.05,1,texts[idc],transform=plt.gca().transAxes,
                         color='r')
            if itime==5:
                axp = ax.get_position()
                cax = [axp.x0, axp.y0-0.03, axp.x1-axp.x0, 1/20*(axp.y1-axp.y0)]
                cax = plt.axes(cax)
                ticks = np.arange(-8e-5, 8.1e-5, 4e-5) if idc<5 else \
                        np.arange(-30, 31, 15)
                hcb = plt.colorbar(
                    hc, cax=cax, extendrect=True, orientation='horizontal',
                    ticks=ticks)
                if idc<5:
                    hcb.formatter.set_powerlimits((0,0))
                    hcb.update_ticks()
    plt.subplots_adjust(wspace=0, hspace=0)
    if alt==300:
        plt.savefig(savepath+'Figure5.eps')
        plt.savefig(savepath+'Figure5.jpeg')
    if alt==600:
        plt.savefig(savepath+'Figure7.eps')
        plt.savefig(savepath+'Figure7.jpeg')
    return

def figure5_7_refer(alt=300):
    times = pd.to_datetime([
        '2003-03-22 00:00:00', '2003-03-22 00:30:00', '2003-03-22 01:00:00',
        '2003-03-22 02:00:00', '2003-03-22 03:00:00', '2003-03-22 06:00:00'])
    plt.close('all')
    plt.figure(figsize=(3, 9.25))
    ylabel = times.strftime('%H:%M')
    qlat, qlon = convert(-90, 0, 0, date=dt.date(2003,3,22), a2g=True)
    for itime, time in enumerate(times):
        strtime = time.strftime('%y%m%d_%H%M')
        fn1 = glob.glob(
            '/home/guod/simulation_output/momentum_analysis/'\
            'run_no_shrink_iondrift_4_1/data/3DALL_t%s*.bin' % strtime)
        fn2 = glob.glob(
            '/home/guod/simulation_output/momentum_analysis/'\
            'run_shrink_iondrift_4_c1/data/3DALL_t%s*.bin' % strtime)
        g1 = gitm.read(fn1[0])
        g2 = gitm.read(fn2[0])

        lon1, lat1, alt1, T1 = \
            g1['Longitude'], g1['Latitude'], g1['Altitude'], g1['Temperature']
        lon2, lat2, alt2, T2 = \
            g2['Longitude'], g2['Latitude'], g2['Altitude'], g2['Temperature']

        alt_ind = np.argmin(np.abs(alt1[0, 0, :]/1000-alt))
        alt_str = '%6.0f' % (alt1[0, 0, alt_ind]/1000)

        lonticklabel=[0, 1, 0, 1]
        if itime == 0:
            lonticklabel[0] = 1
        if itime == 5:
            lonticklabel[2] = 1
        ax, projection = gcc.create_map(
            6, 1, itime+1, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g1, 'polar',  useLT=True),
            coastlines=False, lonticklabel=lonticklabel)
        g1['temp'] = T1-T2
        lon0, lat0, zdata0 = g3ca.contour_data('temp', g1, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(
            lon0, lat0, zdata0, levels=np.linspace(-60, 60, 21),
            transform=ccrs.PlateCarree(), cmap='seismic', extend='both')

        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(
            g1, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(
            g2, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind1-ewind2, nwind1-nwind2,
                plot_type='polar', projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1500,
                scale_units='inches', regrid_shape=20, headwidth=5)
        ax.scatter(
            qlon, qlat, s=10, color='k', transform=ccrs.PlateCarree())
        ax.scatter(
            0, -90, color='w', s=10, transform=ccrs.PlateCarree(),
            zorder=1000)
        plt.text(
            -0.3, 0.5, ylabel[itime], rotation=90,
            transform=ax.transAxes, verticalalignment='center')
        if itime==0:
            plt.title(r'$T_d$ (K)', y=1.15)
        if itime==5:
            axp = ax.get_position()
            cax = [axp.x0, axp.y0-0.03, axp.x1-axp.x0, 1/20*(axp.y1-axp.y0)]
            cax = plt.axes(cax)
            ticks = np.arange(-60, 61, 20)
            hcb = plt.colorbar(
                hc, cax=cax, extendrect=True, orientation='horizontal',
                ticks=ticks)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()
    return

def figure6_8(iialt=2):
    from matplotlib.ticker import MultipleLocator
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 01:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='1min')
    alts = [150, 200, 300, 400, 500, 600]
    lts = [6, 7, 8, 9, 10, 11, 12]
    iilt = 1
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

    fig, ax = plt.subplots(6,1,sharex=True,sharey=True,figsize=[6.07, 9.35])
    plt.subplots_adjust(top=0.95, bottom=0.05)

    #--------------------------------------------------------------------------
    plt.sca(ax[0])
    levels = np.linspace(-20, 20, 41)
    plt.contourf(
        xtime, lat, rho[:,iialt, iilt, :].T,
        levels=levels, cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.set_ticks(levels[::10])
    hcb.set_label(r'$\delta\rho$ (%)', fontsize=14)
    plt.ylabel('Latitude')
    plt.text(0.05,0.8,'(a)',transform=plt.gca().transAxes, color='r')
    galt = g1['Altitude'][0,0,:]
    alt_ind = np.argmin(np.abs(galt/1000-alts[iialt]))
    plt.title('LT = %02d, %5.0f km' % (lts[iilt], galt[alt_ind]/1000))
    print(galt[alt_ind])

    #--------------------------------------------------------------------------
    plt.sca(ax[1])
    levels = np.linspace(-8,8,41)*1e-5
    plt.contourf(
        xtime, lat, divrhov[:,iialt, iilt, :].T,
        levels=levels, cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.set_ticks(levels[::10])
    hcb.formatter.set_powerlimits((0,0))
    hcb.update_ticks()
    hcb.set_label(r'$\frac{\partial \rho}{\rho\partial t}$', fontsize=14)
    plt.ylabel('Latitude')
    plt.text(0.05,0.8,'(b)',transform=plt.gca().transAxes, color='r')

    #--------------------------------------------------------------------------
    plt.sca(ax[2])
    levels = np.linspace(-8,8,41)*1e-5
    plt.contourf(
        xtime, lat, vdivv[:,iialt, iilt, :].T,
        levels=levels, cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.set_ticks(levels[::10])
    hcb.formatter.set_powerlimits((0,0))
    hcb.update_ticks()
    hcb.set_label(r'$-\nabla\cdot w$', fontsize=14)
    plt.ylabel('Latitude')
    plt.text(0.05,0.8,'(c)',transform=plt.gca().transAxes, color='r')

    #--------------------------------------------------------------------------
    plt.sca(ax[3])
    levels = np.linspace(-8,8,41)*1e-5
    plt.contourf(
        xtime, lat, hdivv[:,iialt, iilt, :].T,
        levels=levels, cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.set_ticks(levels[::10])
    hcb.formatter.set_powerlimits((0,0))
    hcb.update_ticks()
    hcb.set_label(r'$-\nabla\cdot u$', fontsize=14)
    plt.ylabel('Latitude')
    plt.text(0.05,0.8,'(d)',transform=plt.gca().transAxes, color='r')

    #--------------------------------------------------------------------------
    plt.sca(ax[4])
    levels = np.linspace(-8,8,41)*1e-5
    plt.contourf(
        xtime, lat, vgradrho[:,iialt, iilt, :].T,
        levels=levels, cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.set_ticks(levels[::10])
    hcb.formatter.set_powerlimits((0,0))
    hcb.update_ticks()
    hcb.set_label(r'$-w\cdot\frac{\nabla\rho}{\rho}$', fontsize=14)
    plt.ylabel('Latitude')
    plt.text(0.05,0.8,'(e)',transform=plt.gca().transAxes, color='r')

    #--------------------------------------------------------------------------
    plt.sca(ax[5])
    levels = np.linspace(-8,8,41)*1e-5
    plt.contourf(
        xtime, lat, hgradrho[:,iialt, iilt, :].T,
        levels=levels, cmap='seismic', extend='both')
    hcb = plt.colorbar(extendrect=True, pad=0.05, aspect=10)
    hcb.set_ticks(levels[::10])
    hcb.formatter.set_powerlimits((0,0))
    hcb.update_ticks()
    hcb.set_label(r'$-u\cdot\frac{\nabla\rho}{\rho}$', fontsize=14)
    plt.ylabel('Latitude')
    plt.xlabel('Time (minute)')
    plt.text(0.05,0.8,'(f)',transform=plt.gca().transAxes, color='r')

    ax[-1].xaxis.set_major_locator(MultipleLocator(10))
    ax[-1].xaxis.set_minor_locator(MultipleLocator(5))
    ax[-1].yaxis.set_major_locator(MultipleLocator(30))
    ax[-1].yaxis.set_minor_locator(MultipleLocator(10))


    plt.ylim(-90, -30)
    if iialt==2:
        plt.savefig(savepath+'Figure6.eps')
        plt.savefig(savepath+'Figure6.jpeg')
    if iialt==5:
        plt.savefig(savepath+'Figure8.eps')
        plt.savefig(savepath+'Figure8.jpeg')

    return


if __name__ == '__main__':
    fn1 = '/home/guod/simulation_output/momentum_analysis/'\
          'run_no_shrink_iondrift_4_1/data/3DALL_t030322_060000.bin'
    fn2 = '/home/guod/simulation_output/momentum_analysis/'\
          'run_shrink_iondrift_4_c1/data/3DALL_t030322_060000.bin'
    # fn1 = '/home/guod/simulation_output/momentum_analysis/'\
    #       'run_no_shrink_no_tides_1_1/UA/data/3DALL_t030322_060001.bin'
    # fn2 = '/home/guod/simulation_output/momentum_analysis/'\
    #       'run_shrink_no_tides_1_c1/UA/data/3DALL_t030322_060000.bin'
    g1 = gitm.read(fn1)
    g2 = gitm.read(fn2)
    ut = g1['time'].hour +g1['time'].minute/60
    nlat, slat= -30, -90
    savepath = '/home/guod/Documents/Work/DensityCell/'\
               'NoIsolatedCellHighAltitude/Paper/'
    # figure1()
    # figure2_3(run=1)
    # figure2_3(run=2)
    # figure4()
    figure5_7(alt=300)
    figure5_7(alt=600)
    # figure5_7_refer(alt=300)
    # figure5_7_refer(alt=600)
    # figure6_8(iialt=2)
    # figure6_8(iialt=5)
