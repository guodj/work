#Global imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from apexpy import Apex
import gitm
import gitm_3D_const_alt as g3ca
import gitm_create_coordinate as gcc
import cartopy.crs as ccrs
from pylab import draw # Use draw()
from spacepy.datamodel import dmarray
import gitm_pressure as gp
from cartopy.util import add_cyclic_point
import matplotlib.animation as animation
import glob
import pandas as pd
import gitm_divergence as gd
import calc_rusanov as cr
sns.set('paper', 'whitegrid')

def plot_animation_den_win(show=False, f107=150, which_alt=600):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_3_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_3/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        if which_alt==400:
            rholevels=np.linspace(0.6, 1.4, 21)*1e-11
        if which_alt==600:
            rholevels=np.linspace(0.32, 1.8, 21)*1e-12
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    fig = plt.figure(figsize=[8,9])
    # read gitm data
    def animate_den_wind(i):
        g1, g2 = [gitm.GitmBin(k[i], varlist=[
            'Rho', 'V!Dn!N (north)', 'V!Dn!N (east)', 'V!Dn!N (up)',
            'V!Di!N (east)', 'V!Di!N (north)']) for k in [fn1, fn2]]
        # create axis
        ax = list(range(6))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            for irun in range(3):
                ax[ins+irun*2], projection[ins+irun*2] = gcc.create_map(
                        3, 2, 1+ins+irun*2, 'polar', nlat=nlat, slat=slat,
                        dlat=10, centrallon=g3ca.calculate_centrallon(
                            g1, 'polar',  useLT=True),
                        coastlines=False)
        # Density
        lon1, lat1, zdata1 = g3ca.contour_data('Rho', g1, alt=which_alt)
        lon2, lat2, zdata2 = g3ca.contour_data('Rho', g2, alt=which_alt)
        hc = [ax[k].contourf(
                lon1, lat1, zdata1, transform=ccrs.PlateCarree(),
                levels=rholevels,
                cmap='jet', extend='both') for k in [0, 1]]
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = [ax[k].contourf(
                lon2, lat2, zdata2, transform=ccrs.PlateCarree(),
                levels=rholevels,
                cmap='jet', extend='both') for k in [2, 3]]
        diffzdata = 100*(zdata2-zdata1)/zdata1
        hc = [ax[k].contourf(
                lon2, lat2, diffzdata, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-30, 30, 21), cmap='seismic',
                extend='both') for k in [4, 5]]
        hc = [ax[k].contour(
                lon2, lat2, diffzdata, [-10], transform=ccrs.PlateCarree(),
                colors='g',linestyles='-') for k in [4, 5]]

        # wind
        lon1, lat1, ewind1, nwind1 = \
                g3ca.vector_data(g1, 'neutral', alt=which_alt)
        lon2, lat2, ewind2, nwind2 = \
                g3ca.vector_data(g2, 'neutral', alt=which_alt)
        for iax in range(6):
            if iax == 0 or iax == 1:
                lon0, lat0, ewind, nwind = (
                        lon1.copy(), lat1.copy(), ewind1.copy(), nwind1.copy())
                lon0, lat0, ewind, nwind = g3ca.convert_vector(
                        lon0, lat0, ewind, nwind, plot_type='polar',
                        projection=projection[iax])
            elif iax == 2 or iax == 3:
                lon0, lat0, ewind, nwind = (
                        lon2.copy(), lat2.copy(), ewind2.copy(), nwind2.copy())
                lon0, lat0, ewind, nwind = g3ca.convert_vector(
                        lon0, lat0, ewind, nwind, plot_type='polar',
                        projection=projection[iax])
            elif iax == 4 or iax == 5:
                lon0, lat0, ewind, nwind = (
                        lon1.copy(), lat1.copy(), ewind2-ewind1, nwind2-nwind1)
                lon0, lat0, ewind, nwind = g3ca.convert_vector(
                        lon0, lat0, ewind, nwind, plot_type='polar',
                        projection=projection[iax])
            hq = ax[iax].quiver(
                    lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
                    color='k', regrid_shape=20)
            # ax.quiverkey(hq, 0.93, 0, 1000, '1000 m/s')
            # hc = plt.colorbar(hc, ticks=np.arange(3, 7)*1e-12)
            # hc.set_label(r'$\rho$ (kg/m$^3$)')
            # ax[iax].scatter(
            #         qlonn, qlatn, color='k', transform=ccrs.PlateCarree())
            # ax[iax].scatter(
            #         qlons, qlats, color='k', transform=ccrs.PlateCarree())
        return
    anim = animation.FuncAnimation(
            fig, animate_den_wind, frames=len(fn1))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, bitrate=1800)
    anim.save(path+'01_den_wind_%d_%d.wmv' % (f107, which_alt),writer=writer)
    return


def plot_animation_vert_wind(show=False, f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_3_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_3/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    fig = plt.figure(figsize=[8,9])
    # read gitm data
    def animate_vert_wind(i):
        g1, g2 = [gitm.GitmBin(k[i]) for k in [fn1, fn2]]
        # create axis
        ax = list(range(6))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            for irun in range(3):
                ax[ins+irun*2], projection[ins+irun*2] = gcc.create_map(
                        3, 2, 1+ins+irun*2, 'polar', nlat=nlat, slat=slat,
                        dlat=10, centrallon=g3ca.calculate_centrallon(
                            g1, 'polar',  useLT=True),
                        coastlines=False)
        # vertical wind
        lon1, lat1, zdata1 = g3ca.contour_data('V!Dn!N (up)', g1, alt=which_alt)
        lon2, lat2, zdata2 = g3ca.contour_data('V!Dn!N (up)', g2, alt=which_alt)
        lon3, lat3, zdata3 = g3ca.contour_data('Rho', g1, alt=which_alt)
        lon4, lat4, zdata4 = g3ca.contour_data('Rho', g2, alt=which_alt)
        hc = [ax[k].contourf(
                lon1, lat1, zdata1, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-20, 20, 21),
                #levels=np.linspace(15e-13, 25e-13, 21),
                cmap='seismic', extend='both') for k in [0, 1]]
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = [ax[k].contourf(
                lon2, lat2, zdata2, 21,transform=ccrs.PlateCarree(),
                levels=np.linspace(-20, 20, 21),
                #levels=np.linspace(30e-13, 60e-13, 21),
                cmap='seismic', extend='both') for k in [2, 3]]
        diffzdata = zdata2-zdata1
        diffrho = 100*(zdata4-zdata3)/zdata3
        hc = [ax[k].contourf(
                lon2, lat2, diffzdata, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-20, 20, 21), cmap='seismic',
                extend='both') for k in [4, 5]]
        hc = [ax[k].contour(
                lon2, lat2, diffrho, [-10], transform=ccrs.PlateCarree(),
                colors='g',linestyles='-') for k in [4, 5]]

        # wind
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=which_alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2, 'neutral', alt=which_alt)
        for iax in range(6):
            if iax == 0 or iax == 1:
                lon0, lat0, ewind, nwind = (
                        lon1.copy(), lat1.copy(), ewind1.copy(), nwind1.copy())
                lon0, lat0, ewind, nwind = g3ca.convert_vector(
                        lon0, lat0, ewind, nwind, plot_type='polar',
                        projection=projection[iax])
            elif iax == 2 or iax == 3:
                lon0, lat0, ewind, nwind = (
                        lon2.copy(), lat2.copy(), ewind2.copy(), nwind2.copy())
                lon0, lat0, ewind, nwind = g3ca.convert_vector(
                        lon0, lat0, ewind, nwind, plot_type='polar',
                        projection=projection[iax])
            elif iax == 4 or iax == 5:
                lon0, lat0, ewind, nwind = (
                        lon1.copy(), lat1.copy(), ewind2-ewind1, nwind2-nwind1)
                lon0, lat0, ewind, nwind = g3ca.convert_vector(
                        lon0, lat0, ewind, nwind, plot_type='polar',
                        projection=projection[iax])
            hq = ax[iax].quiver(
                    lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
                    color='k', regrid_shape=20)
            # ax.quiverkey(hq, 0.93, 0, 1000, '1000 m/s')
            # hc = plt.colorbar(hc, ticks=np.arange(3, 7)*1e-12)
            # hc.set_label(r'$\rho$ (kg/m$^3$)')
            # ax[iax].scatter(
            #         qlonn, qlatn, color='k', transform=ccrs.PlateCarree())
            # ax[iax].scatter(
            #         qlons, qlats, color='k', transform=ccrs.PlateCarree())
        return
    anim = animation.FuncAnimation(
            fig, animate_vert_wind, interval=200, frames=len(fn1))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, bitrate=1800)
    anim.save(path+'02_vertical_wind_%d_%d.mp4' %(f107,which_alt),writer=writer)
    return


def plot_animation_divrhov_rho(show=True, f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_3_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_3/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    fig = plt.figure(figsize=[8,6])
    def animate_divrhov(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]*1000-which_alt))+2
        # create axis
        ax = list(range(2))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            ax[ins], projection[ins] = gcc.create_map(
                    1, 2, 1+ins, 'polar', nlat=nlat, slat=slat,
                    dlat=10, centrallon=g3ca.calculate_centrallon(
                        g1, 'polar',  useLT=True),
                    coastlines=False)

        lon1 = np.array(g1['Longitude'])
        lat1 = np.array(g1['Latitude'])
        alt1 = np.array(g1['Altitude'])
        rho1 = np.array(g1['Rho'])
        omega = (2*np.pi)/(24*3600)
        RR = 6371*1000+alt1
        nwind1 = np.array(g1['V!Dn!N (north)'])
        ewind1 = np.array(g1['V!Dn!N (east)']) + omega*RR*np.cos(lat1)
        uwind1 = np.array(g1['V!Dn!N (up)'])
        div_rhov1 = (cr.calc_div_hozt(lon1, lat1, alt1, rho1*nwind1, rho1*ewind1)\
                   +cr.calc_div_vert(alt1, rho1*uwind1))/rho1

        lon2 = np.array(g2['Longitude'])
        lat2 = np.array(g2['Latitude'])
        alt2 = np.array(g2['Altitude'])
        rho2 = np.array(g2['Rho'])
        omega = (2*np.pi)/(24*3600)
        RR = 6371*1000+alt2
        nwind2 = np.array(g2['V!Dn!N (north)'])
        ewind2 = np.array(g2['V!Dn!N (east)']) + omega*RR*np.cos(lat2)
        uwind2 = np.array(g2['V!Dn!N (up)'])
        div_rhov2 = (cr.calc_div_hozt(lon2, lat2, alt2, rho2*nwind2, rho2*ewind2)\
                   +cr.calc_div_vert(alt2, rho2*uwind2))/rho2
        g1['divrhov'] = div_rhov1 - div_rhov2

        lon0, lat0, div_rhov_d = g3ca.contour_data('divrhov', g1, alt=which_alt)

        hc = ax[0].contourf(
                lon0, lat0, div_rhov_d, np.linspace(-1,1,21)*1e-4,
                transform=ccrs.PlateCarree(),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon0, lat0, div_rhov_d, np.linspace(-1,1,21)*1e-4,
                transform=ccrs.PlateCarree(),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        lon3, lat3, zdata3 = g3ca.contour_data('Rho', g1, alt=which_alt)
        lon4, lat4, zdata4 = g3ca.contour_data('Rho', g2, alt=which_alt)
        diffrho = 100*(zdata4-zdata3)/zdata3
        hc = [ax[k].contour(
                lon0, lat0, diffrho, [-10], transform=ccrs.PlateCarree(),
                colors='g',linestyles='-') for k in [0, 1]]

        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=which_alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2, 'neutral', alt=which_alt)
        lon0, lat0, ewind, nwind = (
                lon1.copy(), lat1.copy(), ewind2-ewind1, nwind2-nwind1)
        for ins in [0, 1]:
            lon0t, lat0t, ewindt, nwindt = g3ca.convert_vector(
                    lon0, lat0, ewind, nwind, plot_type='polar',
                    projection=projection[ins])
            hq = ax[ins].quiver(
                    lon0t, lat0t, ewindt, nwindt, scale=1500,
                    scale_units='inches', color='k', regrid_shape=20)
    anim = animation.FuncAnimation(
            fig, animate_divrhov, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, bitrate=1800)
    anim.save(path+'03_divrhov_rho_%d_%d.mp4' %(f107, which_alt),writer=writer)
    return


def plot_animation_vert_divv(f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_3_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_3/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    fig = plt.figure(figsize=[8,6])
    def animate_vert_divv(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]/1000-which_alt))+2
        # create axis
        ax = list(range(2))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            ax[ins], projection[ins] = gcc.create_map(
                    1, 2, 1+ins, 'polar', nlat=nlat, slat=slat,
                    dlat=10, centrallon=g3ca.calculate_centrallon(
                        g1, 'polar',  useLT=True),
                    coastlines=False)

        # density change (shrink)
        alt1 = np.array(g1['Altitude'])
        uwind1 = np.array(g1['V!Dn!N (up)'])
        divv1 = cr.calc_div_vert(alt1,uwind1)

        # density change (no shrink)
        alt2 = np.array(g2['Altitude'])
        uwind2 = np.array(g2['V!Dn!N (up)'])
        divv2 = cr.calc_div_vert(alt2,uwind2)

        g1['divv'] = divv1-divv2
        lon0,lat0,divv = g3ca.contour_data('divv',g1,alt=which_alt)

        hc = ax[0].contourf(
                lon0, lat0, divv, np.linspace(-1,1,21)*1e-4,
                transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon0, lat0, divv, np.linspace(-1,1,21)*1e-4,
                transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        lon3, lat3, zdata3 = g3ca.contour_data('Rho', g1, alt=which_alt)
        lon4, lat4, zdata4 = g3ca.contour_data('Rho', g2, alt=which_alt)
        diffrho = 100*(zdata4-zdata3)/zdata3
        hc = [ax[k].contour(
                lon0, lat0, diffrho, [-10], transform=ccrs.PlateCarree(),
                colors='g',linestyles='-') for k in [0, 1]]

        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=which_alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2, 'neutral', alt=which_alt)
        lon0, lat0, ewind, nwind = (
                lon1.copy(), lat1.copy(), ewind2-ewind1, nwind2-nwind1)
        for ins in [0, 1]:
            lon0t, lat0t, ewindt, nwindt = g3ca.convert_vector(
                    lon0, lat0, ewind, nwind, plot_type='polar',
                    projection=projection[ins])
            hq = ax[ins].quiver(
                    lon0t, lat0t, ewindt, nwindt, scale=1500,
                    scale_units='inches', color='k', regrid_shape=20)
    anim = animation.FuncAnimation(
            fig, animate_vert_divv, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, bitrate=1800)
    anim.save(path+'04_vert_divv_%d_%d.mp4' %(f107, which_alt),writer=writer)
    return


def plot_animation_hozt_divv(f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_3_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_3/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    fig = plt.figure(figsize=[8,6])
    def animate_hozt_divv(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]/1000-which_alt))+2
        # create axis
        ax = list(range(2))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            ax[ins], projection[ins] = gcc.create_map(
                    1, 2, 1+ins, 'polar', nlat=nlat, slat=slat,
                    dlat=10, centrallon=g3ca.calculate_centrallon(
                        g1, 'polar',  useLT=True),
                    coastlines=False)

        lon1 = np.array(g1['Longitude'])
        lat1 = np.array(g1['Latitude'])
        alt1 = np.array(g1['Altitude'])
        nwind1 = np.array(g1['V!Dn!N (north)'])
        RR = 6371*1000+alt1
        omega = (2*np.pi)/(24*3600)
        ewind1 = np.array(g1['V!Dn!N (east)'])+omega*RR*np.cos(lat1)
        uwind1 = np.array(g1['V!Dn!N (up)'])
        divv1 = cr.calc_div_hozt(lon1,lat1,alt1,nwind1,ewind1)

        lon2 = np.array(g2['Longitude'])
        lat2 = np.array(g2['Latitude'])
        alt2 = np.array(g2['Altitude'])
        nwind2 = np.array(g2['V!Dn!N (north)'])
        RR = 6371*1000+alt2
        omega = (2*np.pi)/(24*3600)
        ewind2 = np.array(g2['V!Dn!N (east)'])+omega*RR*np.cos(lat2)
        uwind2 = np.array(g2['V!Dn!N (up)'])
        divv2 = cr.calc_div_hozt(lon2,lat2,alt2,nwind2,ewind2)

        g1['divv'] = divv1-divv2
        lon0,lat0,divv = g3ca.contour_data('divv',g1,alt=which_alt)

        hc = ax[0].contourf(
                lon0, lat0, divv, np.linspace(-1,1,21)*1e-4,
                transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon0, lat0, divv, np.linspace(-1,1,21)*1e-4,
                transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        lon3, lat3, zdata3 = g3ca.contour_data('Rho', g1, alt=which_alt)
        lon4, lat4, zdata4 = g3ca.contour_data('Rho', g2, alt=which_alt)
        diffrho = 100*(zdata4-zdata3)/zdata3
        hc = [ax[k].contour(
                lon0, lat0, diffrho, [-10], transform=ccrs.PlateCarree(),
                colors='g',linestyles='-') for k in [0, 1]]

        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=which_alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2, 'neutral', alt=which_alt)
        lon0, lat0, ewind, nwind = (
                lon1.copy(), lat1.copy(), ewind2-ewind1, nwind2-nwind1)
        for ins in [0, 1]:
            lon0t, lat0t, ewindt, nwindt = g3ca.convert_vector(
                    lon0, lat0, ewind, nwind, plot_type='polar',
                    projection=projection[ins])
            hq = ax[ins].quiver(
                    lon0t, lat0t, ewindt, nwindt, scale=1500,
                    scale_units='inches', color='k', regrid_shape=20)
    anim = animation.FuncAnimation(
            fig, animate_hozt_divv, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, bitrate=1800)
    anim.save(path+'05_hozt_divv_%d_%d.mp4' %(f107, which_alt),writer=writer)
    return


def plot_animation_vert_vgradrho_rho(f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_3_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_3/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    fig = plt.figure(figsize=[8,6])
    def animate_vgradrho_rho(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]/1000-which_alt))+2
        # create axis
        ax = list(range(2))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            ax[ins], projection[ins] = gcc.create_map(
                    1, 2, 1+ins, 'polar', nlat=nlat, slat=slat,
                    dlat=10, centrallon=g3ca.calculate_centrallon(
                        g1, 'polar',  useLT=True),
                    coastlines=False)

        lon1 = np.array(g1['Longitude'])
        lat1 = np.array(g1['Latitude'])
        alt1 = np.array(g1['Altitude'])
        nwind1 = np.array(g1['V!Dn!N (north)'])
        RR = 6371*1000+alt1
        omega = (2*np.pi)/(24*3600)
        ewind1 = np.array(g1['V!Dn!N (east)'])+omega*RR*np.cos(lat1)
        uwind1 = np.array(g1['V!Dn!N (up)'])
        rho1 = np.array(g1['Rho'])
        vgradrho1 = uwind1*cr.calc_rusanov_alts_ausm(alt1,rho1)/rho1

        lon2 = np.array(g2['Longitude'])
        lat2 = np.array(g2['Latitude'])
        alt2 = np.array(g2['Altitude'])
        nwind2 = np.array(g2['V!Dn!N (north)'])
        RR = 6371*1000+alt2
        omega = (2*np.pi)/(24*3600)
        ewind2 = np.array(g2['V!Dn!N (east)'])+omega*RR*np.cos(lat2)
        uwind2 = np.array(g2['V!Dn!N (up)'])
        rho2 = np.array(g2['Rho'])
        vgradrho2 = uwind2*cr.calc_rusanov_alts_ausm(alt2,rho2)/rho2

        g1['vgradrho'] = vgradrho1-vgradrho2
        lon0,lat0,vgradrho = g3ca.contour_data('vgradrho',g1,alt=which_alt)

        hc = ax[0].contourf(
                lon0, lat0, vgradrho, np.linspace(-1,1,21)*1e-4,
                transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon0, lat0, vgradrho, np.linspace(-1,1,21)*1e-4,
                transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        lon3, lat3, zdata3 = g3ca.contour_data('Rho', g1, alt=which_alt)
        lon4, lat4, zdata4 = g3ca.contour_data('Rho', g2, alt=which_alt)
        diffrho = 100*(zdata4-zdata3)/zdata3
        hc = [ax[k].contour(
                lon0, lat0, diffrho, [-10], transform=ccrs.PlateCarree(),
                colors='g',linestyles='-') for k in [0, 1]]

        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=which_alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2, 'neutral', alt=which_alt)
        lon0, lat0, ewind, nwind = (
                lon1.copy(), lat1.copy(), ewind2-ewind1, nwind2-nwind1)
        for ins in [0, 1]:
            lon0t, lat0t, ewindt, nwindt = g3ca.convert_vector(
                    lon0, lat0, ewind, nwind, plot_type='polar',
                    projection=projection[ins])
            hq = ax[ins].quiver(
                    lon0t, lat0t, ewindt, nwindt, scale=1500,
                    scale_units='inches', color='k', regrid_shape=20)
    anim = animation.FuncAnimation(
            fig, animate_vgradrho_rho, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, bitrate=1800)
    anim.save(path+'06_vert_vgradrho_rho_%d_%d.mp4' %(f107, which_alt),writer=writer)
    return


def plot_animation_hozt_vgradrho_rho(f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_3_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_3/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    fig = plt.figure(figsize=[8,6])
    def animate_vgradrho_rho(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, 2:-2]/1000-which_alt))+2
        # create axis
        ax = list(range(2))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            ax[ins], projection[ins] = gcc.create_map(
                    1, 2, 1+ins, 'polar', nlat=nlat, slat=slat,
                    dlat=10, centrallon=g3ca.calculate_centrallon(
                        g1, 'polar',  useLT=True),
                    coastlines=False)

        lon1 = np.array(g1['Longitude'])
        lat1 = np.array(g1['Latitude'])
        alt1 = np.array(g1['Altitude'])
        nwind1 = np.array(g1['V!Dn!N (north)'])
        RR = 6371*1000+alt1
        omega = (2*np.pi)/(24*3600)
        ewind1 = np.array(g1['V!Dn!N (east)'])+omega*RR*np.cos(lat1)
        uwind1 = np.array(g1['V!Dn!N (up)'])
        rho1 = np.array(g1['Rho'])
        vgradrho1 = nwind1*cr.calc_rusanov_lats(lat1,alt1,rho1)/rho1 \
                   +ewind1*cr.calc_rusanov_lons(lon1,lat1,alt1,rho1)/rho1

        lon2 = np.array(g2['Longitude'])
        lat2 = np.array(g2['Latitude'])
        alt2 = np.array(g2['Altitude'])
        nwind2 = np.array(g2['V!Dn!N (north)'])
        RR = 6371*1000+alt2
        omega = (2*np.pi)/(24*3600)
        ewind2 = np.array(g2['V!Dn!N (east)'])+omega*RR*np.cos(lat2)
        uwind2 = np.array(g2['V!Dn!N (up)'])
        rho2 = np.array(g2['Rho'])
        vgradrho2 = nwind2*cr.calc_rusanov_lats(lat2,alt2,rho2)/rho2 \
                   +ewind2*cr.calc_rusanov_lons(lon2,lat2,alt2,rho2)/rho2

        g1['vgradrho'] = vgradrho1-vgradrho2
        lon0,lat0,vgradrho = g3ca.contour_data('vgradrho',g1,alt=which_alt)

        hc = ax[0].contourf(
                lon0, lat0, vgradrho, np.linspace(-1,1,21)*1e-4,
                transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon0, lat0, vgradrho, np.linspace(-1,1,21)*1e-4,
                transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        lon3, lat3, zdata3 = g3ca.contour_data('Rho', g1, alt=which_alt)
        lon4, lat4, zdata4 = g3ca.contour_data('Rho', g2, alt=which_alt)
        diffrho = 100*(zdata4-zdata3)/zdata3
        hc = [ax[k].contour(
                lon0, lat0, diffrho, [-10], transform=ccrs.PlateCarree(),
                colors='g',linestyles='-') for k in [0, 1]]

        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=which_alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2, 'neutral', alt=which_alt)
        lon0, lat0, ewind, nwind = (
                lon1.copy(), lat1.copy(), ewind2-ewind1, nwind2-nwind1)
        for ins in [0, 1]:
            lon0t, lat0t, ewindt, nwindt = g3ca.convert_vector(
                    lon0, lat0, ewind, nwind, plot_type='polar',
                    projection=projection[ins])
            hq = ax[ins].quiver(
                    lon0t, lat0t, ewindt, nwindt, scale=1500,
                    scale_units='inches', color='k', regrid_shape=20)
    anim = animation.FuncAnimation(
            fig, animate_vgradrho_rho, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, bitrate=1800)
    anim.save(path+'07_hozt_vgradrho_rho_%d_%d.mp4' %(f107, which_alt),writer=writer)
    return


if __name__=='__main__':
    import gc
    plt.close('all')
    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
           'why_no_low_density_cell_at_high_latitude/animation/'
    plot_animation_den_win(f107=150,which_alt=400)
    gc.collect()
