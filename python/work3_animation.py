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
sns.set('paper', 'whitegrid')

def plot_animation_den_win(show=False, f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_2_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_2/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        rholevels = np.linspace(6e-12, 10e-12, 21)
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        rholevels = np.linspace(0.5e-12, 1.5e-12, 21)

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
           + 'why_no_low_density_cell_at_high_latitude/animation/'
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
                lon1, lat1, zdata1, 21, transform=ccrs.PlateCarree(),
                levels=rholevels,
                #levels=np.linspace(15e-13, 25e-13, 21),
                cmap='jet', extend='both') for k in [0, 1]]
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = [ax[k].contourf(
                lon2, lat2, zdata2, 21,transform=ccrs.PlateCarree(),
                levels=rholevels,
                #levels=np.linspace(30e-13, 60e-13, 21),
                cmap='jet', extend='both') for k in [2, 3]]
        diffzdata = 100*(zdata2-zdata1)/zdata1
        hc = [ax[k].contourf(
                lon2, lat2, diffzdata, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-30, 30, 21), cmap='seismic',
                extend='both') for k in [4, 5]]

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
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'den_wind_70.wmv',writer=writer)
    return


def plot_animation_vert_wind(show=False):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_2_continue/data/'
    fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fp2 = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_2/data/'
    fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
           + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,9])
    # read gitm data
    def animate_vert_wind(i):
        g1, g2 = [gitm.GitmBin(k[i], varlist=[
            'Rho', 'V!Dn!N (north)', 'V!Dn!N (east)', 'V!Dn!N (up)',
            'V!Di!N (east)', 'V!Di!N (north)']) for k in [fn1, fn2]]
        which_alt=400
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
        hc = [ax[k].contourf(
                lon2, lat2, diffzdata, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-20, 20, 21), cmap='seismic',
                extend='both') for k in [4, 5]]

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
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vertical_wind.mp4',writer=writer)
    return


def plot_animation_all_forces(show=True):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_no_shrink_iondrift_2/data/'
    fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
           + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])
    # read gitm data
    def animate_all_forces(i):
        g1 = gitm.GitmBin(fn1[i])
        alt = 400
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-alt))
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
        # all forces
        lon0 = np.array(g1['dLon'][2:-2, 0, 0])
        lat0 = np.array(g1['dLat'][0, 2:-2, 0])
        nallf = np.array((g1['NeuPressureGrad (north)'] +
                          g1['IonDragForce (north)'] +
                          g1['CoriolisForce (north)'] +
                          g1['ViscosityForce (north)'] +
                          g1['VGradV (north)'] +
                          g1['SpheGeomForce (north)'] +
                          g1['CentriForce (north)']
                          )[2:-2, 2:-2, alt_ind])
        eallf = np.array((g1['NeuPressureGrad (east)'] +
                          g1['IonDragForce (east)'] +
                          g1['CoriolisForce (east)'] +
                          g1['ViscosityForce (east)'] +
                          g1['VGradV (east)'] +
                          g1['SpheGeomForce (east)']
                          )[2:-2, 2:-2, alt_ind])
        nallf, lon0 = add_cyclic_point(nallf.T, coord=lon0, axis=1)
        eallf = add_cyclic_point(eallf.T, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)

        lon0t, lat0t, eallft, nallft = g3ca.convert_vector(
                lon0, lat0, eallf, nallf, plot_type='polar',
                projection=projection[0])
        hq = ax[0].quiver(
                lon0t, lat0t, eallft, nallft, scale=0.3, scale_units='inches',
                regrid_shape=20)
        ax[0].quiverkey(hq, 0.93, -0.1, 0.2, '0.2 $m/s^2$')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

        lon0t, lat0t, eallft, nallft = g3ca.convert_vector(
                lon0, lat0, eallf, nallf, plot_type='polar',
                projection=projection[1])
        hq = ax[1].quiver(
                lon0t, lat0t, eallft, nallft, scale=0.3, scale_units='inches',
                regrid_shape=20)
        ax[1].quiverkey(hq, 0.93, -0.1, 0.2, '0.2 $m/s^2$')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    anim = animation.FuncAnimation(
            fig, animate_all_forces, interval=200, frames=len(fn1))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'all_forces.mp4',writer=writer)
    return


def plot_animation_density_change_shrink(show=True):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_2_continue/data/'
    fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_density_change(i):
        g1 = gitm.GitmBin(fn1[i])
        which_alt = 400
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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
        div_rhov1 = (
                1.0/(RR**2)
              * np.gradient((RR**2)*rho1*uwind1, axis=2) / np.gradient(alt1, axis=2)
              + 1.0/(RR*np.cos(lat1))
              * (np.gradient(rho1*nwind1*np.cos(lat1), axis=1) / np.gradient(lat1, axis=1)
                 + np.gradient(rho1*ewind1, axis=0) / np.gradient(lon1, axis=0)))
        lon1 = np.array(g1['dLon'][2:-2, 0, 0])
        lat1 = np.array(g1['dLat'][0, 2:-2, 0])
        div_rhov1 = div_rhov1[2:-2, 2:-2, alt_ind]
        div_rhov1, lon1 = add_cyclic_point(div_rhov1.T, coord=lon1, axis=1)
        lon1, lat1 = np.meshgrid(lon1, lat1)
        div_rhov1 = - div_rhov1

        hc = ax[0].contourf(
                lon1, lat1, div_rhov1, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon1, lat1, div_rhov1, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    anim = animation.FuncAnimation(
            fig, animate_density_change, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'density_change_shrink.mp4',writer=writer)
    return


def plot_animation_density_change():
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_shrink_iondrift_2_continue/data/'
    fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fp2 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_no_shrink_iondrift_2/data/'
    fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_density_change(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        which_alt = 400
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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
        div_rhov2 = (
                1.0/(RR**2)
              * np.gradient((RR**2)*rho2*uwind2, axis=2) / np.gradient(alt2, axis=2)
              + 1.0/(RR*np.cos(lat2))
              * (np.gradient(rho2*nwind2*np.cos(lat2), axis=1) / np.gradient(lat2, axis=1)
                 + np.gradient(rho2*ewind2, axis=0) / np.gradient(lon2, axis=0)))
        lon2 = np.array(g2['dLon'][2:-2, 0, 0])
        lat2 = np.array(g2['dLat'][0, 2:-2, 0])
        div_rhov2 = div_rhov2[2:-2, 2:-2, alt_ind]
        div_rhov2, lon2 = add_cyclic_point(div_rhov2.T, coord=lon2, axis=1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        div_rhov2 = - div_rhov2

        hc = ax[0].contourf(
                lon2, lat2, div_rhov2, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon2, lat2, div_rhov2, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

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
            fig, animate_density_change, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'density_change.mp4',writer=writer)
    return


def plot_animation_density_change_diff(show=True):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_2_continue/data/'
    fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fp2 = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_2/data/'
    fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_density_change_diff(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        which_alt = 400
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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
        div_rhov1 = (
                1.0/(RR**2)
              * np.gradient((RR**2)*rho1*uwind1, axis=2) / np.gradient(alt1, axis=2)
              + 1.0/(RR*np.cos(lat1))
              * (np.gradient(rho1*nwind1*np.cos(lat1), axis=1) / np.gradient(lat1, axis=1)
                 + np.gradient(rho1*ewind1, axis=0) / np.gradient(lon1, axis=0)))
        lon1 = np.array(g1['dLon'][2:-2, 0, 0])
        lat1 = np.array(g1['dLat'][0, 2:-2, 0])
        div_rhov1 = div_rhov1[2:-2, 2:-2, alt_ind]
        div_rhov1, lon1 = add_cyclic_point(div_rhov1.T, coord=lon1, axis=1)
        lon1, lat1 = np.meshgrid(lon1, lat1)
        div_rhov1 = - div_rhov1

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
        div_rhov2 = (
                1.0/(RR**2)
              * np.gradient((RR**2)*rho2*uwind2, axis=2) / np.gradient(alt2, axis=2)
              + 1.0/(RR*np.cos(lat2))
              * (np.gradient(rho2*nwind2*np.cos(lat2), axis=1) / np.gradient(lat2, axis=1)
                 + np.gradient(rho2*ewind2, axis=0) / np.gradient(lon2, axis=0)))
        lon2 = np.array(g2['dLon'][2:-2, 0, 0])
        lat2 = np.array(g2['dLat'][0, 2:-2, 0])
        div_rhov2 = div_rhov2[2:-2, 2:-2, alt_ind]
        div_rhov2, lon2 = add_cyclic_point(div_rhov2.T, coord=lon2, axis=1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        div_rhov2 = - div_rhov2

        div_rhov_d = div_rhov2 - div_rhov1

        hc = ax[0].contourf(
                lon2, lat2, div_rhov_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon2, lat2, div_rhov_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

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
            fig, animate_density_change_diff, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'density_change_diff.mp4',writer=writer)
    return


def plot_animation_rhodivv(show=True):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_2_continue/data/'
    fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fp2 = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_2/data/'
    fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_rhodivv(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        which_alt = 400
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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
        rhodivv2 = rho2*(
                1.0/(RR**2)
              * np.gradient((RR**2)*uwind2, axis=2) / np.gradient(alt2, axis=2)
              + 1.0/(RR*np.cos(lat2))
              * (np.gradient(nwind2*np.cos(lat2), axis=1) / np.gradient(lat2, axis=1)
                 + np.gradient(ewind2, axis=0) / np.gradient(lon2, axis=0)))
        lon2 = np.array(g2['dLon'][2:-2, 0, 0])
        lat2 = np.array(g2['dLat'][0, 2:-2, 0])
        rhodivv2 = rhodivv2[2:-2, 2:-2, alt_ind]
        rhodivv2, lon2 = add_cyclic_point(rhodivv2.T, coord=lon2, axis=1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        rhodivv2 = - rhodivv2

        hc = ax[0].contourf(
                lon2, lat2, rhodivv2, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon2, lat2, rhodivv2, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

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
            fig, animate_rhodivv, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'rhodivv.mp4',writer=writer)
    return


def plot_animation_rhodivv_diff(f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_2_continue/data/'
    fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fp2 = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_2/data/'
    fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_rhodivv_diff(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        which_alt = 400
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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
        rhodivv1 = rho1*(
                1.0/(RR**2)
              * np.gradient((RR**2)*uwind1, axis=2) / np.gradient(alt1, axis=2)
              + 1.0/(RR*np.cos(lat1))
              * (np.gradient(nwind1*np.cos(lat1), axis=1) / np.gradient(lat1, axis=1)
                 + np.gradient(ewind1, axis=0) / np.gradient(lon1, axis=0)))
        lon1 = np.array(g1['dLon'][2:-2, 0, 0])
        lat1 = np.array(g1['dLat'][0, 2:-2, 0])
        rhodivv1 = rhodivv1[2:-2, 2:-2, alt_ind]
        rhodivv1, lon1 = add_cyclic_point(rhodivv1.T, coord=lon1, axis=1)
        lon1, lat1 = np.meshgrid(lon1, lat1)
        rhodivv1 = - rhodivv1

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
        rhodivv2 = rho2*(
                1.0/(RR**2)
              * np.gradient((RR**2)*uwind2, axis=2) / np.gradient(alt2, axis=2)
              + 1.0/(RR*np.cos(lat2))
              * (np.gradient(nwind2*np.cos(lat2), axis=1) / np.gradient(lat2, axis=1)
                 + np.gradient(ewind2, axis=0) / np.gradient(lon2, axis=0)))
        lon2 = np.array(g2['dLon'][2:-2, 0, 0])
        lat2 = np.array(g2['dLat'][0, 2:-2, 0])
        rhodivv2 = rhodivv2[2:-2, 2:-2, alt_ind]
        rhodivv2, lon2 = add_cyclic_point(rhodivv2.T, coord=lon2, axis=1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        rhodivv2 = - rhodivv2

        rhodivv_d = rhodivv2-rhodivv1

        hc = ax[0].contourf(
                lon2, lat2, rhodivv_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon2, lat2, rhodivv_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

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
            fig, animate_rhodivv_diff, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'rhodivv_diff.mp4',writer=writer)
    return


def plot_animation_vert_rhodivv_diff(f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_2_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_2/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_vert_rhodivv_diff(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        which_alt = 400
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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
        rhodivv1 = rho1*(1.0/(RR**2) * np.gradient((RR**2)*uwind1, axis=2)\
                         / np.gradient(alt1, axis=2))
        lon1 = np.array(g1['dLon'][2:-2, 0, 0])
        lat1 = np.array(g1['dLat'][0, 2:-2, 0])
        rhodivv1 = rhodivv1[2:-2, 2:-2, alt_ind]
        rhodivv1, lon1 = add_cyclic_point(rhodivv1.T, coord=lon1, axis=1)
        lon1, lat1 = np.meshgrid(lon1, lat1)
        rhodivv1 = - rhodivv1

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
        rhodivv2 = rho2*(1.0/(RR**2) * np.gradient((RR**2)*uwind2, axis=2) \
                         / np.gradient(alt2, axis=2))
        lon2 = np.array(g2['dLon'][2:-2, 0, 0])
        lat2 = np.array(g2['dLat'][0, 2:-2, 0])
        rhodivv2 = rhodivv2[2:-2, 2:-2, alt_ind]
        rhodivv2, lon2 = add_cyclic_point(rhodivv2.T, coord=lon2, axis=1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        rhodivv2 = - rhodivv2

        rhodivv_d = rhodivv2-rhodivv1

        hc = ax[0].contourf(
                lon2, lat2, rhodivv_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon2, lat2, rhodivv_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

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
            fig, animate_vert_rhodivv_diff, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'rhodivv_vert_diff.mp4',writer=writer)
    return


def plot_animation_hozt_rhodivv_diff(f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_2_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_2/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_hozt_rhodivv_diff(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        which_alt = 400
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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
        rhodivv1 = rho1 * (1.0/(RR*np.cos(lat1))
                * (np.gradient(nwind1*np.cos(lat1), axis=1)
                / np.gradient(lat1, axis=1)
                + np.gradient(ewind1, axis=0) / np.gradient(lon1, axis=0)))
        lon1 = np.array(g1['dLon'][2:-2, 0, 0])
        lat1 = np.array(g1['dLat'][0, 2:-2, 0])
        rhodivv1 = rhodivv1[2:-2, 2:-2, alt_ind]
        rhodivv1, lon1 = add_cyclic_point(rhodivv1.T, coord=lon1, axis=1)
        lon1, lat1 = np.meshgrid(lon1, lat1)
        rhodivv1 = - rhodivv1

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
        rhodivv2 = rho2 * (1.0/(RR*np.cos(lat2))
                * (np.gradient(nwind2*np.cos(lat2), axis=1)
                / np.gradient(lat2, axis=1)
                + np.gradient(ewind2, axis=0) / np.gradient(lon2, axis=0)))
        lon2 = np.array(g2['dLon'][2:-2, 0, 0])
        lat2 = np.array(g2['dLat'][0, 2:-2, 0])
        rhodivv2 = rhodivv2[2:-2, 2:-2, alt_ind]
        rhodivv2, lon2 = add_cyclic_point(rhodivv2.T, coord=lon2, axis=1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        rhodivv2 = - rhodivv2

        rhodivv_d = rhodivv2-rhodivv1

        hc = ax[0].contourf(
                lon2, lat2, rhodivv_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon2, lat2, rhodivv_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

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
            fig, animate_hozt_rhodivv_diff, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'rhodivv_hozt_diff.mp4',writer=writer)
    return


def plot_animation_vgradrho(show=True):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_2_continue/data/'
    fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fp2 = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_2/data/'
    fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_vgradrho(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        which_alt = 400
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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
        vgradrho2 = \
                uwind2 * (np.gradient(rho2, axis=2)
                          / np.gradient(alt2, axis=2))\
              + nwind2 * ((1.0/RR)*np.gradient(rho2, axis=1)
                          / np.gradient(lat2, axis=1))\
              + ewind2 * ((1.0/(RR*np.cos(lat2)))*np.gradient(rho2, axis=0)
                          / np.gradient(lon2, axis=0))
        lon2 = np.array(g2['dLon'][2:-2, 0, 0])
        lat2 = np.array(g2['dLat'][0, 2:-2, 0])
        vgradrho2 = vgradrho2[2:-2, 2:-2, alt_ind]
        vgradrho2, lon2 = add_cyclic_point(vgradrho2.T, coord=lon2, axis=1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        vgradrho2 = - vgradrho2

        hc = ax[0].contourf(
                lon2, lat2, vgradrho2, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon2, lat2, vgradrho2, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

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
            fig, animate_vgradrho, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vgradrho.mp4',writer=writer)
    return


def plot_animation_vgradrho_diff(show=True):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_shrink_iondrift_2_continue/data/'
    fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fp2 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_no_shrink_iondrift_2/data/'
    fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_vgradrho_diff(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        which_alt = 400
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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

        # density change (no shrink)
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
        vgradrho1 = \
                uwind1 * (np.gradient(rho1, axis=2)
                          / np.gradient(alt1, axis=2))\
              + nwind1 * ((1.0/RR)*np.gradient(rho1, axis=1)
                          / np.gradient(lat1, axis=1))\
              + ewind1 * ((1.0/(RR*np.cos(lat1)))*np.gradient(rho1, axis=0)
                          / np.gradient(lon1, axis=0))
        lon1 = np.array(g1['dLon'][2:-2, 0, 0])
        lat1 = np.array(g1['dLat'][0, 2:-2, 0])
        vgradrho1 = vgradrho1[2:-2, 2:-2, alt_ind]
        vgradrho1, lon1 = add_cyclic_point(vgradrho1.T, coord=lon1, axis=1)
        lon1, lat1 = np.meshgrid(lon1, lat1)
        vgradrho1 = - vgradrho1

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
        vgradrho2 = \
                uwind2 * (np.gradient(rho2, axis=2)
                          / np.gradient(alt2, axis=2))\
              + nwind2 * ((1.0/RR)*np.gradient(rho2, axis=1)
                          / np.gradient(lat2, axis=1))\
              + ewind2 * ((1.0/(RR*np.cos(lat2)))*np.gradient(rho2, axis=0)
                          / np.gradient(lon2, axis=0))
        lon2 = np.array(g2['dLon'][2:-2, 0, 0])
        lat2 = np.array(g2['dLat'][0, 2:-2, 0])
        vgradrho2 = vgradrho2[2:-2, 2:-2, alt_ind]
        vgradrho2, lon2 = add_cyclic_point(vgradrho2.T, coord=lon2, axis=1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        vgradrho2 = - vgradrho2

        vgradrho_d = vgradrho2-vgradrho1

        hc = ax[0].contourf(
                lon2, lat2, vgradrho_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon2, lat2, vgradrho_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

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
            fig, animate_vgradrho_diff, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vgradrho_diff.mp4',writer=writer)
    return


def plot_animation_vgradrho_vert_diff(f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_2_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_2/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_vgradrho_vert_diff(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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
        vgradrho1 = uwind1 * (np.gradient(rho1, axis=2)
                    / np.gradient(alt1, axis=2))
        lon1 = np.array(g1['dLon'][2:-2, 0, 0])
        lat1 = np.array(g1['dLat'][0, 2:-2, 0])
        vgradrho1 = vgradrho1[2:-2, 2:-2, alt_ind]
        vgradrho1, lon1 = add_cyclic_point(vgradrho1.T, coord=lon1, axis=1)
        lon1, lat1 = np.meshgrid(lon1, lat1)
        vgradrho1 = - vgradrho1

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
        vgradrho2 = uwind2 * (np.gradient(rho2, axis=2)
                    / np.gradient(alt2, axis=2))
        lon2 = np.array(g2['dLon'][2:-2, 0, 0])
        lat2 = np.array(g2['dLat'][0, 2:-2, 0])
        vgradrho2 = vgradrho2[2:-2, 2:-2, alt_ind]
        vgradrho2, lon2 = add_cyclic_point(vgradrho2.T, coord=lon2, axis=1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        vgradrho2 = - vgradrho2

        vgradrho_d = vgradrho2-vgradrho1

        hc = ax[0].contourf(
                lon2, lat2, vgradrho_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon2, lat2, vgradrho_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

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
            fig, animate_vgradrho_vert_diff, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vgradrho_vert_diff_%d.mp4' % f107,writer=writer)
    return


def plot_animation_vgradrho_hozt_diff(f107=150, which_alt=400):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    if f107==150:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_iondrift_2_continue/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_iondrift_2/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
    else:
        fp1 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_shrink_70/data/'
        fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]
        fp2 = '/home/guod/simulation_output/momentum_analysis/'\
              + 'run_no_shrink_70/data/'
        fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
               for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_vgradrho_hozt_diff(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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
        vgradrho1 = nwind1 * ((1.0/RR)*np.gradient(rho1, axis=1)
                              / np.gradient(lat1, axis=1))\
                  + ewind1 * ((1.0/(RR*np.cos(lat1)))*np.gradient(rho1, axis=0)
                              / np.gradient(lon1, axis=0))
        lon1 = np.array(g1['dLon'][2:-2, 0, 0])
        lat1 = np.array(g1['dLat'][0, 2:-2, 0])
        vgradrho1 = vgradrho1[2:-2, 2:-2, alt_ind]
        vgradrho1, lon1 = add_cyclic_point(vgradrho1.T, coord=lon1, axis=1)
        lon1, lat1 = np.meshgrid(lon1, lat1)
        vgradrho1 = - vgradrho1

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
        vgradrho2 = nwind2 * ((1.0/RR)*np.gradient(rho2, axis=1)
                              / np.gradient(lat2, axis=1))\
                  + ewind2 * ((1.0/(RR*np.cos(lat2)))*np.gradient(rho2, axis=0)
                              / np.gradient(lon2, axis=0))
        lon2 = np.array(g2['dLon'][2:-2, 0, 0])
        lat2 = np.array(g2['dLat'][0, 2:-2, 0])
        vgradrho2 = vgradrho2[2:-2, 2:-2, alt_ind]
        vgradrho2, lon2 = add_cyclic_point(vgradrho2.T, coord=lon2, axis=1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        vgradrho2 = - vgradrho2

        vgradrho_d = vgradrho2-vgradrho1

        hc = ax[0].contourf(
                lon2, lat2, vgradrho_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon2, lat2, vgradrho_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

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
            fig, animate_vgradrho_hozt_diff, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vgradrho_hozt_diff_%d.mp4' % f107,writer=writer)
    return


def plot_animation_vgradrho_rhodivv(show=True):
    # A test function
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_2_continue/data/'
    fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fp2 = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_2/data/'
    fn2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
         + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])

    # read gitm data
    def animate_vgradrho_rhodivv(i):
        g1 = gitm.GitmBin(fn1[i])
        g2 = gitm.GitmBin(fn2[i])
        which_alt = 400
        alt_ind = np.argmin(np.abs(g1['Altitude'][0, 0, :]/1000-which_alt))
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
        vgradrho2 = \
                uwind2 * (np.gradient(rho2, axis=2)
                          / np.gradient(alt2, axis=2))\
              + nwind2 * ((1.0/RR)*np.gradient(rho2, axis=1)
                          / np.gradient(lat2, axis=1))\
              + ewind2 * ((1.0/(RR*np.cos(lat2)))*np.gradient(rho2, axis=0)
                          / np.gradient(lon2, axis=0))
        rhodivv2 = rho2*(
                1.0/(RR**2)
              * np.gradient((RR**2)*uwind2, axis=2) / np.gradient(alt2, axis=2)
              + 1.0/(RR*np.cos(lat2))
              * (np.gradient(nwind2*np.cos(lat2), axis=1) / np.gradient(lat2, axis=1)
                 + np.gradient(ewind2, axis=0) / np.gradient(lon2, axis=0)))
        lon2 = np.array(g2['dLon'][2:-2, 0, 0])
        lat2 = np.array(g2['dLat'][0, 2:-2, 0])
        vgradrho2 = vgradrho2[2:-2, 2:-2, alt_ind]
        rhodivv2 = rhodivv2[2:-2, 2:-2, alt_ind]
        vgradrho2, lon2t = add_cyclic_point(vgradrho2.T, coord=lon2, axis=1)
        rhodivv2, lon2 = add_cyclic_point(rhodivv2.T, coord=lon2, axis=1)
        lon2, lat2 = np.meshgrid(lon2, lat2)
        vgradrho2 = - vgradrho2
        rhodivv2 = - rhodivv2

        rho_change = vgradrho2 + rhodivv2

        hc = ax[0].contourf(
                lon2, lat2, rho_change, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon2, lat2, rho_change, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)

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
            fig, animate_vgradrho_rhodivv, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vgradrho_rhodivv.mp4',writer=writer)
    return

def plot_animation_vert_all_forces():
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 03:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')

    fp1 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_shrink_iondrift_3_continue/data/'
    fn1a = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn1m = [glob.glob(fp1+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    fp2 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_no_shrink_iondrift_3/data/'
    fn2a = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn2m = [glob.glob(fp2+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
           + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])
    # read gitm data
    def animate_vert_all_forces(i):
        g1a = gitm.GitmBin(fn1a[i])
        g1m = gitm.GitmBin(fn1m[i])
        g2a = gitm.GitmBin(fn2a[i])
        g2m = gitm.GitmBin(fn2m[i])
        alt = 400
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        # create axis
        ax = list(range(2))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            ax[ins], projection[ins] = gcc.create_map(
                    1, 2, 1+ins, 'polar', nlat=nlat, slat=slat,
                    dlat=10, centrallon=g3ca.calculate_centrallon(
                        g1a, 'polar',  useLT=True),
                    coastlines=False)
        # all forces
        lon0 = np.array(g1a['dLon'][2:-2, 0, 0])
        lat0 = np.array(g1a['dLat'][0, 2:-2, 0])

        all_forces1 = g1m['NeuPressureGrad (up)']\
                    + g1m['GravityForce']\
                    + g1m['VGradV (up)']\
                    + g1m['CentriForce (up)']\
                    + g1m['CoriolisForce (up)']\
                    + g1m['SpheGeomForce (up)']
                    #+ g1m['ViscosityForce (up)']
        all_forces2 = g2m['NeuPressureGrad (up)']\
                    + g2m['GravityForce']\
                    + g2m['VGradV (up)']\
                    + g2m['CentriForce (up)']\
                    + g2m['CoriolisForce (up)']\
                    + g2m['SpheGeomForce (up)']
                    #+ g1m['ViscosityForce (up)']
        all_forces = np.array((all_forces2-all_forces1)[2:-2,2:-2,alt_ind])
        zdata0, lon0 = add_cyclic_point(all_forces.T, coord=lon0, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        for ins in range(2):
            ax[ins].contourf(lon0,lat0,zdata0,
                    levels=np.linspace(-0.01, 0.01, 21),
                    transform=ccrs.PlateCarree(),
                    cmap='jet',extend='both')

        ax[1].set_title(g1a['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    anim = animation.FuncAnimation(
            fig, animate_vert_all_forces, interval=200, frames=len(fn1a))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vert_all_forces.mp4',writer=writer)
    return


def plot_animation_vert_presgrad():
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 03:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')

    fp1 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_shrink_iondrift_3_continue/data/'
    fn1a = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn1m = [glob.glob(fp1+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    fp2 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_no_shrink_iondrift_3/data/'
    fn2a = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn2m = [glob.glob(fp2+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
           + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])
    # read gitm data
    def animate_vert_presgrad(i):
        g1a = gitm.GitmBin(fn1a[i])
        g1m = gitm.GitmBin(fn1m[i])
        g2a = gitm.GitmBin(fn2a[i])
        g2m = gitm.GitmBin(fn2m[i])
        alt = 400
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        # create axis
        ax = list(range(2))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            ax[ins], projection[ins] = gcc.create_map(
                    1, 2, 1+ins, 'polar', nlat=nlat, slat=slat,
                    dlat=10, centrallon=g3ca.calculate_centrallon(
                        g1a, 'polar',  useLT=True),
                    coastlines=False)
        # all forces
        lon0 = np.array(g1a['dLon'][2:-2, 0, 0])
        lat0 = np.array(g1a['dLat'][0, 2:-2, 0])

        all_forces1 = g1m['NeuPressureGrad (up)']
                    #+ g1m['GravityForce']\
                    #+ g1m['VGradV (up)']\
                    #+ g1m['CentriForce (up)']\
                    #+ g1m['CoriolisForce (up)']\
                    #+ g1m['SpheGeomForce (up)']
                    #+ g1m['ViscosityForce (up)']
        all_forces2 = g2m['NeuPressureGrad (up)']
                    #+ g2m['GravityForce']\
                    #+ g2m['VGradV (up)']\
                    #+ g2m['CentriForce (up)']\
                    #+ g2m['CoriolisForce (up)']\
                    #+ g2m['SpheGeomForce (up)']
                    #+ g1m['ViscosityForce (up)']
        all_forces = np.array((all_forces2-all_forces1)[2:-2,2:-2,alt_ind])
        zdata0, lon0 = add_cyclic_point(all_forces.T, coord=lon0, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        for ins in range(2):
            ax[ins].contourf(lon0,lat0,zdata0,
                    levels=np.linspace(-0.03, 0.03, 21),
                    transform=ccrs.PlateCarree(),
                    cmap='seismic',extend='both')

        ax[1].set_title(g1a['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    anim = animation.FuncAnimation(
            fig, animate_vert_presgrad, interval=200, frames=len(fn1a))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vert_presgrad.mp4',writer=writer)
    return


def plot_animation_vert_vgradv():
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')

    fp1 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_shrink_iondrift_3_continue/data/'
    fn1a = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn1m = [glob.glob(fp1+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    fp2 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_no_shrink_iondrift_3/data/'
    fn2a = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn2m = [glob.glob(fp2+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
           + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])
    # read gitm data
    def animate_vert_vgradv(i):
        g1a = gitm.GitmBin(fn1a[i])
        g1m = gitm.GitmBin(fn1m[i])
        g2a = gitm.GitmBin(fn2a[i])
        g2m = gitm.GitmBin(fn2m[i])
        alt = 400
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        # create axis
        ax = list(range(2))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            ax[ins], projection[ins] = gcc.create_map(
                    1, 2, 1+ins, 'polar', nlat=nlat, slat=slat,
                    dlat=10, centrallon=g3ca.calculate_centrallon(
                        g1a, 'polar',  useLT=True),
                    coastlines=False)
        # all forces
        lon0 = np.array(g1a['dLon'][2:-2, 0, 0])
        lat0 = np.array(g1a['dLat'][0, 2:-2, 0])

        all_forces1 = g1m['VGradV (up)']
        all_forces2 = g2m['VGradV (up)']
        all_forces = np.array((all_forces2-all_forces1)[2:-2,2:-2,alt_ind])
        zdata0, lon0 = add_cyclic_point(all_forces.T, coord=lon0, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        for ins in range(2):
            ax[ins].contourf(lon0,lat0,zdata0,
                    levels=np.linspace(-0.03, 0.03, 21),
                    transform=ccrs.PlateCarree(),
                    cmap='seismic',extend='both')

        ax[1].set_title(g1a['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    anim = animation.FuncAnimation(
            fig, animate_vert_vgradv, interval=200, frames=len(fn1a))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vert_vgradv.mp4',writer=writer)
    return


def plot_animation_vert_sphegeom():
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')

    fp1 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_shrink_iondrift_3_continue/data/'
    fn1a = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn1m = [glob.glob(fp1+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    fp2 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_no_shrink_iondrift_3/data/'
    fn2a = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn2m = [glob.glob(fp2+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
           + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])
    # read gitm data
    def animate_vert_sphegeom(i):
        g1a = gitm.GitmBin(fn1a[i])
        g1m = gitm.GitmBin(fn1m[i])
        g2a = gitm.GitmBin(fn2a[i])
        g2m = gitm.GitmBin(fn2m[i])
        alt = 400
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        # create axis
        ax = list(range(2))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            ax[ins], projection[ins] = gcc.create_map(
                    1, 2, 1+ins, 'polar', nlat=nlat, slat=slat,
                    dlat=10, centrallon=g3ca.calculate_centrallon(
                        g1a, 'polar',  useLT=True),
                    coastlines=False)
        # all forces
        lon0 = np.array(g1a['dLon'][2:-2, 0, 0])
        lat0 = np.array(g1a['dLat'][0, 2:-2, 0])

        all_forces1 = g1m['SpheGeomForce (up)']
                    #+ g1m['ViscosityForce (up)']
        all_forces2 = g2m['SpheGeomForce (up)']
                    #+ g1m['ViscosityForce (up)']
        all_forces = np.array((all_forces2-all_forces1)[2:-2,2:-2,alt_ind])
        zdata0, lon0 = add_cyclic_point(all_forces.T, coord=lon0, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        for ins in range(2):
            ax[ins].contourf(lon0,lat0,zdata0,
                    levels=np.linspace(-0.03, 0.03, 21),
                    transform=ccrs.PlateCarree(),
                    cmap='seismic',extend='both')

        ax[1].set_title(g1a['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    anim = animation.FuncAnimation(
            fig, animate_vert_sphegeom, interval=200, frames=len(fn1a))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vert_sphegeom.mp4',writer=writer)
    return


def plot_animation_vert_centrifugal():
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')

    fp1 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_shrink_iondrift_3_continue/data/'
    fn1a = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn1m = [glob.glob(fp1+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    fp2 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_no_shrink_iondrift_3/data/'
    fn2a = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn2m = [glob.glob(fp2+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
           + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])
    # read gitm data
    def animate_vert_centrifugal(i):
        g1a = gitm.GitmBin(fn1a[i])
        g1m = gitm.GitmBin(fn1m[i])
        g2a = gitm.GitmBin(fn2a[i])
        g2m = gitm.GitmBin(fn2m[i])
        alt = 400
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        # create axis
        ax = list(range(2))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            ax[ins], projection[ins] = gcc.create_map(
                    1, 2, 1+ins, 'polar', nlat=nlat, slat=slat,
                    dlat=10, centrallon=g3ca.calculate_centrallon(
                        g1a, 'polar',  useLT=True),
                    coastlines=False)
        # all forces
        lon0 = np.array(g1a['dLon'][2:-2, 0, 0])
        lat0 = np.array(g1a['dLat'][0, 2:-2, 0])

        all_forces1 = g1m['CentriForce (up)']
        all_forces2 = g2m['CentriForce (up)']
        all_forces = np.array((all_forces2-all_forces1)[2:-2,2:-2,alt_ind])
        zdata0, lon0 = add_cyclic_point(all_forces.T, coord=lon0, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        for ins in range(2):
            ax[ins].contourf(lon0,lat0,zdata0,
                    levels=np.linspace(-0.03, 0.03, 21),
                    transform=ccrs.PlateCarree(),
                    cmap='seismic',extend='both')

        ax[1].set_title(g1a['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    anim = animation.FuncAnimation(
            fig, animate_vert_centrifugal, interval=200, frames=len(fn1a))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vert_centrifugal.mp4',writer=writer)
    return


def plot_animation_vert_coriolis():
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 03:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')

    fp1 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_shrink_iondrift_3_continue/data/'
    fn1a = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn1m = [glob.glob(fp1+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    fp2 = '/home/guod/simulation_output/momentum_analysis/'\
          + 'run_no_shrink_iondrift_3/data/'
    fn2a = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fn2m = [glob.glob(fp2+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
           + 'why_no_low_density_cell_at_high_latitude/iondrift_with_or_not/'
    fig = plt.figure(figsize=[8,6])
    # read gitm data
    def animate_vert_coriolis(i):
        g1a = gitm.GitmBin(fn1a[i])
        g1m = gitm.GitmBin(fn1m[i])
        g2a = gitm.GitmBin(fn2a[i])
        g2m = gitm.GitmBin(fn2m[i])
        alt = 400
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        # create axis
        ax = list(range(2))
        projection = ax.copy()
        for ins in range(2):
            nlat, slat = [90, 40] if ins==0 else [-40, -90]
            ax[ins], projection[ins] = gcc.create_map(
                    1, 2, 1+ins, 'polar', nlat=nlat, slat=slat,
                    dlat=10, centrallon=g3ca.calculate_centrallon(
                        g1a, 'polar',  useLT=True),
                    coastlines=False)
        # all forces
        lon0 = np.array(g1a['dLon'][2:-2, 0, 0])
        lat0 = np.array(g1a['dLat'][0, 2:-2, 0])

        all_forces1 = g1m['CoriolisForce (up)']\
                    + g1m['SpheGeomForce (up)']
        all_forces2 = g2m['CoriolisForce (up)']\
                    + g2m['SpheGeomForce (up)']
        all_forces = np.array((all_forces2-all_forces1)[2:-2,2:-2,alt_ind])
        zdata0, lon0 = add_cyclic_point(all_forces.T, coord=lon0, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        for ins in range(2):
            ax[ins].contourf(lon0,lat0,zdata0,
                    levels=np.linspace(-0.03, 0.03, 21),
                    transform=ccrs.PlateCarree(),
                    cmap='seismic',extend='both')

        ax[1].set_title(g1a['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    anim = animation.FuncAnimation(
            fig, animate_vert_coriolis, interval=200, frames=len(fn1a))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'vert_coriolis.mp4',writer=writer)
    return


if __name__=='__main__':
    import gc
    plt.close('all')
    plot_animation_den_win()
    gc.collect()
