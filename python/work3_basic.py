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
import gitm_divergence_new as gd
import gitm_gradient_new as gg
import calc_rusanov as cr
sns.set('paper', 'whitegrid')

def plot_den_win():
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    g = g2a
    alts = [130, 400]
    for ialt, alt in enumerate(alts):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, 2:-2]/1000-alt))+2
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                1, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        # density
        lon0, lat0, zdata0 = g3ca.contour_data('Rho', g, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, 21,
                         transform=ccrs.PlateCarree(), cmap='jet',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$\rho$ (kg/m$^3$)')
        # wind
        lon0, lat0, ewind, nwind = g3ca.vector_data(g, 'neutral', alt=alt)
        lon0, lat0, ewind, nwind = g3ca.convert_vector(
                lon0, lat0, ewind, nwind, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
                regrid_shape=20)
        ax.quiverkey(hq, 0.93, -0.1, 1000, '1000 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+titletime, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    plt.show()
    return


def plot_den_win_diff():
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    #plt.figure(figsize=(7.26, 9.25))
    plt.figure(figsize=(8.97, 3.45))
    alts = [130, 400]
    for ialt, alt in enumerate(alts):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                1, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        # density difference
        lon1, lat1, zdata1 = g3ca.contour_data('Rho', g1a, alt=alt)
        lon2, lat2, zdata2 = g3ca.contour_data('Rho', g2a, alt=alt)
        fp = (lat1[:,0]>slat) & (lat1[:,0]<nlat)
        lon0, lat0, zdata0 = lon1[fp, :], lat1[fp,:], \
                             (100*(zdata2-zdata1)/zdata1)[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, np.linspace(-20,20,21),
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$100*\frac{\rho2-\rho1}{\rho1}$ (%)')
        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=500, scale_units='inches',
                regrid_shape=20,headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 500, '500m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        # rho difference
        # lon3, lat3, zdata3 = g3ca.contour_data('Rho', g1a, alt=alt)
        # lon4, lat4, zdata4 = g3ca.contour_data('Rho', g2a, alt=alt)
        # diffrho = 100*(zdata4-zdata3)/zdata3
        # hc = ax.contour(
        #         lon4, lat4, diffrho, [-10], transform=ccrs.PlateCarree(),
        #         colors='g',linestyles='-')

        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+titletime, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    plt.show()
    return


def plot_temperature():
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    g = g2a
    alts = [130, 400]
    for ialt, alt in enumerate(alts):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, 2:-2]/1000-alt))+2
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                1, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('Temperature', g, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, 21,
                         transform=ccrs.PlateCarree(), cmap='jet',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label('Temperature (K)')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+titletime, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    plt.show()
    return


def plot_temperature_diff():
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    #plt.figure(figsize=(7.26, 9.25))
    plt.figure(figsize=(10.3,4.3))
    alts = [130, 400]
    for ialt, alt in enumerate(alts):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                1, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        # temperature diff
        lon1, lat1, zdata1 = g3ca.contour_data('Temperature', g1a, alt=alt)
        lon2, lat2, zdata2 = g3ca.contour_data('Temperature', g2a, alt=alt)
        fp = (lat1[:,0]>slat) & (lat1[:,0]<nlat)
        lon0, lat0, zdata0 = lon1[fp, :], lat1[fp,:], (zdata2-zdata1)[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, np.linspace(-80,80,21),
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$T_2-T_1$ (K)')
        # geomagnetic pole
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        # wind difference
        # lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        # lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        # lon0, lat0 = lon1, lat1
        # lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
        #     lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
        #     projection=projection)
        # hq = ax.quiver(
        #     lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
        #     regrid_shape=20, headwidth=5)
        # ax.quiverkey(hq, 0.93, -0.1, 300, '500 m/s')
        # rho difference
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+titletime, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    plt.show()
    return


def plot_vert_vgradrho_rho_diff(show=True, save=True):
    rho1 = np.array(g1a['Rho'])
    vgradrho1 = \
        g1a['V!Dn!N (up)']\
        *cr.calc_rusanov_alts_ausm(g1a['Altitude'],rho1)/g1a['Rho']

    rho2 = np.array(g2a['Rho'])
    vgradrho2 = \
        g2a['V!Dn!N (up)']\
        *cr.calc_rusanov_alts_ausm(g2a['Altitude'],rho2)/g2a['Rho']

    g1a['vgradrho_diff'] = vgradrho1-vgradrho2

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    #plt.figure(figsize=(7.26, 9.25))
    plt.figure(figsize=(8.97, 3.45))
    alts = [130, 400]
    for ialt, alt in enumerate(alts):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
            1, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
            coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('vgradrho_diff', g1a, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(
            lon0, lat0, zdata0, np.linspace(-1,1,21)*1e-4,
            transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        hcb = plt.colorbar(hc, pad=0.17)
        hcb.formatter.set_powerlimits((0,0))
        hcb.update_ticks()
        hcb.set_label(r'$-\vec{u}\cdot\frac{\nabla\rho}{\rho}$ (up)')

        # lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        # lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        # lon0, lat0 = lon1, lat1
        # lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
        #         lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
        #         projection=projection)
        # hq = ax.quiver(
        #         lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
        #         regrid_shape=20, headwidth=5)
        # ax.quiverkey(hq, 0.93, -0.1, 300, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+titletime, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    plt.show()
    return


if __name__=='__main__':
    import gc
    import gitm
    import gitm_create_coordinate as gcc
    from PyPDF2 import PdfFileMerger, PdfFileReader
    Re = 6371*1000 # Earth radius, unit: m

    time = '2003-03-24 00:00:00'
    pdtime = pd.Timestamp(time)
    fntime = pdtime.strftime('%y%m%d_%H%M')
    titletime = pdtime.strftime('%H%M')

    filename = glob.glob(
        '/home/guod/simulation_output/momentum_analysis/'\
        'run_shrink_iondrift_4_all_day/data/3DALL_t%s*.bin'%fntime)[0]
    g1a = gitm.GitmBin(filename)
    filename = glob.glob(
        '/home/guod/simulation_output/momentum_analysis/'\
        'run_no_shrink_iondrift_4_all_day/data/3DALL_t%s*.bin'%fntime)[0]
    g2a = gitm.GitmBin(filename)

    nlat, slat = -30, -90
    plot_den_win_diff()
