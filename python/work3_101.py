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

Re = 6371*1000 # Earth radius, unit: m
filename = '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_1'\
           '/data/3DALL_t030322_000000.bin'
path =  '/home/guod/Documents/work/fig/density_cell/'\
        'why_no_low_density_cell_at_high_latitude/150_december_north/'
path =  '/home/guod/tmp/3DALL_t030322_000000.bin'
g = gitm.GitmBin(filename)
gp.calc_pressure(g)
nlat, slat = -40, -90


def plot_den_win(show=True, save=False):
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    # 150 equinox
    # cbl = [np.linspace(0.8, 1.4, 21)*1e-8, np.linspace(2.1, 2.7, 21)*1e-10,
    #        np.linspace(2.6, 4, 21)*1e-11, np.linspace(5.5, 10, 21)*1e-12,
    #        np.linspace(1.2, 3.5, 21)*1e-12, np.linspace(3.5, 9, 21)*1e-13]
    # 150 june
    # cbl = [np.linspace(0.7, 1.4, 21)*1e-8, np.linspace(1.8, 2.5, 21)*1e-10,
    #        np.linspace(2.3, 3.5, 21)*1e-11, np.linspace(3, 7, 21)*1e-12,
    #        np.linspace(0.7, 3, 21)*1e-12, np.linspace(1.5, 6, 21)*1e-13]
    # 150 december
    # cbl = [np.linspace(0.7, 1.4, 21)*1e-8, np.linspace(1.8, 3, 21)*1e-10,
    #        np.linspace(2.8, 4, 21)*1e-11, np.linspace(4, 9, 21)*1e-12,
    #        np.linspace(1, 4, 21)*1e-12, np.linspace(4, 9, 21)*1e-13]
    # 70 equinox
    # cbl = [np.linspace(0.8, 1.4, 21)*1e-8, np.linspace(1.150, 1.6, 21)*1e-10,
    #        np.linspace(0.5, 1, 21)*1e-11, np.linspace(0.5, 1.5, 21)*1e-12,
    #        np.linspace(0.04, 0.2, 21)*1e-12, np.linspace(0.3, 2, 21)*1e-13]
    # ..............
    cbl = [np.linspace(0.8, 1.4, 21)*1e-8, np.linspace(2.3, 2.7, 21)*1e-10,
           np.linspace(3, 4, 21)*1e-11, np.linspace(6, 10, 21)*1e-12,
           np.linspace(1.5, 3.5, 21)*1e-12, np.linspace(4, 9, 21)*1e-13]
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('Rho', g, alt=alt)
        hc = ax.contourf(lon0, lat0, zdata0, cbl[ialt],
                         transform=ccrs.PlateCarree(), cmap='jet',
                         extend='both')
        lon0, lat0, ewind, nwind = g3ca.vector_data(g, 'neutral', alt=alt)
        lon0, lat0, ewind, nwind = g3ca.convert_vector(
                lon0, lat0, ewind, nwind, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
                alpha=0.5, regrid_shape=20)
        ax.quiverkey(hq, 0.93, -0.1, 1000, '1000 m/s')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$\rho$ (kg/m$^3$)')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    if show:
        plt.show()
    if save:
        plt.savefig(path+'den_win.pdf')
    return


def plot_ion_drift():
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, ewind, nwind = g3ca.vector_data(g, 'ion', alt=alt)
        lon0, lat0, ewind, nwind = g3ca.convert_vector(
                lon0, lat0, ewind, nwind, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 1000, '1000 m/s')
        # ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.show()
    return


def plot_temperature(show=True, save=False):
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    # 150 equinox
    # cbl = [np.linspace(400, 500, 21), np.linspace(1000, 1300, 21),
    #        np.linspace(1200, 1550, 21), np.linspace(1200, 1600, 21),
    #        np.linspace(1200, 1600, 21), np.linspace(1200, 1600, 21)]
    # 150 june
    # cbl = [np.linspace(400, 500, 21), np.linspace(800, 1100, 21),
    #        np.linspace(900, 1300, 21), np.linspace(900, 1300, 21),
    #        np.linspace(900, 1300, 21), np.linspace(900, 1300, 21)]
    # 150 december
    cbl = [np.linspace(400, 500, 21), np.linspace(1100, 1400, 21),
           np.linspace(1400, 1800, 21), np.linspace(1500, 1950, 21),
           np.linspace(1500, 1950, 21), np.linspace(1500, 1950, 21)]
    # 70 equinox
    # cbl = [np.linspace(300, 400, 21), np.linspace(700, 900, 21),
    #        np.linspace(800, 1000, 21), np.linspace(800, 1000, 21),
    #        np.linspace(800, 1000, 21), np.linspace(800, 1000, 21)]
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('Temperature', g, alt=alt)
        hc = ax.contourf(lon0, lat0, zdata0, cbl[ialt],
                         transform=ccrs.PlateCarree(), cmap='jet',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label('Temperature (K)')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    if show:
        plt.show()
    if save:
        plt.savefig(path+'temperature.pdf')
    return


def plot_vertical_wind(show=True):
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('V!Dn!N (up)', g, alt=alt)
        hc = ax.contourf(lon0, lat0, zdata0, np.linspace(-20, 20, 21),
                         transform=ccrs.PlateCarree(), cmap='jet',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label('Vertical Wind (m/s)')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    if show:
        plt.show()
    plt.savefig(path+'vertical_wind.pdf')
    return


def plot_pressure(show=True, save=False):
    global g
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    # 150 equinox
    # cbl = [np.linspace(1.3, 2, 21)*1e-3, np.linspace(9, 13, 21)*1e-5,
    #        np.linspace(15, 24, 21)*1e-6, np.linspace(3.5,9, 21)*1e-6,
    #        np.linspace(12, 22, 21)*1e-7, np.linspace(25, 80, 21)*1e-8]
    # 150 june
    # cbl = [np.linspace(1.2, 1.8, 21)*1e-3, np.linspace(7, 11, 21)*1e-5,
    #        np.linspace(11, 24, 21)*1e-6, np.linspace(1.7, 6, 21)*1e-6,
    #        np.linspace(5, 20, 21)*1e-7, np.linspace(15, 60, 21)*1e-8]
    # 150 december
    cbl = [np.linspace(1.2, 1.8, 21)*1e-3, np.linspace(8, 13, 21)*1e-5,
           np.linspace(16, 28, 21)*1e-6, np.linspace(3, 6, 21)*1e-6,
           np.linspace(10, 25, 21)*1e-7, np.linspace(40, 80, 21)*1e-8]
    # 70 equinox
    # cbl = [np.linspace(1, 1.8, 21)*1e-3, np.linspace(3.5, 5, 21)*1e-5,
    #        np.linspace(1, 6, 21)*1e-6, np.linspace(0.2,1.2, 21)*1e-6,
    #        np.linspace(0.3, 1.2, 21)*1e-7, np.linspace(2, 15, 21)*1e-8]
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('pressure', g, alt=alt)
        hc = ax.contourf(lon0, lat0, zdata0, cbl[ialt],
                         transform=ccrs.PlateCarree(), cmap='jet',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17, format='%.1e')
        hc.set_label('pressure (Pa)')
        #ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    if show:
        plt.show()
    if save:
        plt.savefig(path+'pressure.pdf')
    return


def plot_ave_m(show=True, save=False):
    global g
    Re = 6371*1000 # Earth radius, unit: m
    ave_m = ((g['O(!U3!NP)']*16 + g['O!D2!N']*32 + g['N!D2!N']*28 +
              g['N(!U4!NS)']*14 + g['NO']*30 + g['N(!U2!ND)']*14 +
              g['N(!U2!NP)']*14 + g['H']*2 + g['He']*4 + g['CO!D2!N']*44 +
              g['O(!U1!ND)']*16) /
             (g['O(!U3!NP)']+g['O!D2!N']+g['N!D2!N']+g['N(!U4!NS)']+g['NO']+
              g['N(!U2!ND)']+g['N(!U2!NP)']+g['H']+g['He']+g['CO!D2!N']+
              g['O(!U1!ND)']))
    g['ave_m'] = dmarray(
            ave_m, attrs={'units':'', 'scale':'linear', 'name':'ave_m'})
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    # 150 equinox
    # cbl = [np.linspace(26, 27, 21), np.linspace(20, 23, 21),
    #        np.linspace(17, 20, 21), np.linspace(16,18, 21),
    #        np.linspace(16, 17, 21), np.linspace(15, 16.3, 21)]
    # 150 june
    # cbl = [np.linspace(25, 26, 21), np.linspace(18, 22, 21),
    #        np.linspace(15, 18, 21), np.linspace(12,17, 21),
    #        np.linspace(10, 17, 21), np.linspace(9, 15, 21)]
    # 150 december
    cbl = [np.linspace(27, 28, 21), np.linspace(23, 25, 21),
           np.linspace(20, 23, 21), np.linspace(18,20, 21),
           np.linspace(17, 18, 21), np.linspace(16, 17, 21)]
    # 70 equinox
    # cbl = [np.linspace(26.5, 27.5, 21), np.linspace(20, 23, 21),
    #        np.linspace(17, 18, 21), np.linspace(15,16, 21),
    #        np.linspace(13, 14, 21), np.linspace(13, 14, 21)]
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('ave_m', g, alt=alt)
        hc = ax.contourf(lon0, lat0, zdata0, #cbl[ialt],
                         transform=ccrs.PlateCarree(), cmap='jet',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label('Ave M')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    if show:
        plt.show()
    if save:
        plt.savefig(path+'ave_m.pdf')
    return


def plot_pressure_gradient(show=True, save=False):
    global g
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    # cbl = [np.linspace(23, 25, 21), np.linspace(20, 23, 21),
    #        np.linspace(17, 20, 21), np.linspace(16,18, 21),
    #        np.linspace(16, 17, 21), np.linspace(16, 16.3, 21)]
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0 = np.array(g['dLon'][2:-2, 0, 0])
        lat0 = np.array(g['dLat'][0, 2:-2, 0])
        npgf = np.array(g['NeuPressureGrad (north)'][2:-2, 2:-2, alt_ind])
        epgf = np.array(g['NeuPressureGrad (east)'][2:-2, 2:-2, alt_ind])
        npgf, lon0 = add_cyclic_point(npgf.T, coord=lon0, axis=1)
        epgf = add_cyclic_point(epgf.T, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        lon0, lat0, epgf, npgf = g3ca.convert_vector(
                lon0, lat0, epgf, npgf, plot_type='polar',
                projection=projection)
        # hq = ax.quiver(
        #         lon0, lat0, epgf, npgf, scale=0.3, scale_units='inches',
        #         regrid_shape=15)
        hq = ax.quiver(
                lon0, lat0, epgf, npgf, scale=0.1, scale_units='inches',
                regrid_shape=15, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 0.2, r'0.2 m/$s^2$')
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Pressure Gradient Force', fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(path+'pressure_grad.pdf')
    return


def plot_coriolis(show=True, save=False):
    global g
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    # cbl = [np.linspace(23, 25, 21), np.linspace(20, 23, 21),
    #        np.linspace(17, 20, 21), np.linspace(16,18, 21),
    #        np.linspace(16, 17, 21), np.linspace(16, 16.3, 21)]
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0 = np.array(g['dLon'][2:-2, 0, 0])
        lat0 = np.array(g['dLat'][0, 2:-2, 0])
        ncf = np.array(g['CoriolisForce (north)'][2:-2, 2:-2, alt_ind])
        ecf = np.array(g['CoriolisForce (east)'][2:-2, 2:-2, alt_ind])
        ncf, lon0 = add_cyclic_point(ncf.T, coord=lon0, axis=1)
        ecf = add_cyclic_point(ecf.T, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        lon0, lat0, ecf, ncf = g3ca.convert_vector(
                lon0, lat0, ecf, ncf, plot_type='polar',
                projection=projection)
        # hq = ax.quiver(
        #         lon0, lat0, ecf, ncf, scale=0.3, scale_units='inches',
        #         regrid_shape=15)
        hq = ax.quiver(
                lon0, lat0, ecf, ncf, scale=0.1, scale_units='inches',
                regrid_shape=15, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 0.2, r'0.2 m/$s^2$')
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Coriolis', fontsize=15, horizontalalignment='center',
             transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(path+'coriolis.pdf')
    return


def plot_ion_drag(show=True, save=False):
    global g
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    # cbl = [np.linspace(23, 25, 21), np.linspace(20, 23, 21),
    #        np.linspace(17, 20, 21), np.linspace(16,18, 21),
    #        np.linspace(16, 17, 21), np.linspace(16, 16.3, 21)]
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0 = np.array(g['dLon'][2:-2, 0, 0])
        lat0 = np.array(g['dLat'][0, 2:-2, 0])
        nidf = np.array(g['IonDragForce (north)'][2:-2, 2:-2, alt_ind])
        eidf = np.array(g['IonDragForce (east)'][2:-2, 2:-2, alt_ind])
        nidf, lon0 = add_cyclic_point(nidf.T, coord=lon0, axis=1)
        eidf = add_cyclic_point(eidf.T, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        lon0, lat0, eidf, nidf = g3ca.convert_vector(
                lon0, lat0, eidf, nidf, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, eidf, nidf, scale=0.3, scale_units='inches',
                regrid_shape=15)
        ax.quiverkey(hq, 0.93, -0.1, 0.2, '0.2 $m/s^2$')
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Ion Drag', fontsize=15, horizontalalignment='center',
             transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(path+'ion_drag.pdf')
    return


def plot_viscosity(show=True, save=False):
    global g
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    # cbl = [np.linspace(23, 25, 21), np.linspace(20, 23, 21),
    #        np.linspace(17, 20, 21), np.linspace(16,18, 21),
    #        np.linspace(16, 17, 21), np.linspace(16, 16.3, 21)]
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0 = np.array(g['dLon'][2:-2, 0, 0])
        lat0 = np.array(g['dLat'][0, 2:-2, 0])
        nvf = np.array(g['ViscosityForce (north)'][2:-2, 2:-2, alt_ind])
        evf = np.array(g['ViscosityForce (east)'][2:-2, 2:-2, alt_ind])
        nvf, lon0 = add_cyclic_point(nvf.T, coord=lon0, axis=1)
        evf = add_cyclic_point(evf.T, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        lon0, lat0, evf, nvf = g3ca.convert_vector(
                lon0, lat0, evf, nvf, plot_type='polar',
                projection=projection)
        # hq = ax.quiver(
        #         lon0, lat0, evf, nvf, scale=0.3, scale_units='inches',
        #         regrid_shape=15)
        hq = ax.quiver(
                lon0, lat0, evf, nvf, scale=0.1, scale_units='inches',
                regrid_shape=15, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 0.2, '0.2 $m/s^2$')
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Viscosity Force', fontsize=15, horizontalalignment='center',
             transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(path+'viscosity.pdf')
    return


def plot_vgradv(show=True, save=False):
    global g
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    # cbl = [np.linspace(23, 25, 21), np.linspace(20, 23, 21),
    #        np.linspace(17, 20, 21), np.linspace(16,18, 21),
    #        np.linspace(16, 17, 21), np.linspace(16, 16.3, 21)]
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0 = np.array(g['dLon'][2:-2, 0, 0])
        lat0 = np.array(g['dLat'][0, 2:-2, 0])
        nf = np.array((g['VGradV (north)']+g['SpheGeomForce (north)']+g['CentriForce (north)'])
                [2:-2, 2:-2, alt_ind])
        ef = np.array((g['VGradV (east)']+g['SpheGeomForce (east)'])
                [2:-2, 2:-2, alt_ind])
        nf = np.array((g['CentriForce (north)'])[2:-2, 2:-2, alt_ind])
        ef = np.array((g['SpheGeomForce (east)'])[2:-2, 2:-2, alt_ind])*0
        nf, lon0 = add_cyclic_point(nf.T, coord=lon0, axis=1)
        ef = add_cyclic_point(ef.T, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        lon0, lat0, ef, nf = g3ca.convert_vector(
                lon0, lat0, ef, nf, plot_type='polar',
                projection=projection)
        # hq = ax.quiver(
        #         lon0, lat0, ef, nf, scale=0.3, scale_units='inches',
        #         regrid_shape=15)
        hq = ax.quiver(
                lon0, lat0, ef, nf, scale=0.1, scale_units='inches',
                regrid_shape=15, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 0.2, '0.2 $m/s^2$')
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'VGradV', fontsize=15, horizontalalignment='center',
             transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(path+'vgradv.pdf')
    return


def plot_all_forces(show=True, save=False):
    global g
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    # cbl = [np.linspace(23, 25, 21), np.linspace(20, 23, 21),
    #        np.linspace(17, 20, 21), np.linspace(16,18, 21),
    #        np.linspace(16, 17, 21), np.linspace(16, 16.3, 21)]
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0 = np.array(g['dLon'][2:-2, 0, 0])
        lat0 = np.array(g['dLat'][0, 2:-2, 0])
        nallf = np.array((g['NeuPressureGrad (north)'] +
                          g['IonDragForce (north)'] +
                          g['CoriolisForce (north)'] +
                          g['ViscosityForce (north)'] +
                          g['VGradV (north)'] +
                          g['SpheGeomForce (north)'] +
                          g['CentriForce (north)']
                          )[2:-2, 2:-2, alt_ind])
        eallf = np.array((g['NeuPressureGrad (east)'] +
                          g['IonDragForce (east)'] +
                          g['CoriolisForce (east)'] +
                          g['ViscosityForce (east)'] +
                          g['VGradV (east)'] +
                          g['SpheGeomForce (east)']
                          )[2:-2, 2:-2, alt_ind])
        nallf, lon0 = add_cyclic_point(nallf.T, coord=lon0, axis=1)
        eallf = add_cyclic_point(eallf.T, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        lon0, lat0, eallf, nallf = g3ca.convert_vector(
                lon0, lat0, eallf, nallf, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, eallf, nallf, scale=0.3, scale_units='inches',
                regrid_shape=15)
        ax.quiverkey(hq, 0.93, -0.1, 0.2, '0.2 $m/s^2$')
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'All Forces', fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(path+'all_forces.pdf')
    return


def plot_corio_vgradv(show=True):
    global g
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    # cbl = [np.linspace(23, 25, 21), np.linspace(20, 23, 21),
    #        np.linspace(17, 20, 21), np.linspace(16,18, 21),
    #        np.linspace(16, 17, 21), np.linspace(16, 16.3, 21)]
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0 = np.array(g['dLon'][2:-2, 0, 0])
        lat0 = np.array(g['dLat'][0, 2:-2, 0])
        nallf = np.array((#g['NeuPressureGrad (north)'] +
                          g['IonDragForce (north)'] +
                          g['CoriolisForce (north)'] +
                          #g['ViscosityForce (north)'] +
                          g['VGradV (north)'] +
                          g['SpheGeomForce (north)'] +
                          g['CentriForce (north)']
                          )[2:-2, 2:-2, alt_ind])
        eallf = np.array((#g['NeuPressureGrad (east)'] +
                          g['IonDragForce (east)'] +
                          g['CoriolisForce (east)'] +
                          #g['ViscosityForce (east)'] +
                          g['VGradV (east)'] +
                          g['SpheGeomForce (east)']
                          )[2:-2, 2:-2, alt_ind])
        nallf, lon0 = add_cyclic_point(nallf.T, coord=lon0, axis=1)
        eallf = add_cyclic_point(eallf.T, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        lon0, lat0, eallf, nallf = g3ca.convert_vector(
                lon0, lat0, eallf, nallf, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, eallf, nallf, scale=0.3, scale_units='inches',
                regrid_shape=15)
        ax.quiverkey(hq, 0.93, -0.1, 0.2, '0.2 $m/s^2$')
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'All Forces', fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    plt.savefig(path+'corio_vgradv.pdf')
    return


def presgrad_70vs150():
    fn70 = '/home/guod/simulation_output/momentum_analysis/'\
           'run_70_equinox/data/3DALL_t030322_000000.bin'
    fn150 = '/home/guod/simulation_output/momentum_analysis/'\
            'run_150_equinox/data/3DALL_t030322_000001.bin'
    g70 = gitm.GitmBin(fn70)
    g150 = gitm.GitmBin(fn150)
    lat = g70['Latitude']
    lat1 = g70['dLat'][0, :, 0]
    alt70 = g70['Altitude'][0,0,2:-2]/1000
    alt150 = g150['Altitude'][0,0,2:-2]/1000

    pg70n = g70['NeuPressureGrad (north)'][:, lat1<=-40, 2:-2]
    pg70e = g70['NeuPressureGrad (east)'][:, lat1<=-40, 2:-2]
    pg70 = np.sqrt(pg70n**2+pg70e**2)
    pg150n = g150['NeuPressureGrad (north)'][:, lat1<=-40, 2:-2]
    pg150e = g150['NeuPressureGrad (east)'][:, lat1<=-40, 2:-2]
    pg150 = np.sqrt(pg150n**2+pg150e**2)
    lat = lat[:, lat1>=40, 2:-2]

    avepg70 = np.ones(lat.shape[2])
    avepg150 = np.ones(lat.shape[2])
    for ialt, valt in enumerate(alt70):
        avepg70[ialt] = np.mean(pg70[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                        np.mean(np.cos(lat[:,:,ialt]))
        avepg150[ialt] = np.mean(pg150[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
    plt.close('all')
    plt.figure()
    plt.plot(avepg70, alt70,label='F107 = 70')
    plt.plot(avepg150, alt150, label='F107 = 150')
    plt.legend()
    plt.show()
    return


def iondrag_70vs150():
    fn70 = '/home/guod/simulation_output/momentum_analysis/'\
           'run_70_equinox/data/3DALL_t030322_000000.bin'
    fn150 = '/home/guod/simulation_output/momentum_analysis/'\
            'run_150_equinox/data/3DALL_t030322_000001.bin'
    g70 = gitm.GitmBin(fn70)
    g150 = gitm.GitmBin(fn150)
    lat = g70['Latitude']
    lat1 = g70['dLat'][0, :, 0]
    alt70 = g70['Altitude'][0,0,2:-2]/1000
    alt150 = g150['Altitude'][0,0,2:-2]/1000

    iondrag70n = g70['IonDragForce (north)'][:, lat1<=-40, 2:-2]
    iondrag70e = g70['IonDragForce (east)'][:, lat1<=-40, 2:-2]
    iondrag70 = np.sqrt(iondrag70n**2+iondrag70e**2)
    iondrag150n = g150['IonDragForce (north)'][:, lat1<=-40, 2:-2]
    iondrag150e = g150['IonDragForce (east)'][:, lat1<=-40, 2:-2]
    iondrag150 = np.sqrt(iondrag150n**2+iondrag150e**2)
    lat = lat[:, lat1>=40, 2:-2]

    aveiondrag70 = np.ones(lat.shape[2])
    aveiondrag150 = np.ones(lat.shape[2])
    for ialt, valt in enumerate(alt70):
        aveiondrag70[ialt] = np.mean(iondrag70[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                        np.mean(np.cos(lat[:,:,ialt]))
        aveiondrag150[ialt] = np.mean(iondrag150[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
    plt.close('all')
    plt.figure()
    plt.plot(aveiondrag70, alt70,label='F107 = 70')
    plt.plot(aveiondrag150, alt150, label='F107 = 150')
    plt.legend()
    plt.show()
    return


def coriolis_70vs150():
    fn70 = '/home/guod/simulation_output/momentum_analysis/'\
           'run_70_equinox/data/3DALL_t030322_000000.bin'
    fn150 = '/home/guod/simulation_output/momentum_analysis/'\
            'run_150_equinox/data/3DALL_t030322_000001.bin'
    g70 = gitm.GitmBin(fn70)
    g150 = gitm.GitmBin(fn150)
    lat = g70['Latitude']
    lat1 = g70['dLat'][0, :, 0]
    alt70 = g70['Altitude'][0,0,2:-2]/1000
    alt150 = g150['Altitude'][0,0,2:-2]/1000

    coriolis70n = g70['CoriolisForce (north)'][:, lat1<=-40, 2:-2]
    coriolis70e = g70['CoriolisForce (east)'][:, lat1<=-40, 2:-2]
    coriolis70 = np.sqrt(coriolis70n**2+coriolis70e**2)
    coriolis150n = g150['CoriolisForce (north)'][:, lat1<=-40, 2:-2]
    coriolis150e = g150['CoriolisForce (east)'][:, lat1<=-40, 2:-2]
    coriolis150 = np.sqrt(coriolis150n**2+coriolis150e**2)
    lat = lat[:, lat1>=40, 2:-2]

    avecoriolis70 = np.ones(lat.shape[2])
    avecoriolis150 = np.ones(lat.shape[2])
    for ialt, valt in enumerate(alt70):
        avecoriolis70[ialt] = np.mean(coriolis70[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                        np.mean(np.cos(lat[:,:,ialt]))
        avecoriolis150[ialt] = np.mean(coriolis150[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
    plt.close('all')
    plt.figure()
    plt.plot(avecoriolis70, alt70,label='F107 = 70')
    plt.plot(avecoriolis150, alt150, label='F107 = 150')
    plt.legend()
    plt.show()
    return


def vgradv_70vs150():
    fn70 = '/home/guod/simulation_output/momentum_analysis/'\
           'run_70_equinox/data/3DALL_t030322_000000.bin'
    fn150 = '/home/guod/simulation_output/momentum_analysis/'\
            'run_150_equinox/data/3DALL_t030322_000001.bin'
    g70 = gitm.GitmBin(fn70)
    g150 = gitm.GitmBin(fn150)
    lat = g70['Latitude']
    lat1 = g70['dLat'][0, :, 0]
    alt70 = g70['Altitude'][0,0,2:-2]/1000
    alt150 = g150['Altitude'][0,0,2:-2]/1000

    vgradv70n = (g70['VGradV (north)']+g70['SpheGeomForce (north)']+
                 g70['CentriForce (north)'])[:, lat1<=-40, 2:-2]
    vgradv70e = (g70['VGradV (east)']+g70['SpheGeomForce (east)'])\
                [:, lat1<=-40, 2:-2]
    vgradv70 = np.sqrt(vgradv70n**2+vgradv70e**2)
    vgradv150n = (g150['VGradV (north)']+g150['SpheGeomForce (north)']+
                 g150['CentriForce (north)'])[:, lat1<=-40, 2:-2]
    vgradv150e = (g150['VGradV (east)']+g150['SpheGeomForce (east)'])\
                 [:, lat1<=-40, 2:-2]
    vgradv150 = np.sqrt(vgradv150n**2+vgradv150e**2)
    lat = lat[:, lat1>=40, 2:-2]

    avevgradv70 = np.ones(lat.shape[2])
    avevgradv150 = np.ones(lat.shape[2])
    for ialt, valt in enumerate(alt70):
        avevgradv70[ialt] = np.mean(vgradv70[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                        np.mean(np.cos(lat[:,:,ialt]))
        avevgradv150[ialt] = np.mean(vgradv150[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
    plt.close('all')
    plt.figure()
    plt.plot(avevgradv70, alt70,label='F107 = 70')
    plt.plot(avevgradv150, alt150, label='F107 = 150')
    plt.legend()
    plt.show()
    return


def viscosity_70vs150():
    fn70 = '/home/guod/simulation_output/momentum_analysis/'\
           'run_70_equinox/data/3DALL_t030322_000000.bin'
    fn150 = '/home/guod/simulation_output/momentum_analysis/'\
            'run_150_equinox/data/3DALL_t030322_000001.bin'
    g70 = gitm.GitmBin(fn70)
    g150 = gitm.GitmBin(fn150)
    lat = g70['Latitude']
    lat1 = g70['dLat'][0, :, 0]
    alt70 = g70['Altitude'][0,0,2:-2]/1000
    alt150 = g150['Altitude'][0,0,2:-2]/1000

    vis70n = g70['ViscosityForce (north)'][:, lat1<=-40, 2:-2]
    vis70e = g70['ViscosityForce (east)'][:, lat1<=-40, 2:-2]
    vis70 = np.sqrt(vis70n**2+vis70e**2)
    vis150n = g150['ViscosityForce (north)'][:, lat1<=-40, 2:-2]
    vis150e = g150['ViscosityForce (east)'][:, lat1<=-40, 2:-2]
    vis150 = np.sqrt(vis150n**2+vis150e**2)
    lat = lat[:, lat1>=40, 2:-2]

    avevis70 = np.ones(lat.shape[2])
    avevis150 = np.ones(lat.shape[2])
    for ialt, valt in enumerate(alt70):
        avevis70[ialt] = np.mean(vis70[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
        avevis150[ialt] = np.mean(vis150[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                          np.mean(np.cos(lat[:,:,ialt]))
    plt.close('all')
    plt.figure()
    plt.plot(avevis70, alt70, '*-', label='F107 = 70')
    plt.plot(avevis150, alt150, '*-', label='F107 = 150')
    plt.legend()
    plt.show()
    return


def allforces_70vs150():
    fn70 = '/home/guod/simulation_output/momentum_analysis/'\
           'run_70_equinox/data/3DALL_t030322_000000.bin'
    fn150 = '/home/guod/simulation_output/momentum_analysis/'\
            'run_150_equinox/data/3DALL_t030322_000001.bin'
    g70 = gitm.GitmBin(fn70)
    g150 = gitm.GitmBin(fn150)
    lat = g70['Latitude']
    lat1 = g70['dLat'][0, :, 0]
    lt1 = g70['LT'][:, 0, 0]
    alt70 = g70['Altitude'][0,0,2:-2]/1000
    alt150 = g150['Altitude'][0,0,2:-2]/1000
    lt1p = (lt1>=0) & (lt1<=24)
    lat1p = lat1<-40

    #viscosity
    vis70n = g70['ViscosityForce (north)'][lt1p][:, lat1p, 2:-2]
    vis70e = g70['ViscosityForce (east)'][lt1p][:, lat1p, 2:-2]
    vis70 = np.sqrt(vis70n**2+vis70e**2)
    vis150n = g150['ViscosityForce (north)'][lt1p][:, lat1p, 2:-2]
    vis150e = g150['ViscosityForce (east)'][lt1p][:, lat1p, 2:-2]
    vis150 = np.sqrt(vis150n**2+vis150e**2)
    #pressure gradient
    pg70n = g70['NeuPressureGrad (north)'][lt1p][:, lat1p, 2:-2]
    pg70e = g70['NeuPressureGrad (east)'][lt1p][:, lat1p, 2:-2]
    pg70 = np.sqrt(pg70n**2+pg70e**2)
    pg150n = g150['NeuPressureGrad (north)'][lt1p][:, lat1p, 2:-2]
    pg150e = g150['NeuPressureGrad (east)'][lt1p][:, lat1p, 2:-2]
    pg150 = np.sqrt(pg150n**2+pg150e**2)
    #ion drag
    iondrag70n = g70['IonDragForce (north)'][lt1p][:, lat1p, 2:-2]
    iondrag70e = g70['IonDragForce (east)'][lt1p][:, lat1p, 2:-2]
    iondrag70 = np.sqrt(iondrag70n**2+iondrag70e**2)
    iondrag150n = g150['IonDragForce (north)'][lt1p][:, lat1p, 2:-2]
    iondrag150e = g150['IonDragForce (east)'][lt1p][:, lat1p, 2:-2]
    iondrag150 = np.sqrt(iondrag150n**2+iondrag150e**2)
    #vgradv
    vgradv70n = (g70['VGradV (north)']+g70['SpheGeomForce (north)']+
            g70['CentriForce (north)'])[lt1p][:, lat1p, 2:-2]
    vgradv70e = (g70['VGradV (east)']+g70['SpheGeomForce (east)'])\
            [lt1p][:, lat1p, 2:-2]
    vgradv70 = np.sqrt(vgradv70n**2+vgradv70e**2)
    vgradv150n = (g150['VGradV (north)']+g150['SpheGeomForce (north)']+
            g150['CentriForce (north)'])[lt1p][:, lat1p, 2:-2]
    vgradv150e = (g150['VGradV (east)']+g150['SpheGeomForce (east)'])\
            [lt1p][:, lat1p, 2:-2]
    vgradv150 = np.sqrt(vgradv150n**2+vgradv150e**2)
    #coriolis
    coriolis70n = g70['CoriolisForce (north)'][lt1p][:, lat1p, 2:-2]
    coriolis70e = g70['CoriolisForce (east)'][lt1p][:, lat1p, 2:-2]
    coriolis70 = np.sqrt(coriolis70n**2+coriolis70e**2)
    coriolis150n = g150['CoriolisForce (north)'][lt1p][:, lat1p, 2:-2]
    coriolis150e = g150['CoriolisForce (east)'][lt1p][:, lat1p, 2:-2]
    coriolis150 = np.sqrt(coriolis150n**2+coriolis150e**2)
    lat = lat[lt1p][:, lat1>=40, 2:-2]

    avevis70 = np.ones(lat.shape[2])
    avevis150 = np.ones(lat.shape[2])
    avepg70 = np.ones(lat.shape[2])
    avepg150 = np.ones(lat.shape[2])
    aveiondrag70 = np.ones(lat.shape[2])
    aveiondrag150 = np.ones(lat.shape[2])
    avevgradv70 = np.ones(lat.shape[2])
    avevgradv150 = np.ones(lat.shape[2])
    avecoriolis70 = np.ones(lat.shape[2])
    avecoriolis150 = np.ones(lat.shape[2])
    for ialt, valt in enumerate(alt70):
        avevis70[ialt] = np.mean(vis70[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
        avevis150[ialt] = np.mean(vis150[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                          np.mean(np.cos(lat[:,:,ialt]))
        avepg70[ialt] = np.mean(pg70[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
        avepg150[ialt] = np.mean(pg150[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                          np.mean(np.cos(lat[:,:,ialt]))
        aveiondrag70[ialt] = np.mean(iondrag70[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
        aveiondrag150[ialt] = np.mean(iondrag150[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                          np.mean(np.cos(lat[:,:,ialt]))
        avevgradv70[ialt] = np.mean(vgradv70[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
        avevgradv150[ialt] = np.mean(vgradv150[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                          np.mean(np.cos(lat[:,:,ialt]))
        avecoriolis70[ialt] = np.mean(coriolis70[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
        avecoriolis150[ialt] = np.mean(coriolis150[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                          np.mean(np.cos(lat[:,:,ialt]))
    plt.close('all')
    plt.figure()
    plt.plot(avepg70, alt70, 'ro-', label='F107 = 70, PresGrad')
    #plt.plot(avepg150, alt150, 'r*-', label='F107 = 150, PresGrad')
    plt.plot(avevis70, alt70, 'bo-', label='F107 = 70, Viscosity')
    #plt.plot(avevis150, alt150, 'b*-', label='F107 = 150, Viscosity')
    #plt.plot(aveiondrag70, alt70, 'ko-', label='F107 = 70, Ion drag')
    #plt.plot(aveiondrag150, alt150, 'k*-', label='F107 = 150, Ion drag')
    plt.plot(avevgradv70, alt70, 'yo-', label='F107 = 70, VGradV')
    #plt.plot(avevgradv150, alt150, 'y*-', label='F107 = 150, VGradV')
    plt.plot(avecoriolis70, alt70, 'go-', label='F107 = 70, Coriolis')
    #plt.plot(avecoriolis150, alt150, 'g*-', label='F107 = 150, VGradV')
    #    plt.xlim([0,0.1])
    #    plt.ylim([0,700])
    #    plt.xlabel(r'Forces (m/$s^2$)')
    #    plt.ylabel('Altitude (km)')
    #    plt.legend()
    #    plt.show()
    #    plt.savefig('/home/guod/Desktop/why_no_low_density_cell_at_high_latitude/'\
    #                'allforces_70.pdf')
    #    plt.figure()
    #plt.plot(avepg70, alt70, 'ro-', label='F107 = 70, PresGrad')
    plt.plot(avepg150, alt150, 'r*--', label='F107 = 150, PresGrad')
    #plt.plot(avevis70, alt70, 'bo-', label='F107 = 70, Viscosity')
    plt.plot(avevis150, alt150, 'b*--', label='F107 = 150, Viscosity')
    #plt.plot(aveiondrag70, alt70, 'ko-', label='F107 = 70, Ion drag')
    #plt.plot(aveiondrag150, alt150, 'k*-', label='F107 = 150, Ion drag')
    #plt.plot(avevgradv70, alt70, 'yo-', label='F107 = 70, VGradV')
    plt.plot(avevgradv150, alt150, 'y*--', label='F107 = 150, VGradV')
    #plt.plot(avecoriolis70, alt70, 'go-', label='F107 = 70, VGradV')
    plt.plot(avecoriolis150, alt150, 'g*--', label='F107 = 150, Coriolis')
    plt.xlim([0,0.15])
    plt.ylim([0,700])
    plt.xlabel(r'Forces (m/$s^2$)')
    plt.ylabel('Altitude (km)')
    plt.legend()
    plt.show()
    plt.savefig('/home/guod/Desktop/why_no_low_density_cell_at_high_latitude/'\
                'allforces_70vs150.pdf')
    return


def allforces_dec_june():
    fndec = '/home/guod/simulation_output/momentum_analysis/'\
           'run_150_december/data/3DALL_t031222_000000.bin'
    fnjun = '/home/guod/simulation_output/momentum_analysis/'\
            'run_150_june/data/3DALL_t030622_000000.bin'
    gdec = gitm.GitmBin(fndec)
    gjun = gitm.GitmBin(fnjun)
    lat = gdec['Latitude']
    lat1 = gdec['dLat'][0, :, 0]
    lt1 = gdec['LT'][:, 0, 0]
    altdec = gdec['Altitude'][0,0,2:-2]/1000
    altjun = gjun['Altitude'][0,0,2:-2]/1000
    lt1p = (lt1>=0) & (lt1<=24)
    lat1p = lat1<-40

    #viscosity
    visdecn = gdec['ViscosityForce (north)'][lt1p][:, lat1p, 2:-2]
    visdece = gdec['ViscosityForce (east)'][lt1p][:, lat1p, 2:-2]
    visdec = np.sqrt(visdecn**2+visdece**2)
    visjunn = gjun['ViscosityForce (north)'][lt1p][:, lat1p, 2:-2]
    visjune = gjun['ViscosityForce (east)'][lt1p][:, lat1p, 2:-2]
    visjun = np.sqrt(visjunn**2+visjune**2)
    #pressure gradient
    pgdecn = gdec['NeuPressureGrad (north)'][lt1p][:, lat1p, 2:-2]
    pgdece = gdec['NeuPressureGrad (east)'][lt1p][:, lat1p, 2:-2]
    pgdec = np.sqrt(pgdecn**2+pgdece**2)
    pgjunn = gjun['NeuPressureGrad (north)'][lt1p][:, lat1p, 2:-2]
    pgjune = gjun['NeuPressureGrad (east)'][lt1p][:, lat1p, 2:-2]
    pgjun = np.sqrt(pgjunn**2+pgjune**2)
    #ion drag
    iondragdecn = gdec['IonDragForce (north)'][lt1p][:, lat1p, 2:-2]
    iondragdece = gdec['IonDragForce (east)'][lt1p][:, lat1p, 2:-2]
    iondragdec = np.sqrt(iondragdecn**2+iondragdece**2)
    iondragjunn = gjun['IonDragForce (north)'][lt1p][:, lat1p, 2:-2]
    iondragjune = gjun['IonDragForce (east)'][lt1p][:, lat1p, 2:-2]
    iondragjun = np.sqrt(iondragjunn**2+iondragjune**2)
    #vgradv
    vgradvdecn = (gdec['VGradV (north)']+gdec['SpheGeomForce (north)']+
            gdec['CentriForce (north)'])[lt1p][:, lat1p, 2:-2]
    vgradvdece = (gdec['VGradV (east)']+gdec['SpheGeomForce (east)'])\
            [lt1p][:, lat1p, 2:-2]
    vgradvdec = np.sqrt(vgradvdecn**2+vgradvdece**2)
    vgradvjunn = (gjun['VGradV (north)']+gjun['SpheGeomForce (north)']+
            gjun['CentriForce (north)'])[lt1p][:, lat1p, 2:-2]
    vgradvjune = (gjun['VGradV (east)']+gjun['SpheGeomForce (east)'])\
            [lt1p][:, lat1p, 2:-2]
    vgradvjun = np.sqrt(vgradvjunn**2+vgradvjune**2)
    #coriolis
    coriolisdecn = gdec['CoriolisForce (north)'][lt1p][:, lat1p, 2:-2]
    coriolisdece = gdec['CoriolisForce (east)'][lt1p][:, lat1p, 2:-2]
    coriolisdec = np.sqrt(coriolisdecn**2+coriolisdece**2)
    coriolisjunn = gjun['CoriolisForce (north)'][lt1p][:, lat1p, 2:-2]
    coriolisjune = gjun['CoriolisForce (east)'][lt1p][:, lat1p, 2:-2]
    coriolisjun = np.sqrt(coriolisjunn**2+coriolisjune**2)
    lat = lat[lt1p][:, lat1>=40, 2:-2]

    avevisdec = np.ones(lat.shape[2])
    avevisjun = np.ones(lat.shape[2])
    avepgdec = np.ones(lat.shape[2])
    avepgjun = np.ones(lat.shape[2])
    aveiondragdec = np.ones(lat.shape[2])
    aveiondragjun = np.ones(lat.shape[2])
    avevgradvdec = np.ones(lat.shape[2])
    avevgradvjun = np.ones(lat.shape[2])
    avecoriolisdec = np.ones(lat.shape[2])
    avecoriolisjun = np.ones(lat.shape[2])
    for ialt, valt in enumerate(altdec):
        avevisdec[ialt] = np.mean(visdec[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
        avevisjun[ialt] = np.mean(visjun[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                          np.mean(np.cos(lat[:,:,ialt]))
        avepgdec[ialt] = np.mean(pgdec[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
        avepgjun[ialt] = np.mean(pgjun[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                          np.mean(np.cos(lat[:,:,ialt]))
        aveiondragdec[ialt] = np.mean(iondragdec[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
        aveiondragjun[ialt] = np.mean(iondragjun[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                          np.mean(np.cos(lat[:,:,ialt]))
        avevgradvdec[ialt] = np.mean(vgradvdec[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
        avevgradvjun[ialt] = np.mean(vgradvjun[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                          np.mean(np.cos(lat[:,:,ialt]))
        avecoriolisdec[ialt] = np.mean(coriolisdec[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                         np.mean(np.cos(lat[:,:,ialt]))
        avecoriolisjun[ialt] = np.mean(coriolisjun[:,:,ialt]*np.cos(lat[:,:,ialt]))/ \
                          np.mean(np.cos(lat[:,:,ialt]))
    plt.close('all')
    plt.figure()
    plt.plot(avepgdec, altdec, 'ro-', label='Dec, PresGrad')
    #plt.plot(avepgjun, altjun, 'r*-', label='F107 = jun, PresGrad')
    plt.plot(avevisdec, altdec, 'bo-', label='Dec, Viscosity')
    #plt.plot(avevisjun, altjun, 'b*-', label='F107 = jun, Viscosity')
    #plt.plot(aveiondragdec, altdec, 'ko-', label='F107 = dec, Ion drag')
    #plt.plot(aveiondragjun, altjun, 'k*-', label='F107 = jun, Ion drag')
    plt.plot(avevgradvdec, altdec, 'yo-', label='Dec, VGradV')
    #plt.plot(avevgradvjun, altjun, 'y*-', label='F107 = jun, VGradV')
    plt.plot(avecoriolisdec, altdec, 'go-', label='Dec, Coriolis')
    #plt.plot(avecoriolisjun, altjun, 'g*-', label='F107 = jun, VGradV')

    #plt.plot(avepgdec, altdec, 'ro-', label='F107 = dec, PresGrad')
    plt.plot(avepgjun, altjun, 'r*--', label='Jun, PresGrad')
    #plt.plot(avevisdec, altdec, 'bo-', label='F107 = dec, Viscosity')
    plt.plot(avevisjun, altjun, 'b*--', label='Jun, Viscosity')
    #plt.plot(aveiondragdec, altdec, 'ko-', label='F107 = dec, Ion drag')
    #plt.plot(aveiondragjun, altjun, 'k*-', label='F107 = jun, Ion drag')
    #plt.plot(avevgradvdec, altdec, 'yo-', label='F107 = dec, VGradV')
    plt.plot(avevgradvjun, altjun, 'y*--', label='Jun, VGradV')
    #plt.plot(avecoriolisdec, altdec, 'go-', label='F107 = dec, VGradV')
    plt.plot(avecoriolisjun, altjun, 'g*--', label='Jun, Coriolis')
    plt.xlim([0,0.15])
    plt.ylim([0,700])
    plt.xlabel(r'Forces (m/$s^2$)')
    plt.ylabel('Altitude (km)')
    plt.legend()
    plt.show()
    plt.savefig('/home/guod/Desktop/why_no_low_density_cell_at_high_latitude/'\
                'allforces_dec_jun.pdf')
    return


def fig1(show=True):
    # rho, wind and temperature
    plt.close('all')
    plt.figure(figsize=(9, 6.5))
    filename = '/home/guod/simulation_output/momentum_analysis/'\
               'run_70_equinox/data/3DALL_t030322_000000.bin'
    path =  '/home/guod/Documents/work/fig/density_cell/'\
            'why_no_low_density_cell_at_high_latitude/'
    g = gitm.GitmBin(filename)
    gp.calc_pressure(g)
    # density and wind
    cbl = [np.linspace(1.1, 1.6, 21)*1e-10, np.linspace(0.5, 1.5, 21)*1e-12]
    tl = [np.arange(1.1, 1.7, 0.1)*1e-10, np.arange(0.5, 1.6, 0.2)*1e-12]
    for ialt, alt in enumerate([200, 400]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                2, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False, lonticklabel=[0, 0, 1, 1])
        lon0, lat0, zdata0 = g3ca.contour_data('Rho', g, alt=alt)
        hc = ax.contourf(
                lon0, lat0, zdata0, cbl[ialt], transform=ccrs.PlateCarree(),
                cmap='jet', extend='both')
        lon0, lat0, ewind, nwind = g3ca.vector_data(g, 'neutral', alt=alt)
        lon0, lat0, ewind, nwind = g3ca.convert_vector(
                lon0, lat0, ewind, nwind, plot_type='polar',
                projection=projection)
        if ialt==0:
            plt.text(-0.25, 0.5, 'Density and wind',transform=plt.gca().transAxes,
                     rotation=90, verticalalignment='center', fontsize=14)
        hq = ax.quiver(
                lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
                regrid_shape=20, headwidth=7)
        ax.quiverkey(hq, 0.1, -0.05, 1000, '1000 $m/s$')
        hc = plt.colorbar(hc, pad=0.1, ticks=tl[ialt])
        hc.set_label(r'$\rho$ (kg/m$^3$)')
        plt.title('%s km' % alt_str, fontsize=14)
    # pressure
    cbl = [np.linspace(3.5, 5, 21)*1e-5, np.linspace(0.2,1.2, 21)*1e-6]
    tl = [np.arange(3.5, 5.1, 0.5)*1e-5, np.arange(0.2, 1.3, 0.2)*1e-6]
    for ialt, alt in enumerate([200, 400]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                2, 2, ialt+3, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False, lonticklabel=[0, 0, 1, 1])
        lon0, lat0, zdata0 = g3ca.contour_data('pressure', g, alt=alt)
        hc = ax.contourf(
                lon0, lat0, zdata0, cbl[ialt], transform=ccrs.PlateCarree(),
                cmap='jet', extend='both')
        if ialt==0:
            plt.text(-0.25, 0.5, 'Pressure',transform=plt.gca().transAxes,
                     rotation=90, verticalalignment='center', fontsize=14)
        formatter = matplotlib.ticker.ScalarFormatter()
        formatter.set_powerlimits((0,0))
        hc = plt.colorbar(hc, pad=0.1, ticks=tl[ialt],format=formatter)
        hc.set_label('Pressure (Pa)')
        #ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('Pressure at %s km' % alt_str)
    plt.subplots_adjust(hspace=0.25)
    if show:
        plt.show()
    plt.savefig(path+'fig1.pdf')
    return


def fig2(show=True):
    # PresGrad, ion drag, coriolis, vgradv,viscosity
    plt.close('all')
    plt.figure(figsize=(5.44, 9.7))
    filename = '/home/guod/simulation_output/momentum_analysis/'\
               'run_70_equinox/data/3DALL_t030322_000000.bin'
    path =  '/home/guod/Documents/work/fig/density_cell/'\
            'why_no_low_density_cell_at_high_latitude/'
    g = gitm.GitmBin(filename)
    g['vgradvp (north)'] = g['VGradV (north)']+g['SpheGeomForce (north)']+\
                           g['CentriForce (north)']
    g['vgradvp (east)']= g['VGradV (east)']+g['SpheGeomForce (east)']
    nefl =[['NeuPressureGrad (north)', 'NeuPressureGrad (east)'],
           ['IonDragForce (north)', 'IonDragForce (east)'],
           ['CoriolisForce (north)', 'CoriolisForce (east)'],
           ['vgradvp (north)', 'vgradvp (east)'],
           ['ViscosityForce (north)', 'ViscosityForce (east)']]
    textl = [r'$-\frac{\nabla p}{\rho}$', r'$\nu_{ni}(\vec{V}_i-\vec{U}_n)$',
             r'$-2\vec{\omega}\times\vec{U}_n$',
             r'$\vec{U}_n\cdot\nabla\vec{U}_n$',
             r'$\frac{\partial}{\partial r}\eta\frac{\partial\vec{U}_n}{\partial r}$']
    for ialt, alt in enumerate([200, 400]):
        for iforce, force in enumerate(nefl):
            alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
            alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
            ax, projection = gcc.create_map(
                    5, 2, ialt+iforce*2+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                    centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                    coastlines=False, lonticklabel=[0,0,1,1])
            lon0 = np.array(g['dLon'][2:-2, 0, 0])
            lat0 = np.array(g['dLat'][0, 2:-2, 0])
            nf = np.array(g[force[0]][2:-2, 2:-2, alt_ind])
            ef = np.array(g[force[1]][2:-2, 2:-2, alt_ind])
            nf, lon0 = add_cyclic_point(nf.T, coord=lon0, axis=1)
            ef = add_cyclic_point(ef.T, axis=1)
            lon0, lat0 = np.meshgrid(lon0, lat0)
            lon0, lat0, ef, nf = g3ca.convert_vector(
                    lon0, lat0, ef, nf, plot_type='polar',
                    projection=projection)
            hq = ax.quiver(
                    lon0, lat0, ef, nf, scale=0.3, scale_units='inches',
                    regrid_shape=20, headwidth=8)
            ax.quiverkey(hq, 0.05, -0.1, 0.1, '0.1 $m/s^2$')
            if iforce==0:
                plt.title('%s km' % alt_str, y=1.05, fontsize=14)
            if ialt==0:
                plt.text(-0.35, 0.50, textl[iforce], fontsize=14,
                         verticalalignment='center', rotation=90,
                         transform=plt.gca().transAxes)

    plt.subplots_adjust(bottom=0.05, top =0.95)
    # plt.text(0.5, 0.95, 'SpheGeom', fontsize=15, horizontalalignment='center',
    #          transform=plt.gcf().transFigure)
    plt.show()
    plt.savefig(path+'fig2.pdf')
    return


def fig3(show=True):
    # rho, wind and temperature
    plt.close('all')
    plt.figure(figsize=(9.36, 5.23))
    filename = '/home/guod/simulation_output/momentum_analysis/'\
               'run_150_equinox/data/3DALL_t030322_000001.bin'
    path =  '/home/guod/Documents/work/fig/density_cell/'\
            'why_no_low_density_cell_at_high_latitude/'
    g = gitm.GitmBin(filename)
    gp.calc_pressure(g)
    # density and wind
    cbl = [np.linspace(2.1, 2.7, 21)*1e-10, np.linspace(6, 10, 21)*1e-12,
           np.linspace(4, 9, 21)*1e-13]
    tl = [np.arange(2.1, 2.8, 0.1)*1e-10, np.arange(6, 11, 1)*1e-12,
          np.arange(4,10,1)*1e-13]
    for ialt, alt in enumerate([200, 400, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                2, 3, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False, lonticklabel=[0, 0, 1, 1])
        lon0, lat0, zdata0 = g3ca.contour_data('Rho', g, alt=alt)
        hc = ax.contourf(
                lon0, lat0, zdata0, cbl[ialt], transform=ccrs.PlateCarree(),
                cmap='jet', extend='both')
        lon0, lat0, ewind, nwind = g3ca.vector_data(g, 'neutral', alt=alt)
        lon0, lat0, ewind, nwind = g3ca.convert_vector(
                lon0, lat0, ewind, nwind, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.1, -0.05, 1000, '1000 $m/s$')
        if ialt==0:
            plt.text(-0.25, 0.5, 'Density and wind',transform=plt.gca().transAxes,
                     rotation=90, verticalalignment='center', fontsize=14)
        hc = plt.colorbar(hc, pad=0.1, ticks=tl[ialt])
        hc.set_label(r'$\rho$ (kg/m$^3$)')
        plt.title('%s km' % alt_str, fontsize=14)
    # pressure
    cbl = [np.linspace(9, 13, 21)*1e-5, np.linspace(4,9, 21)*1e-6,
           np.linspace(30, 80)*1e-8]
    tl = [np.arange(9, 14, 1)*1e-5, np.arange(4, 10, 1)*1e-6,
          np.arange(30, 90, 10)*1e-8]
    for ialt, alt in enumerate([200, 400, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                2, 3, ialt+4, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False, lonticklabel=[0, 0, 1, 1])
        lon0, lat0, zdata0 = g3ca.contour_data('pressure', g, alt=alt)
        hc = ax.contourf(
                lon0, lat0, zdata0, cbl[ialt], transform=ccrs.PlateCarree(),
                cmap='jet', extend='both')
        if ialt==0:
            plt.text(-0.25, 0.5, 'Pressure',transform=plt.gca().transAxes,
                     rotation=90, verticalalignment='center', fontsize=14)
        formatter = matplotlib.ticker.ScalarFormatter()
        formatter.set_powerlimits((0,0))
        hc = plt.colorbar(hc, pad=0.1, ticks=tl[ialt],format=formatter)
        hc.set_label('Pressure (Pa)')
        #ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
    plt.subplots_adjust(left=0.06, bottom=0.05, right=0.95)
    if show:
        plt.show()
    plt.savefig(path+'fig3.pdf')
    return


def fig4(show=True):
    # PresGrad, ion drag, coriolis, vgradv,viscosity
    plt.close('all')
    plt.figure(figsize=(7.13, 9.46))
    filename = '/home/guod/simulation_output/momentum_analysis/'\
               'run_150_equinox/data/3DALL_t030322_000001.bin'
    path =  '/home/guod/Documents/work/fig/density_cell/'\
            'why_no_low_density_cell_at_high_latitude/'
    g = gitm.GitmBin(filename)
    g['vgradvp (north)'] = g['VGradV (north)']+g['SpheGeomForce (north)']+\
                           g['CentriForce (north)']
    g['vgradvp (east)']= g['VGradV (east)']+g['SpheGeomForce (east)']
    nefl =[['NeuPressureGrad (north)', 'NeuPressureGrad (east)'],
           ['IonDragForce (north)', 'IonDragForce (east)'],
           ['CoriolisForce (north)', 'CoriolisForce (east)'],
           ['vgradvp (north)', 'vgradvp (east)'],
           ['ViscosityForce (north)', 'ViscosityForce (east)']]
    textl = [r'$-\frac{\nabla p}{\rho}$', r'$\nu_{ni}(\vec{V}_i-\vec{U}_n)$',
             r'$-2\vec{\omega}\times\vec{U}_n$',
             r'$\vec{U}_n\cdot\nabla\vec{U}_n$',
             r'$\frac{\partial}{\partial r}\eta\frac{\partial\vec{U}_n}{\partial r}$']
    for ialt, alt in enumerate([200, 400, 600]):
        for iforce, force in enumerate(nefl):
            alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
            alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
            ax, projection = gcc.create_map(
                    5, 3, ialt+iforce*3+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                    centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                    coastlines=False, lonticklabel=[0,0,1,1])
            lon0 = np.array(g['dLon'][2:-2, 0, 0])
            lat0 = np.array(g['dLat'][0, 2:-2, 0])
            nf = np.array(g[force[0]][2:-2, 2:-2, alt_ind])
            ef = np.array(g[force[1]][2:-2, 2:-2, alt_ind])
            nf, lon0 = add_cyclic_point(nf.T, coord=lon0, axis=1)
            ef = add_cyclic_point(ef.T, axis=1)
            lon0, lat0 = np.meshgrid(lon0, lat0)
            lon0, lat0, ef, nf = g3ca.convert_vector(
                    lon0, lat0, ef, nf, plot_type='polar',
                    projection=projection)
            hq = ax.quiver(
                    lon0, lat0, ef, nf, scale=0.3, scale_units='inches',
                    regrid_shape=20, headwidth=8)
            ax.quiverkey(hq, 0.05, -0.1, 0.1, '0.1 $m/s^2$')
            if iforce==0:
                plt.title('%s km' % alt_str, y=1.05, fontsize=14)
            if ialt==0:
                plt.text(-0.35, 0.50, textl[iforce], fontsize=14,
                         verticalalignment='center', rotation=90,
                         transform=plt.gca().transAxes)
    plt.subplots_adjust(bottom=0.05, top =0.95)
    # plt.text(0.5, 0.95, 'SpheGeom', fontsize=15, horizontalalignment='center',
    #          transform=plt.gcf().transFigure)
    plt.show()
    plt.savefig(path+'fig4.pdf')
    return


def plot_animation_den_win(show=False):
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
        lon1, lat1, zdata1 = g3ca.contour_data('Rho', g1, alt=400)
        lon2, lat2, zdata2 = g3ca.contour_data('Rho', g2, alt=400)
        hc = [ax[k].contourf(
                lon1, lat1, zdata1, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(6e-12, 10e-12, 21),
                #levels=np.linspace(15e-13, 25e-13, 21),
                cmap='jet', extend='both') for k in [0, 1]]
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = [ax[k].contourf(
                lon2, lat2, zdata2, 21,transform=ccrs.PlateCarree(),
                levels=np.linspace(6e-12, 10e-12, 21),
                #levels=np.linspace(30e-13, 60e-13, 21),
                cmap='jet', extend='both') for k in [2, 3]]
        diffzdata = 100*(zdata2-zdata1)/zdata1
        hc = [ax[k].contourf(
                lon2, lat2, diffzdata, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-30, 30, 21), cmap='seismic',
                extend='both') for k in [4, 5]]

        # wind
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=400)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2, 'neutral', alt=400)
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
            fig, animate_den_wind, interval=200, frames=len(fn1))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'den_wind.mp4',writer=writer)
    return


def plot_animation_all_forces(show=True):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_2/data/'
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
        lon0, lat0, eallf, nallf = g3ca.convert_vector(
                lon0, lat0, eallf, nallf, plot_type='polar',
                projection=projection[0])
        hq = ax[0].quiver(
                lon0, lat0, eallf, nallf, scale=0.3, scale_units='inches',
                regrid_shape=20)
        ax[0].quiverkey(hq, 0.93, -0.1, 0.2, '0.2 $m/s^2$')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hq = ax[1].quiver(
                lon0, lat0, eallf, nallf, scale=0.3, scale_units='inches',
                regrid_shape=20)
        ax[1].quiverkey(hq, 0.93, -0.1, 0.2, '0.2 $m/s^2$')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    anim = animation.FuncAnimation(
            fig, animate_all_forces, interval=200, frames=len(fn1))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'all_forces.mp4',writer=writer)
    return


def plot_animation_density_change(show=True):
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
        # density change (shrink)
        lon1 = np.array(g1['Longitude'])
        lat1 = np.array(g1['Latitude'])
        alt1 = np.array(g1['Altitude'])
        Re = 6371*1000 # Earth radius, unit: m
        RR = Re+alt1
        omega = 2*np.pi/24
        rho1 = np.array(g1['Rho'])
        nwind1 = np.array(g1['V!Dn!N (north)'])
        ewind1 = np.array(g1['V!Dn!N (east)'])# + omega*RR*np.cos(lat1)
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
        omega = 2*np.pi/24
        rho2 = np.array(g2['Rho'])
        nwind2 = np.array(g2['V!Dn!N (north)'])
        ewind2 = np.array(g2['V!Dn!N (east)'])# + omega*RR*np.cos(lat2)
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
                lon1, lat1, div_rhov_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon1, lat1, div_rhov_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    anim = animation.FuncAnimation(
            fig, animate_density_change, interval=200, frames=len(fn1))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'density_change.mp4',writer=writer)
    return


def plot_animation_rhodivv(show=True):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 00:30:00')
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
        # density change (shrink)
        lon1 = np.array(g1['Longitude'])
        lat1 = np.array(g1['Latitude'])
        alt1 = np.array(g1['Altitude'])
        Re = 6371*1000 # Earth radius, unit: m
        RR = Re+alt1
        omega = 2*np.pi/24
        rho1 = np.array(g1['Rho'])
        nwind1 = np.array(g1['V!Dn!N (north)'])
        ewind1 = np.array(g1['V!Dn!N (east)'])# + omega*RR*np.cos(lat1)
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
        omega = 2*np.pi/24
        rho2 = np.array(g2['Rho'])
        nwind2 = np.array(g2['V!Dn!N (north)'])
        ewind2 = np.array(g2['V!Dn!N (east)'])# + omega*RR*np.cos(lat2)
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
                lon1, lat1, div_rhov_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[0].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
        hc = ax[1].contourf(
                lon1, lat1, div_rhov_d, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(-5e-16, 5e-16, 21),
                cmap='seismic', extend='both')
        ax[1].set_title(g1['time'].strftime('%d-%b-%y %H:%M')+' UT', y=1.05)
    anim = animation.FuncAnimation(
            fig, animate_density_change, frames=len(fn1), repeat=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    #anim.save(path+'density_change.mp4',writer=writer)
    return


def plot_animation_den_win_1(show=False):
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_2/data/'
    fn1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]

    # save path
    path = '/home/guod/Documents/work/fig/density_cell/' \
           + 'why_no_low_density_cell_at_high_latitude/ion_drag_with_or_not/'
    fig = plt.figure(figsize=[8,6])
    # read gitm data
    def animate_den_wind(i):
        g1 = gitm.GitmBin(fn1[i], varlist=[
            'Rho', 'V!Dn!N (north)', 'V!Dn!N (east)', 'V!Dn!N (up)',
            'V!Di!N (east)', 'V!Di!N (north)'])
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
        # Density
        lon1, lat1, zdata1 = g3ca.contour_data('Rho', g1, alt=400)
        hc = [ax[k].contourf(
                lon1, lat1, zdata1, 21, transform=ccrs.PlateCarree(),
                levels=np.linspace(6e-12, 10e-12, 21),
                cmap='jet', extend='both') for k in [0, 1]]

        # wind
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1, 'neutral', alt=400)
        for ins in [0, 1]:
            lon0, lat0, ewind, nwind = g3ca.convert_vector(
                    lon1, lat1, ewind1, nwind1, plot_type='polar',
                    projection=projection[ins])
            hq = ax[ins].quiver(
                    lon0, lat0, ewind, nwind, scale=1500, scale_units='inches',
                    color='k', regrid_shape=20)
            ax[ins].quiverkey(hq, 0.93, 0, 1000, '1000 m/s')
        return
    anim = animation.FuncAnimation(
            fig, animate_den_wind, interval=200, frames=len(fn1))
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save(path+'den_wind.mp4',writer=writer)
    return


def plot_all_figures():
    plot_den_win(show=False)
    plot_pressure(show=False)
    plot_temperature(show=False)
    plot_ave_m(show=False)
    plot_pressure_gradient(show=False)
    plot_ion_drag(show=False)
    plot_coriolis(show=False)
    plot_vgradv(show=False)
    plot_viscosity(show=False)
    plot_all_forces(show=False)
    return


if __name__=='__main__':
    plot_animation_density_change()
