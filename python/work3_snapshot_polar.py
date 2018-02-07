#Global imports
import numpy as np
import matplotlib.pyplot as plt
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

def plot_den_win(show=True, save=True):
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    g = g2a
    for ialt, alt in enumerate([130, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, 2:-2]/1000-alt))+2
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
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

        # rho difference
        # lon3, lat3, zdata3 = g3ca.contour_data('Rho', g1a, alt=alt)
        # lon4, lat4, zdata4 = g3ca.contour_data('Rho', g2a, alt=alt)
        # diffrho = 100*(zdata4-zdata3)/zdata3
        # hc = ax.contour(
        #         lon4, lat4, diffrho, [-10], transform=ccrs.PlateCarree(),
        #         colors='g',linestyles='-')
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'01_den_win_run2_%s%s.pdf' %(tstrday,tstring))
    return


def plot_den_win_diff(show=True, save=True):
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([130, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
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
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 500, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        # rho difference
        #lon3, lat3, zdata3 = g3ca.contour_data('Rho', g1a, alt=alt)
        #lon4, lat4, zdata4 = g3ca.contour_data('Rho', g2a, alt=alt)
        #diffrho = 100*(zdata4-zdata3)/zdata3
        #hc = ax.contour(
        #        lon4, lat4, diffrho, [-10], transform=ccrs.PlateCarree(),
        #        colors='g',linestyles='-')
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
        horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'01_den_win_diff_%s%s.pdf' %(tstrday,tstring))
    return


def plot_temperature(show=True, save=True):
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    g = g2a
    for ialt, alt in enumerate([130, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, 2:-2]/1000-alt))+2
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('Temperature', g, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(
            lon0, lat0, zdata0, 21,
            transform=ccrs.PlateCarree(), cmap='jet', extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label('Temperature (K)')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
        horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'02_temperature_run2_%s%s.pdf' % (tstrday,tstring))
    return


def plot_temperature_diff(show=True, save=True):
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([130, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
            3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
            coastlines=False)
        # temperature diff
        lon1, lat1, zdata1 = g3ca.contour_data('Temperature', g1a, alt=alt)
        lon2, lat2, zdata2 = g3ca.contour_data('Temperature', g2a, alt=alt)
        fp = (lat1[:,0]>slat) & (lat1[:,0]<nlat)
        lon0, lat0, zdata0 = lon1[fp, :], lat1[fp,:], (zdata2-zdata1)[fp,:]
        hc = ax.contourf(
            lon0, lat0, zdata0, np.linspace(-80,80,21),
            transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$T_2-T_1$ (K)')
        # geomagnetic pole
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())

        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
            lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
            projection=projection)
        hq = ax.quiver(
            lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
            regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 500, '500 m/s')

        # rho difference
        #lon3, lat3, zdata3 = g3ca.contour_data('Rho', g1a, alt=alt)
        #lon4, lat4, zdata4 = g3ca.contour_data('Rho', g2a, alt=alt)
        #diffrho = 100*(zdata4-zdata3)/zdata3
        #hc = ax.contour(
        #        lon4, lat4, diffrho, [-10], transform=ccrs.PlateCarree(),
        #        colors='g',linestyles='-')
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'02_temperature_diff_%s%s.pdf' % (tstrday,tstring))
    return


def plot_ion_drift(show=True,save=True):
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)
    g = g2a

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([130, 200, 300, 400, 500, 600]):
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
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'03_ion_drift_%s%s.pdf' % (tstrday,tstring))
    return


def plot_vert_divv_diff(show=True, save=True):
    velr = np.array(g1a['V!Dn!N (up)'])
    divv1 = cr.calc_div_vert(g1a['Altitude'], velr)

    velr = np.array(g2a['V!Dn!N (up)'])
    divv2 = cr.calc_div_vert(g2a['Altitude'], velr)

    g1a['vert_divv_diff'] = divv1-divv2

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([130, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
            3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
            coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('vert_divv_diff', g1a, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(
            lon0, lat0, zdata0, levels=np.linspace(-1,1,21)*1e-4,
            transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        hcb = plt.colorbar(hc, pad=0.17)
        hcb.set_label(r'$-\nabla\cdot\vec{u}$ (up)')
        hcb.formatter.set_powerlimits((-2,2))
        hcb.update_ticks()

        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
            lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
            projection=projection)
        hq = ax.quiver(
            lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
            regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 500, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'04_vert_divv_diff_%s%s.pdf' % (tstrday,tstring))
    return


def plot_hozt_divv_diff(show=True, save=True):

    lat1 = np.array(g1a['Latitude'])
    alt1 = np.array(g1a['Altitude'])
    Re = 6371*1000 # Earth radius, unit: m
    RR = Re+alt1
    omega = 2*np.pi/(24*3600)
    nwind1 = np.array(g1a['V!Dn!N (north)'])
    ewind1 = np.array(g1a['V!Dn!N (east)'])# + omega*RR*np.cos(lat1)
    divv1 = cr.calc_div_hozt(
        g1a['Longitude'], g1a['Latitude'], g1a['Altitude'], nwind1, ewind1)

    lat2 = np.array(g2a['Latitude'])
    alt2 = np.array(g2a['Altitude'])
    Re = 6371*1000 # Earth radius, unit: m
    RR = Re+alt2
    omega = 2*np.pi/(24*3600)
    nwind2 = np.array(g2a['V!Dn!N (north)'])
    ewind2 = np.array(g2a['V!Dn!N (east)']) #+ omega*RR*np.cos(lat2)
    divv2 = cr.calc_div_hozt(
        g2a['Longitude'], g2a['Latitude'], g2a['Altitude'], nwind2, ewind2)

    g1a['hozt_divv_diff'] = divv1-divv2

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([130, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
            3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
            coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('hozt_divv_diff', g1a, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(
            lon0, lat0, zdata0, levels=np.linspace(-1,1,21)*1e-4,
            transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        hcb = plt.colorbar(hc, pad=0.17)
        hcb.formatter.set_powerlimits((0,0))
        hcb.update_ticks()
        hcb.set_label(r'$-\nabla\cdot\vec{u}$ (horizontal)')
        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
            lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
            projection=projection)
        hq = ax.quiver(
            lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
            regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 500, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'05_hozt_divv_diff_%s%s.pdf' % (tstrday,tstring))
    return


def plot_vert_vgradrho_rho_diff(show=True, save=True):
    rho1 = np.array(g1a['Rho'])
    vgradrho1 = \
        g1a['V!Dn!N (up)'] \
        * cr.calc_rusanov_alts_ausm(g1a['Altitude'],rho1)/g1a['Rho']

    rho2 = np.array(g2a['Rho'])
    vgradrho2 = \
        g2a['V!Dn!N (up)']\
        * cr.calc_rusanov_alts_ausm(g2a['Altitude'],rho2)/g2a['Rho']

    g1a['vgradrho_diff'] = vgradrho1-vgradrho2

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([130, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
            3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
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
        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
            lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
            projection=projection)
        hq = ax.quiver(
            lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
            regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 500, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'06_vert_vgradrho_rho_diff_%s%s.pdf' % (tstrday,tstring))
    return


def plot_hozt_vgradrho_rho_diff(show=True, save=True):
    lon1 = np.array(g1a['Longitude'])
    lat1 = np.array(g1a['Latitude'])
    alt1 = np.array(g1a['Altitude'])
    Re = 6371*1000 # Earth radius, unit: m
    RR = Re+alt1
    omega = 2*np.pi/(24*3600)
    rho1 = np.array(g1a['Rho'])
    nwind1 = np.array(g1a['V!Dn!N (north)'])
    ewind1 = np.array(g1a['V!Dn!N (east)']) + omega*RR*np.cos(lat1)
    vgradrho1 = \
        (nwind1*cr.calc_rusanov_lats(lat1,alt1,rho1)\
        +ewind1*cr.calc_rusanov_lons(lon1,lat1,alt1,rho1))/rho1

    lon2 = np.array(g2a['Longitude'])
    lat2 = np.array(g2a['Latitude'])
    alt2 = np.array(g2a['Altitude'])
    Re = 6371*1000 # Earth radius, unit: m
    RR = Re+alt2
    omega = 2*np.pi/(24*3600)
    rho2 = np.array(g2a['Rho'])
    nwind2 = np.array(g2a['V!Dn!N (north)'])
    ewind2 = np.array(g2a['V!Dn!N (east)']) + omega*RR*np.cos(lat2)
    vgradrho2 = \
        (nwind2*cr.calc_rusanov_lats(lat2,alt2,rho2)\
         +ewind2*cr.calc_rusanov_lons(lon2,lat2,alt2,rho2))/rho2
    g1a['hozt_vgradrho_diff'] = vgradrho1-vgradrho2

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([130, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
            3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
            coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('hozt_vgradrho_diff', g1a, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(
            lon0, lat0, zdata0, np.linspace(-1,1,21)*1e-4,
            transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.formatter.set_powerlimits((0,0))
        hc.update_ticks()
        hc.set_label(r'-$\vec{u}\cdot\frac{\nabla\rho}{\rho}$ (horizontal)')
        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 500, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'07_hozt_vgradrho_rho_diff_%s%s.pdf' % (tstrday,tstring))
    return


def plot_divrhov_rho_diff(show=True, save=True):
    # density change (shrink)
    lon1 = np.array(g1a['Longitude'])
    lat1 = np.array(g1a['Latitude'])
    alt1 = np.array(g1a['Altitude'])
    Re = 6371*1000 # Earth radius, unit: m
    RR = Re+alt1
    omega = 2*np.pi/(24*3600)
    rho1 = np.array(g1a['Rho'])
    nwind1 = np.array(g1a['V!Dn!N (north)'])
    ewind1 = np.array(g1a['V!Dn!N (east)']) + omega*RR*np.cos(lat1)
    uwind1 = np.array(g1a['V!Dn!N (up)'])
    div_rhov1 = \
        (cr.calc_div_hozt(lon1, lat1, alt1, rho1*nwind1, rho1*ewind1)\
        +cr.calc_div_vert(alt1, rho1*uwind1))/rho1

    # density change (no shrink)
    lon2 = np.array(g2a['Longitude'])
    lat2 = np.array(g2a['Latitude'])
    alt2 = np.array(g2a['Altitude'])
    Re = 6371*1000 # Earth radius, unit: m
    RR = Re+alt2
    omega = 2*np.pi/(24*3600)
    rho2 = np.array(g2a['Rho'])
    nwind2 = np.array(g2a['V!Dn!N (north)'])
    ewind2 = np.array(g2a['V!Dn!N (east)']) + omega*RR*np.cos(lat2)
    uwind2 = np.array(g2a['V!Dn!N (up)'])
    div_rhov2 = \
        (cr.calc_div_hozt(lon2, lat2, alt2, rho2*nwind2, rho2*ewind2)\
        +cr.calc_div_vert(alt2, rho2*uwind2))/rho2
    g1a['divrhov_diff'] = div_rhov1-div_rhov2

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([130, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
            3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
            centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
            coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('divrhov_diff', g1a, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(
            lon0, lat0, zdata0, np.linspace(-1,1,21)*1e-4,
            transform=ccrs.PlateCarree(), cmap='seismic', extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.formatter.set_powerlimits((0,0))
        hc.update_ticks()
        hc.set_label(r'$-\frac{\nabla\cdot(\rho\vec{u})}{\rho}$')
        # wind difference
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 500, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'08_divrhov_rho_diff_%s%s.pdf' % (tstrday,tstring))
    return


if __name__=='__main__':
    import gc
    import gitm
    import gitm_create_coordinate as gcc
    Re = 6371*1000 # Earth radius, unit: m
    trange = pd.date_range(
            '2003-03-22 00:00:00', '2003-03-22 04:00:00', freq='10min')
    spath = '/home/guod/Documents/work/fig/density_cell/'\
            'why_no_low_density_cell_at_high_latitude/snapshot_polar/'\
            'snapshot_polar_01/'
    den_change_label = \
        [np.linspace(-3.5, 3.5,21)*1e-13,np.linspace(-1.2,1.2,21)*1e-14,
         np.linspace(-1.6, 1.6,21)*1e-15,np.linspace(-5,5,21)*1e-16,
         np.linspace(-2,2,21)*1e-16, np.linspace(-8,8,21)*1e-17]
    gradp_label = \
        [np.linspace(-0.035, 0.035,21),np.linspace(-0.03,0.03,21),
         np.linspace(-0.05, 0.05,21),np.linspace(-0.05,0.05,21),
         np.linspace(-0.05,0.05,21), np.linspace(-0.05,0.05,21)]
    for t in trange:
        tstring = t.strftime('%H%M')
        tstrday = t.strftime('%d')
        filename = glob.glob(
            '/home/guod/simulation_output/momentum_analysis/'\
            'run_shrink_iondrift_4_c1'\
            '/data/3DALL_t0303%s_%s*.bin' % (tstrday,tstring))[0]
        g1a = gitm.GitmBin(filename)
        filename = glob.glob(
            '/home/guod/simulation_output/momentum_analysis/'\
            'run_no_shrink_iondrift_4_1'\
            '/data/3DALL_t0303%s_%s*.bin' % (tstrday,tstring))[0]
        g2a = gitm.GitmBin(filename)

        nlat, slat = -40, -90

        #plot_den_win(show=False)
        #plot_den_win_diff(show=False)
        #plot_temperature(show=False)
        #plot_temperature_diff(show=False)
        #plot_ion_drift(show=False)
        #plot_vert_divv_diff(show=False)
        #plot_hozt_divv_diff(show=False)
        #plot_vert_vgradrho_rho_diff(show=False)
        #plot_hozt_vgradrho_rho_diff(show=False)
        plot_divrhov_rho_diff(show=False)
    gc.collect()
