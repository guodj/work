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

def plot_den_win(show=True, save=True):
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
        lon0, lat0, zdata0 = g3ca.contour_data('Rho', g, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, 21,
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
        plt.savefig(spath+'den_win.pdf')
    return


def plot_den_win_diff(show=True, save=True):
    apex = Apex(date=2010)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        lon1, lat1, zdata1 = g3ca.contour_data('Rho', g1a, alt=alt)
        lon2, lat2, zdata2 = g3ca.contour_data('Rho', g2a, alt=alt)
        fp = (lat1[:,0]>slat) & (lat1[:,0]<nlat)
        lon0, lat0, zdata0 = lon1[fp, :], lat1[fp,:], \
                             (100*(zdata2-zdata1)/zdata1)[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, np.linspace(-20,20,21),
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 300, '500 m/s')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$100*\frac{\rho2-\rho1}{\rho1}$ (%)')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'01_den_win_diff_%s.pdf' %tstring)
    return


def plot_ion_drift(show=True,save=True):
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
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'ion_drift.pdf')
    return


def plot_temperature(show=True, save=True):
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
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'temperature.pdf')
    return


def plot_temperature_diff(show=True, save=True):
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        lon1, lat1, zdata1 = g3ca.contour_data('Temperature', g1a, alt=alt)
        lon2, lat2, zdata2 = g3ca.contour_data('Temperature', g2a, alt=alt)
        fp = (lat1[:,0]>slat) & (lat1[:,0]<nlat)
        lon0, lat0, zdata0 = lon1[fp, :], lat1[fp,:], \
                             (100*(zdata2-zdata1)/zdata1)[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, np.linspace(-20,20,21),
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$100*\frac{T_2-T_1}{T_1}$ (%)')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'02_temperature_diff_%s.pdf' % tstring)
    return


def plot_vertical_wind(show=True, save=True):
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
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, np.linspace(-10, 10, 21),
                         transform=ccrs.PlateCarree(), cmap='jet',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label('Vertical Wind (m/s)')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'vertical_wind.pdf')
    return


def plot_pressure(show=True, save=True):
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
        lon0, lat0, zdata0 = g3ca.contour_data('pressure', g, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, 21,
                         transform=ccrs.PlateCarree(), cmap='jet',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17, format='%.1e')
        hc.set_label('pressure (Pa)')
        #ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'pressure.pdf')
    return


def plot_pressure_diff(show=True, save=True):
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        lon1, lat1, zdata1 = g3ca.contour_data('pressure', g1a, alt=alt)
        lon2, lat2, zdata2 = g3ca.contour_data('pressure', g2a, alt=alt)
        fp = (lat1[:,0]>slat) & (lat1[:,0]<nlat)
        lon0, lat0, zdata0 = lon1[fp, :], lat1[fp,:], \
                             (100*(zdata2-zdata1)/zdata1)[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, np.linspace(-20,20,21),
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$100*\frac{P_2-P_1}{P_1}$ (%)')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'03_pressure_diff_%s.pdf' % tstring)
    return


def plot_ave_m(show=True, save=True):
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
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('ave_m', g, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0,21,
                         transform=ccrs.PlateCarree(), cmap='jet',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label('Ave M')
        #ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'ave_m.pdf')
    return


def plot_ave_m_diff(show=True, save=True):
    Re = 6371*1000 # Earth radius, unit: m
    ave_m1 = ((g1a['O(!U3!NP)']*16 + g1a['O!D2!N']*32 + g1a['N!D2!N']*28 +
              g1a['N(!U4!NS)']*14 + g1a['NO']*30 + g1a['N(!U2!ND)']*14 +
              g1a['N(!U2!NP)']*14 + g1a['H']*2 + g1a['He']*4 + g1a['CO!D2!N']*44 +
              g1a['O(!U1!ND)']*16) /
             (g1a['O(!U3!NP)']+g1a['O!D2!N']+g1a['N!D2!N']+g1a['N(!U4!NS)']+g1a['NO']+
              g1a['N(!U2!ND)']+g1a['N(!U2!NP)']+g1a['H']+g1a['He']+g1a['CO!D2!N']+
              g1a['O(!U1!ND)']))
    ave_m2 = ((g2a['O(!U3!NP)']*16 + g2a['O!D2!N']*32 + g2a['N!D2!N']*28 +
              g2a['N(!U4!NS)']*14 + g2a['NO']*30 + g2a['N(!U2!ND)']*14 +
              g2a['N(!U2!NP)']*14 + g2a['H']*2 + g2a['He']*4 + g2a['CO!D2!N']*44 +
              g2a['O(!U1!ND)']*16) /
             (g2a['O(!U3!NP)']+g2a['O!D2!N']+g2a['N!D2!N']+g2a['N(!U4!NS)']+g2a['NO']+
              g2a['N(!U2!ND)']+g2a['N(!U2!NP)']+g2a['H']+g2a['He']+g2a['CO!D2!N']+
              g2a['O(!U1!ND)']))
    g1a['ave_m_diff'] = dmarray(
            100*(ave_m2-ave_m1)/ave_m1,
            attrs={'units':'', 'scale':'linear', 'name':'ave_m'})

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('ave_m_diff', g1a, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, np.linspace(-20,20,21),
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$100*\frac{M_2-M_1}{M_1}$ (%)')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'04_ave_m_diff_%s.pdf' % tstring)
    return


def plot_vert_rhodivv_diff(show=True, save=True):
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
    rhodivv1 = rho1*(1.0/(RR**2) * np.gradient((RR**2)*uwind1, axis=2)\
                     / np.gradient(alt1, axis=2))

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
    rhodivv2 = rho2*(1.0/(RR**2) * np.gradient((RR**2)*uwind2, axis=2) \
                     / np.gradient(alt2, axis=2))
    g1a['vert_rhodivv_diff'] = rhodivv1-rhodivv2

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('vert_rhodivv_diff', g1a, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, levels=den_change_label[ialt],
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$\rho\nabla\cdot\vec{u}$ (up)')
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 300, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'05_vert_rhodivv_diff_%s.pdf' % tstring)
    return


def plot_hozt_rhodivv_diff(show=True, save=True):
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
    rhodivv1 = rho1 * (1.0/(RR*np.cos(lat1))
            * (np.gradient(nwind1*np.cos(lat1), axis=1)
            / np.gradient(lat1, axis=1)
            + np.gradient(ewind1, axis=0) / np.gradient(lon1, axis=0)))

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
    rhodivv2 = rho2 * (1.0/(RR*np.cos(lat2))
            * (np.gradient(nwind2*np.cos(lat2), axis=1)
            / np.gradient(lat2, axis=1)
            + np.gradient(ewind2, axis=0) / np.gradient(lon2, axis=0)))
    g1a['hozt_rhodivv_diff'] = rhodivv1-rhodivv2

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_ind = np.min([alt_ind,g1a['Altitude'].shape[2]-3])
        alt_ind = np.max([alt_ind,2])
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('hozt_rhodivv_diff', g1a, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, levels=den_change_label[ialt],
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$\rho\nabla\cdot\vec{u}$ (horizontal)')
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 300, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'07_hozt_rhodivv_diff_%s.pdf' % tstring)
    return


def plot_vert_vgradrho_diff(show=True, save=True):
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
    vgradrho1 = uwind1 * (np.gradient(rho1, axis=2)
                / np.gradient(alt1, axis=2))

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
    vgradrho2 = uwind2 * (np.gradient(rho2, axis=2)
                / np.gradient(alt2, axis=2))
    g1a['vgradrho_diff'] = vgradrho1-vgradrho2

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('vgradrho_diff', g1a, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, den_change_label[ialt],
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$\vec{u}\cdot\nabla\rho$ (up)')
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 300, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'06_vert_vgradrho_diff_%s.pdf' % tstring)
    return


def plot_hozt_vgradrho_diff(show=True, save=True):
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
    vgradrho1 = nwind1 * ((1.0/RR)*np.gradient(rho1, axis=1)
                          / np.gradient(lat1, axis=1))\
              + ewind1 * ((1.0/(RR*np.cos(lat1)))*np.gradient(rho1, axis=0)
                          / np.gradient(lon1, axis=0))

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
    vgradrho2 = nwind2 * ((1.0/RR)*np.gradient(rho2, axis=1)
                          / np.gradient(lat2, axis=1))\
              + ewind2 * ((1.0/(RR*np.cos(lat2)))*np.gradient(rho2, axis=0)
                          / np.gradient(lon2, axis=0))
    g1a['hozt_vgradrho_diff'] = vgradrho1-vgradrho2

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('hozt_vgradrho_diff', g1a, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, den_change_label[ialt],
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$\vec{u}\cdot\nabla\rho$ (horizontal)')
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 300, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'08_hozt_vgradrho_diff_%s.pdf' % tstring)
    return


def plot_divrhov_diff(show=True, save=True):
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
    div_rhov1 = (
            1.0/(RR**2)
          * np.gradient((RR**2)*rho1*uwind1, axis=2) / np.gradient(alt1, axis=2)
          + 1.0/(RR*np.cos(lat1))
          * (np.gradient(rho1*nwind1*np.cos(lat1), axis=1) / np.gradient(lat1, axis=1)
             + np.gradient(rho1*ewind1, axis=0) / np.gradient(lon1, axis=0)))
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
    div_rhov2 = (
            1.0/(RR**2)
          * np.gradient((RR**2)*rho2*uwind2, axis=2) / np.gradient(alt2, axis=2)
          + 1.0/(RR*np.cos(lat2))
          * (np.gradient(rho2*nwind2*np.cos(lat2), axis=1) / np.gradient(lat2, axis=1)
             + np.gradient(rho2*ewind2, axis=0) / np.gradient(lon2, axis=0)))
    g1a['divrhov_diff'] = div_rhov1-div_rhov2

    apex = Apex(date=2003)
    qlat, qlon = apex.convert(-90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        lon0, lat0, zdata0 = g3ca.contour_data('divrhov_diff', g1a, alt=alt)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, zdata0 = lon0[fp, :], lat0[fp,:], zdata0[fp,:]
        hc = ax.contourf(lon0, lat0, zdata0, den_change_label[ialt],
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        hc.set_label(r'$\nabla\cdot(\rho\vec{u})$')
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 300, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'09_divrhov_diff_%s.pdf' % tstring)
    return


def plot_vert_pressure_gradient(show=True, save=True):
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        lon0 = np.array(g1a['dLon'][2:-2, 0, 0])
        lat0 = np.array(g1a['dLat'][0, 2:-2, 0])
        gradp_up1 = np.array(g1m['NeuPressureGrad (up)'][2:-2,2:-2,alt_ind])
        gradp_up2 = np.array(g2m['NeuPressureGrad (up)'][2:-2,2:-2,alt_ind])
        gradp_diff = gradp_up2-gradp_up1
        gradp_diff, lon0 = add_cyclic_point(gradp_diff.T, coord=lon0, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, gradp_diff = lon0[fp, :], lat0[fp,:], gradp_diff[fp,:]
        hc = ax.contourf(lon0, lat0, gradp_diff, levels=gradp_label[ialt],
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 300, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    #plt.text(0.5, 0.95, 'Pressure Gradient Force', fontsize=15,
    #        horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'10_vert_gradp_diff_%s.pdf' % tstring)
    return


def plot_vert_forces(show=True, save=True):
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

    plt.close('all')
    plt.figure(figsize=(7.26, 9.25))
    for ialt, alt in enumerate([120, 200, 300, 400, 500, 600]):
        alt_ind = np.argmin(np.abs(g1a['Altitude'][0, 0, :]/1000-alt))
        alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
        ax, projection = gcc.create_map(
                3, 2, ialt+1, 'polar', nlat=nlat, slat=slat, dlat=10,
                centrallon=g3ca.calculate_centrallon(g1a, 'polar',  useLT=True),
                coastlines=False)
        lon0 = np.array(g1a['dLon'][2:-2, 0, 0])
        lat0 = np.array(g1a['dLat'][0, 2:-2, 0])
        force_up1 = np.array(g1m['NeuPressureGrad (up)']
                            +g1m['CoriolisForce (up)']
                            +g1m['CentriForce (up)']
                            +g1m['VGradV (up)']
                            +g1m['IonDragForce (up)']
                            +g1m['SpheGeomForce (up)'])[2:-2,2:-2,alt_ind]
        force_up2 = np.array(g2m['NeuPressureGrad (up)']
                            +g2m['CoriolisForce (up)']
                            +g2m['CentriForce (up)']
                            +g2m['VGradV (up)']
                            +g2m['IonDragForce (up)']
                            +g2m['SpheGeomForce (up)'])[2:-2,2:-2,alt_ind]
        force_up = force_up2-force_up1
        force_up, lon0 = add_cyclic_point(force_up.T, coord=lon0, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        fp = (lat0[:,0]>slat) & (lat0[:,0]<nlat)
        lon0, lat0, force_up = lon0[fp, :], lat0[fp,:], force_up[fp,:]
        hc = ax.contourf(lon0, lat0, force_up, levels=gradp_label[ialt],
                         transform=ccrs.PlateCarree(), cmap='seismic',
                         extend='both')
        hc = plt.colorbar(hc, pad=0.17)
        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1a, 'neutral', alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2a, 'neutral', alt=alt)
        lon0, lat0 = lon1, lat1
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
                lon0, lat0, ewind2-ewind1, nwind2-nwind1, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ewind0, nwind0, scale=1000, scale_units='inches',
                regrid_shape=20, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 300, '500 m/s')
        ax.scatter(qlon, qlat, color='k', transform=ccrs.PlateCarree())
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    #plt.text(0.5, 0.95, 'Pressure Gradient Force', fontsize=15,
    #        horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'11_vert_force_diff_%s.pdf' % tstring)
    return


def plot_pressure_gradient(show=True, save=True):
    filename ='/home/guod/simulation_output/momentum_analysis/'\
              'run_no_shrink_iondrift_3/data/3DMOM_t030322_003502.bin'
    g = gitm.GitmBin(filename,varlist=[
        'NeuPressureGrad (east)','NeuPressureGrad (north)',
        'NeuPressureGrad (up)'])
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
        plt.savefig(spath+'pressure_grad.pdf')
    return


def plot_coriolis(show=True, save=True):
    global g
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
        plt.savefig(spath+'coriolis.pdf')
    return


def plot_ion_drag(show=True, save=True):
    global g
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
        plt.savefig(spath+'ion_drag.pdf')
    return


def plot_viscosity(show=True, save=True):
    global g
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
        plt.savefig(spath+'viscosity.pdf')
    return


def plot_vgradv(show=True, save=True):
    global g
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
        lon0 = np.array(g['dLon'][2:-2, 0, 0])
        lat0 = np.array(g['dLat'][0, 2:-2, 0])
        nf = np.array((g['VGradV (north)']+g['SpheGeomForce (north)'])
                [2:-2, 2:-2, alt_ind])
        ef = np.array((g['VGradV (east)']+g['SpheGeomForce (east)'])
                [2:-2, 2:-2, alt_ind])
        #nf = np.array((g['CentriForce (north)'])[2:-2, 2:-2, alt_ind])
        #ef = np.array((g['SpheGeomForce (east)'])[2:-2, 2:-2, alt_ind])*0
        nf, lon0 = add_cyclic_point(nf.T, coord=lon0, axis=1)
        ef = add_cyclic_point(ef.T, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        lon0, lat0, ef, nf = g3ca.convert_vector(
                lon0, lat0, ef, nf, plot_type='polar',
                projection=projection)
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
        plt.savefig(spath+'vgradv.pdf')
    return


def plot_centriforce(show=True, save=True):
    global g
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
        lon0 = np.array(g['dLon'][2:-2, 0, 0])
        lat0 = np.array(g['dLat'][0, 2:-2, 0])
        nf = np.array((g['CentriForce (north)'])[2:-2, 2:-2, alt_ind])
        ef = np.array((g['SpheGeomForce (east)'])[2:-2, 2:-2, alt_ind])*0
        nf, lon0 = add_cyclic_point(nf.T, coord=lon0, axis=1)
        ef = add_cyclic_point(ef.T, axis=1)
        lon0, lat0 = np.meshgrid(lon0, lat0)
        lon0, lat0, ef, nf = g3ca.convert_vector(
                lon0, lat0, ef, nf, plot_type='polar',
                projection=projection)
        hq = ax.quiver(
                lon0, lat0, ef, nf, scale=0.1, scale_units='inches',
                regrid_shape=15, headwidth=5)
        ax.quiverkey(hq, 0.93, -0.1, 0.2, '0.2 $m/s^2$')
        plt.title('%s km' % alt_str, y=1.05)
    plt.text(0.5, 0.95, 'Centrifugal', fontsize=15, horizontalalignment='center',
             transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'centrifugal.pdf')
    return


def plot_all_forces(show=True, save=True):
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
        plt.savefig(spath+'all_forces.pdf')
    return


if __name__=='__main__':
    import gc
    from PyPDF2 import PdfFileMerger, PdfFileReader
    Re = 6371*1000 # Earth radius, unit: m

    trange = pd.date_range(
            '2003-03-22 00:00:00', '2003-03-22 04:00:00', freq='10min')
    spath = '/home/guod/Documents/work/fig/density_cell/'\
            'why_no_low_density_cell_at_high_latitude/snapshot/'
    outpdf_den_win = PdfFileMerger()
    outpdf_temperature = PdfFileMerger()
    outpdf_pressure = PdfFileMerger()
    outpdf_ave_m = PdfFileMerger()
    outpdf_vert_rhodivv = PdfFileMerger()
    outpdf_hozt_rhodivv = PdfFileMerger()
    outpdf_vert_vgradrho = PdfFileMerger()
    outpdf_hozt_vgradrho = PdfFileMerger()
    outpdf_divrhov = PdfFileMerger()
    outpdf_vert_gradp = PdfFileMerger()
    outpdf_vert_force = PdfFileMerger()
    den_change_label = [np.linspace(-3.5, 3.5,21)*1e-13,np.linspace(-1.2,1.2,21)*1e-14,
                        np.linspace(-1.6, 1.6,21)*1e-15,np.linspace(-5,5,21)*1e-16,
                        np.linspace(-2,2,21)*1e-16, np.linspace(-8,8,21)*1e-17]
    gradp_label = [np.linspace(-0.035, 0.035,21),np.linspace(-0.03,0.03,21),
                   np.linspace(-0.05, 0.05,21),np.linspace(-0.05,0.05,21),
                   np.linspace(-0.05,0.05,21), np.linspace(-0.05,0.05,21)]
    for t in trange:
        tstring = t.strftime('%H%M')
        filename = glob.glob(
                '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_3_continue'\
                '/data/3DALL_t030322_%s*.bin' % tstring)[0]
        g1a = gitm.GitmBin(filename)
        filename = glob.glob(
                '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_3_continue'\
                '/data/3DMOM_t030322_%s*.bin' % tstring)[0]
        g1m = gitm.GitmBin(filename)
        filename = glob.glob(
                '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_3'\
                '/data/3DALL_t030322_%s*.bin' % tstring)[0]
        g2a = gitm.GitmBin(filename)
        filename = glob.glob(
                '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_3'\
                '/data/3DMOM_t030322_%s*.bin' % tstring)[0]
        g2m = gitm.GitmBin(filename)

        gp.calc_pressure(g1a)
        gp.calc_pressure(g2a)
        nlat, slat = -40, -90

        # plot_den_win_diff(show=False)
        # outpdf_den_win.append(PdfFileReader(open(spath+'01_den_win_diff_%s.pdf' % tstring, 'rb')))
        # plot_temperature_diff(show=False)
        # outpdf_temperature.append(PdfFileReader(open(spath+'02_temperature_diff_%s.pdf' % tstring, 'rb')))
        # plot_pressure_diff(show=False)
        # outpdf_pressure.append(PdfFileReader(open(spath+'03_pressure_diff_%s.pdf' % tstring, 'rb')))
        # plot_ave_m_diff(show=False)
        # outpdf_ave_m.append(PdfFileReader(open(spath+'04_ave_m_diff_%s.pdf' % tstring, 'rb')))
        # plot_vert_rhodivv_diff(show=False)
        # outpdf_vert_rhodivv.append(PdfFileReader(open(spath+'05_vert_rhodivv_diff_%s.pdf' % tstring, 'rb')))
        # plot_hozt_rhodivv_diff(show=False)
        # outpdf_hozt_rhodivv.append(PdfFileReader(open(spath+'07_hozt_rhodivv_diff_%s.pdf' % tstring, 'rb')))
        # plot_vert_vgradrho_diff(show=False)
        # outpdf_vert_vgradrho.append(PdfFileReader(open(spath+'06_vert_vgradrho_diff_%s.pdf' % tstring, 'rb')))
        # plot_hozt_vgradrho_diff(show=False)
        # outpdf_hozt_vgradrho.append(PdfFileReader(open(spath+'08_hozt_vgradrho_diff_%s.pdf' % tstring, 'rb')))
        # plot_divrhov_diff(show=False)
        # outpdf_divrhov.append(PdfFileReader(open(spath+'09_divrhov_diff_%s.pdf' % tstring, 'rb')))
        # plot_vert_pressure_gradient(show=False)
        # outpdf_vert_gradp.append(PdfFileReader(open(spath+'10_vert_gradp_diff_%s.pdf' % tstring, 'rb')))
        plot_vert_forces(show=False)
        outpdf_vert_force.append(PdfFileReader(open(spath+'11_vert_force_diff_%s.pdf' % tstring, 'rb')))

    #outpdf_den_win.write(spath+'01_den_win_diff.pdf')
    #outpdf_temperature.write(spath+'02_temperature_diff.pdf')
    #outpdf_pressure.write(spath+'03_pressure_diff.pdf')
    #outpdf_ave_m.write(spath+'04_ave_m_diff.pdf')
    #outpdf_vert_rhodivv.write(spath+'05_vert_rhodivv_diff.pdf')
    #outpdf_hozt_rhodivv.write(spath+'07_hozt_rhodivv_diff.pdf')
    #outpdf_vert_vgradrho.write(spath+'06_vert_vgradrho_diff.pdf')
    #outpdf_divrhov.write(spath+'09_divrhov_diff.pdf')
    #outpdf_vert_gradp.write(spath+'10_vert_gradp_diff.pdf')
    outpdf_vert_force.write(spath+'11_vert_force.pdf')
    gc.collect()
