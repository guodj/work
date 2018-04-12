import gitm_create_coordinate as gcc
import matplotlib.pyplot as plt
from aacgmv2 import convert
import datetime as dt
import gitm_3D_const_alt as g3ca
import numpy as np
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import matplotlib.pyplot as plt
import gitm_new as gitm
from matplotlib.ticker import ScalarFormatter
plt.style.use('ggplot')


def calc_rhoi(g):
    g['rhoi'] = 1.66e-27*(
        g['O!D2!U+!N']*32 + g['O_4SP_!U+!N']*16 + g['N!U+!N']*14 +
        g['N!D2!U+!N']*28 + g['NO!U+!N']*30 + g['He!U+!N']*4 +
        g['H!U+!N']*1 + g['O(!U2!NP)!U+!N']*16 +g['O(!U2!ND)!U+!N']*16)
    return g['rhoi']

def plot_polar_contour_vector(
    nrow, ncol, iax, g, nlat, slat, alt, contour=True, zstr='Rho',
    vector=True, neuion='neu', useLT=True, coastlines=False, N=21, levels=None):

    centrallon = g3ca.calculate_centrallon(g, 'pol', useLT)
    ax, projection = gcc.create_map(
            nrow, ncol, iax, 'pol', nlat, slat, centrallon,
            coastlines=coastlines, dlat=10,
            useLT=useLT, lonticklabel=(0, 0, 1, 1))
    hc = []
    hq = []
    # contour
    if contour:
        lon0, lat0, zdata0 = g3ca.contour_data(zstr, g, alt=alt)
        cb = N if levels is None else levels
        hc.append(ax.contourf(lon0, lat0, zdata0, cb,
            transform=ccrs.PlateCarree(),cmap='jet', extend='both'))
    # vector
    if vector:
        lon0, lat0, ewind0, nwind0 = g3ca.vector_data(g,neuion,alt=alt)
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
            lon0, lat0, ewind0, nwind0, 'pol', projection)
        hq.append(ax.quiver(
            lon0,lat0,ewind0,nwind0,scale=1500,scale_units='inches',
            regrid_shape=20))
    return ax, projection, hc, hq

def plot_polar_contour_vector_diff(
    nrow, ncol, iax, g1, g2, nlat, slat, alt, contour=True, zstr='Rho',
    vector=True, neuion='neu', useLT=True, coastlines=False, N=21, levels=None):
    # g1 and g2 have the same time, the same grid

    centrallon = g3ca.calculate_centrallon(g1, 'pol', useLT)
    ax, projection = gcc.create_map(
            nrow, ncol, iax, 'pol', nlat, slat, centrallon,
            coastlines=coastlines, dlat=10,
            useLT=useLT, lonticklabel=(0, 0, 1, 1))
    hc = []
    hq = []
    # contour
    if contour:
        lon0, lat0, zdata1 = g3ca.contour_data(zstr, g1, alt=alt)
        lon0, lat0, zdata2 = g3ca.contour_data(zstr, g2, alt=alt)
        cb = N if levels is None else levels
        hc.append(ax.contourf(lon0, lat0, 100*(zdata2-zdata1)/zdata1, cb,
            transform=ccrs.PlateCarree(),cmap='jet', extend='both'))
    # vector
    if vector:
        lon0, lat0, ewind1, nwind1 = g3ca.vector_data(g1,neuion,alt=alt)
        lon0, lat0, ewind2, nwind2 = g3ca.vector_data(g2,neuion,alt=alt)
        lon0, lat0, ewind0, nwind0 = g3ca.convert_vector(
            lon0, lat0, ewind2-ewind1, nwind2-nwind1, 'pol', projection)
        hq.append(ax.quiver(
            lon0,lat0,ewind0,nwind0,scale=1500,scale_units='inches',
            regrid_shape=20))
    return ax, projection, hc, hq

def plot_06ut_18ut():
    path = '/home/guod/simulation_output/momentum_analysis/run_shrink_dipole_1_c1/UA/data'
    fn1 = path+'/3DALL_t030322_060002.bin'
    fn2 = path+'/3DALL_t030322_180002.bin'
    g11 = [gitm.read(fn1), gitm.read(fn2)] # shrink, dipole

    path = '/home/guod/simulation_output/momentum_analysis/run_shrink_igrf_1_c1/UA/data'
    fn1 = path+'/3DALL_t030322_060002.bin'
    fn2 = path+'/3DALL_t030322_180003.bin'
    g12 = [gitm.read(fn1), gitm.read(fn2)] # shrink, igrf

    path = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_dipole_1_1/UA/data'
    fn1 = path+'/3DALL_t030322_060000.bin'
    fn2 = path+'/3DALL_t030322_180001.bin'
    g21 = [gitm.read(fn1), gitm.read(fn2)] # no shrink, diple

    path = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_igrf_1_1/UA/data'
    fn1 = path+'/3DALL_t030322_060001.bin'
    fn2 = path+'/3DALL_t030322_180002.bin'
    g22 = [gitm.read(fn1), gitm.read(fn2)] # no shrink, igrf

    #  # Rho and wind
    #  plt.figure(1, figsize=(8, 7.55))
    #  nlat, slat, alt = -30, -90, 300
    #  #nlat, slat, alt = 90, 30, 300
    #  glats, glons = convert(-90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  glatn, glonn = convert(90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  for k in range(4):
    #      gtemp = g22 if k in [0, 1] else g21
    #      ax, projection, hc, hq = plot_polar_contour_vector(
    #          2, 2, k+1, gtemp[k%2], nlat, slat, alt, contour=True, zstr='Rho',
    #          vector=True, neuion='neu', useLT=True, coastlines=False,
    #          levels=np.linspace(1e-11, 4e-11, 21))
    #      if k in [0, 1]:
    #          ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
    #          ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
    #      if k in [0, 1]:
    #          plt.title('UT: '+gtemp[k%2]['time'].strftime('%H%M'), y=1.03)
    #      if k in [0, 2]:
    #          text = 'IGRF' if k==0 else 'Dipole'
    #          plt.text(
    #              -0.2, 0.5, text, fontsize=14, verticalalignment='center',
    #              transform=plt.gca().transAxes, rotation=90)
    #  plt.subplots_adjust(bottom=0.13,top=0.94, wspace=0.1, hspace=0.1)
    #  cax = plt.axes([0.3, 0.08, 0.4, 0.02])
    #  plt.colorbar(
    #      hc[0], cax=cax, orientation='horizontal',
    #      ticks=np.arange(1e-11, 4.1e-11, 1e-11))
    #  cax.set_xlabel(r'$\rho (kg/m^3)$')

    #  # Difference Rho and wind
    #  plt.figure(2, figsize=(8, 7.55))
    #  nlat, slat, alt = -30, -90, 300
    #  #nlat, slat, alt = 90, 30, 300
    #  glats, glons = convert(-90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  glatn, glonn = convert( 90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  for k in range(4):
    #      gtemp1 = g12 if k in [0, 1] else g11
    #      gtemp2 = g22 if k in [0, 1] else g21
    #      ax, projection, hc, hq = plot_polar_contour_vector_diff(
    #          2, 2, k+1, gtemp1[k%2], gtemp2[k%2], nlat, slat, alt, contour=True,
    #          zstr='Rho', vector=True, neuion='neu', useLT=True,
    #          coastlines=False, levels=np.linspace(-30, 30, 21))
    #      if k in [0, 1]:
    #          ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
    #          ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
    #      if k in [0, 1]:
    #          plt.title('UT: '+gtemp1[k%2]['time'].strftime('%H%M'), y=1.03)
    #      if k in [0, 2]:
    #          text = 'IGRF' if k==0 else 'Dipole'
    #          plt.text(
    #              -0.2, 0.5, text, fontsize=14, verticalalignment='center',
    #              transform=plt.gca().transAxes, rotation=90)
    #  plt.subplots_adjust(bottom=0.13,top=0.94, wspace=0.1, hspace=0.1)
    #  cax = plt.axes([0.3, 0.08, 0.4, 0.02])
    #  plt.colorbar(
    #      hc[0], cax=cax, orientation='horizontal',
    #      ticks=np.arange(-30, 31, 10))
    #  cax.set_xlabel(r'$\rho_d (\%)$')

    # ion drag
    #  plt.figure(1, figsize=(8, 7.55))
    #  nlat, slat, alt = -30, -90, 300
    #  nlat, slat, alt = 90, 30, 300
    #  glats, glons = convert(-90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  glatn, glonn = convert(90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  for k in range(4):
    #      gtemp = g22 if k in [0, 1] else g21
    #      calc_rhoi(gtemp[0])
    #      calc_rhoi(gtemp[1])
    #      gtemp[0]['iondrag'] = \
    #          gtemp[0]['rhoi']/gtemp[0]['Rho']*gtemp[0]['Collision(in)']*\
    #          np.sqrt((gtemp[0]['V!Di!N (east)']-gtemp[0]['V!Dn!N (east)'])**2 +\
    #                  (gtemp[0]['V!Di!N (north)']-gtemp[0]['V!Dn!N (north)'])**2)
    #      gtemp[1]['iondrag'] = \
    #          gtemp[1]['rhoi']/gtemp[1]['Rho']*gtemp[1]['Collision(in)']*\
    #          np.sqrt((gtemp[1]['V!Di!N (east)']-gtemp[1]['V!Dn!N (east)'])**2 +\
    #                  (gtemp[1]['V!Di!N (north)']-gtemp[1]['V!Dn!N (north)'])**2)
    #      ax, projection, hc, hq = plot_polar_contour_vector(
    #          2, 2, k+1, gtemp[k%2], nlat, slat, alt, contour=True, zstr='iondrag',
    #          vector=False, neuion='neu', useLT=True, coastlines=False,
    #          levels=np.linspace(0, 0.2, 21))
    #      if k in [0, 1]:
    #          ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
    #          ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
    #      if k in [0, 1]:
    #          plt.title('UT: '+gtemp[k%2]['time'].strftime('%H%M'), y=1.03)
    #      if k in [0, 2]:
    #          text = 'IGRF' if k==0 else 'Dipole'
    #          plt.text(
    #              -0.2, 0.5, text, fontsize=14, verticalalignment='center',
    #              transform=plt.gca().transAxes, rotation=90)
    #  plt.subplots_adjust(bottom=0.13,top=0.94, wspace=0.1, hspace=0.1)
    #  cax = plt.axes([0.3, 0.08, 0.4, 0.02])
    #  plt.colorbar(
    #      hc[0], cax=cax, orientation='horizontal',
    #      ticks=np.arange(0, 0.21, 0.05))
    #  cax.set_xlabel(r'Ion Drag $(m/s^2)$')

    #  # vni
    #  plt.figure(1, figsize=(8, 7.55))
    #  nlat, slat, alt = -30, -90, 300
    #  nlat, slat, alt = 90, 30, 300
    #  glats, glons = convert(-90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  glatn, glonn = convert(90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  for k in range(4):
    #      gtemp = g22 if k in [0, 1] else g21
    #      calc_rhoi(gtemp[0])
    #      calc_rhoi(gtemp[1])
    #      gtemp[0]['vni'] = \
    #          gtemp[0]['rhoi']/gtemp[0]['Rho']*gtemp[0]['Collision(in)']
    #      gtemp[1]['vni'] = \
    #          gtemp[1]['rhoi']/gtemp[1]['Rho']*gtemp[1]['Collision(in)']
    #      ax, projection, hc, hq = plot_polar_contour_vector(
    #          2, 2, k+1, gtemp[k%2], nlat, slat, alt, contour=True, zstr='vni',
    #          vector=False, neuion='neu', useLT=True, coastlines=False,
    #          levels=np.linspace(0, 0.001, 21))
    #      if k in [0, 1]:
    #          ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
    #          ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
    #      if k in [0, 1]:
    #          plt.title('UT: '+gtemp[k%2]['time'].strftime('%H%M'), y=1.03)
    #      if k in [0, 2]:
    #          text = 'IGRF' if k==0 else 'Dipole'
    #          plt.text(
    #              -0.2, 0.5, text, fontsize=14, verticalalignment='center',
    #              transform=plt.gca().transAxes, rotation=90)
    #  plt.subplots_adjust(bottom=0.13,top=0.94, wspace=0.1, hspace=0.1)
    #  cax = plt.axes([0.3, 0.08, 0.4, 0.02])
    #  plt.colorbar(
    #      hc[0], cax=cax, orientation='horizontal',
    #      ticks=np.arange(0, 0.0012, 0.0002))
    #  cax.set_xlabel(r'$\nu_{ni}$ $(s^{-1})$')

    #  # ion density
    #  plt.figure(1, figsize=(8, 7.55))
    #  nlat, slat, alt = -30, -90, 300
    #  #nlat, slat, alt = 90, 30, 300
    #  glats, glons = convert(-90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  glatn, glonn = convert(90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  for k in range(4):
    #      gtemp = g22 if k in [0, 1] else g21
    #      calc_rhoi(gtemp[0])
    #      calc_rhoi(gtemp[1])
    #      ax, projection, hc, hq = plot_polar_contour_vector(
    #          2, 2, k+1, gtemp[k%2], nlat, slat, alt, contour=True, zstr='rhoi',
    #          vector=False, neuion='neu', useLT=True, coastlines=False,
    #          levels=np.linspace(0, 3e-14, 21))
    #      if k in [0, 1]:
    #          ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
    #          ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
    #      if k in [0, 1]:
    #          plt.title('UT: '+gtemp[k%2]['time'].strftime('%H%M'), y=1.03)
    #      if k in [0, 2]:
    #          text = 'IGRF' if k==0 else 'Dipole'
    #          plt.text(
    #              -0.2, 0.5, text, fontsize=14, verticalalignment='center',
    #              transform=plt.gca().transAxes, rotation=90)
    #  plt.subplots_adjust(bottom=0.13,top=0.94, wspace=0.1, hspace=0.1)
    #  cax = plt.axes([0.3, 0.08, 0.4, 0.02])
    #  plt.colorbar(
    #      hc[0], cax=cax, orientation='horizontal',
    #      ticks=np.arange(0, 3.1e-14, 5e-15))
    #  cax.set_xlabel(r'$\rho_i$ $(kg/m^3)$')

    #  #  electron density
    #  plt.figure(1, figsize=(8, 7.55))
    #  nlat, slat, alt = -30, -90, 300
    #  nlat, slat, alt = 90, 30, 300
    #  glats, glons = convert(-90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  glatn, glonn = convert(90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    #  for k in range(4):
    #      gtemp = g22 if k in [0, 1] else g21
    #      calc_rhoi(gtemp[0])
    #      calc_rhoi(gtemp[1])
    #      ax, projection, hc, hq = plot_polar_contour_vector(
    #          2, 2, k+1, gtemp[k%2], nlat, slat, alt, contour=True, zstr='e-',
    #          vector=False, neuion='neu', useLT=True, coastlines=False,
    #          levels=np.linspace(0, 1.5e12, 21))
    #      if k in [0, 1]:
    #          ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
    #          ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
    #      if k in [0, 1]:
    #          plt.title('UT: '+gtemp[k%2]['time'].strftime('%H%M'), y=1.03)
    #      if k in [0, 2]:
    #          text = 'IGRF' if k==0 else 'Dipole'
    #          plt.text(
    #              -0.2, 0.5, text, fontsize=14, verticalalignment='center',
    #              transform=plt.gca().transAxes, rotation=90)
    #  plt.subplots_adjust(bottom=0.13,top=0.94, wspace=0.1, hspace=0.1)
    #  cax = plt.axes([0.3, 0.08, 0.4, 0.02])
    #  plt.colorbar(
    #      hc[0], cax=cax, orientation='horizontal',
    #      ticks=np.arange(0, 1.6e12, 0.3e12))
    #  cax.set_xlabel(r'Ne ($m^{-3}$)')

    # ion convection
    plt.figure(1, figsize=(8, 7.55))
    nlat, slat, alt = -30, -90, 300
    nlat, slat, alt = 90, 30, 300
    glats, glons = convert(-90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    glatn, glonn = convert(90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    for k in range(4):
        gtemp = g22 if k in [0, 1] else g21
        gtemp[0]['vi'] = np.sqrt(
            gtemp[0]['V!Di!N (east)']**2 + gtemp[0]['V!Di!N (north)']**2)
        gtemp[1]['vi'] = np.sqrt(
            gtemp[1]['V!Di!N (east)']**2 + gtemp[1]['V!Di!N (north)']**2)
        ax, projection, hc, hq = plot_polar_contour_vector(
            2, 2, k+1, gtemp[k%2], nlat, slat, alt, contour=True, zstr='vi',
            vector=True, neuion='ion', useLT=True, coastlines=False,
            levels=np.linspace(0, 1000, 21))
        if k in [0, 1]:
            ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
            ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
        if k in [0, 1]:
            plt.title('UT: '+gtemp[k%2]['time'].strftime('%H%M'), y=1.03)
        if k in [0, 2]:
            text = 'IGRF' if k==0 else 'Dipole'
            plt.text(
                -0.2, 0.5, text, fontsize=14, verticalalignment='center',
                transform=plt.gca().transAxes, rotation=90)
    plt.subplots_adjust(bottom=0.13,top=0.94, wspace=0.1, hspace=0.1)
    cax = plt.axes([0.3, 0.08, 0.4, 0.02])
    plt.colorbar(
        hc[0], cax=cax, orientation='horizontal',
        ticks=np.arange(0, 1001, 200))
    cax.set_xlabel(r'Vi (m/s)')
    return

def plot_rho_wind_ut():
    path = 'D:/simulation_results/low_density_cell/run_shrink_igrf_1_c1/UA/data'
    fn01 = path+'/3DALL_t030323_000000.bin'
    fn02 = path+'/3DALL_t030323_030002.bin'
    fn03 = path+'/3DALL_t030323_060003.bin'
    fn04 = path+'/3DALL_t030323_090002.bin'
    fn05 = path+'/3DALL_t030323_120000.bin'
    fn06 = path+'/3DALL_t030323_150000.bin'
    fn07 = path+'/3DALL_t030323_180002.bin'
    fn08 = path+'/3DALL_t030323_210001.bin'
    fn09 = path+'/3DALL_t030324_000000.bin'
    fn1 = [fn01, fn02, fn03, fn04, fn05, fn06, fn07, fn08, fn09]
    path = 'D:/simulation_results/low_density_cell/run_no_shrink_igrf_1_1/UA/data'
    fn01 = path+'/3DALL_t030323_000000.bin'
    fn02 = path+'/3DALL_t030323_030000.bin'
    fn03 = path+'/3DALL_t030323_060001.bin'
    fn04 = path+'/3DALL_t030323_090002.bin'
    fn05 = path+'/3DALL_t030323_120000.bin'
    fn06 = path+'/3DALL_t030323_150001.bin'
    fn07 = path+'/3DALL_t030323_180000.bin'
    fn08 = path+'/3DALL_t030323_210002.bin'
    fn09 = path+'/3DALL_t030324_000000.bin'
    fn2 = [fn01, fn02, fn03, fn04, fn05, fn06, fn07, fn08, fn09]

    #nlat, slat = -30, -90
    nlat, slat = 90, 30
    alt = 300
    glats, glons = convert(-90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    glatn, glonn = convert(90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    fig1 = plt.figure(1, figsize=(8,7.75))
    fig2 = plt.figure(2, figsize=(8,7.75))
    fig2 = plt.figure(3, figsize=(8,7.75))
    for k in range(9):
        g1, g2 = gitm.read(fn1[k]), gitm.read(fn2[k])

        title = 'UT: '+g1['time'].strftime('%H%M')
        lon1, lat1, rho1 = g3ca.contour_data('Rho', g1, alt=alt)
        lon2, lat2, rho2 = g3ca.contour_data('Rho', g2, alt=alt)
        lon0, lat0, rho0 = lon1, lat1, 100*(rho2-rho1)/rho1

        lon1, lat1, ewind1, nwind1 = g3ca.vector_data(g1,'neu',alt=alt)
        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2,'neu',alt=alt)

        if k in [0, 3]:
            lonticklabel = [0, 0, 0, 1]
        elif k in [7, 8]:
            lonticklabel = [0, 0, 1, 0]
        elif k==6:
            lonticklabel = [0, 0, 1, 1]
        else:
            lonticklabel = [0, 0, 0, 0]

        plt.figure(1)
        centrallon = g3ca.calculate_centrallon(g1, 'pol', useLT=True)
        ax, projection = gcc.create_map(
            3,3,k+1, 'pol', nlat, slat, centrallon, coastlines=False,
            dlat=10, useLT=True, lonticklabel=lonticklabel)
        ax.contourf(lon2, lat2, rho2, np.linspace(1e-11, 4e-11),
            transform=ccrs.PlateCarree(),cmap='jet', extend='both')
        lon9, lat9, ewind9, nwind9 = g3ca.convert_vector(
            lon2, lat2, ewind2, nwind2, 'pol', projection)
        ax.quiver(
            lon9,lat9,ewind9,nwind9,scale=1500,scale_units='inches',
            regrid_shape=20)
        ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=35, c='k')
        ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=35, c='k')
        plt.title(title)

        plt.figure(2)
        centrallon = g3ca.calculate_centrallon(g1, 'pol', useLT=True)
        ax, projection = gcc.create_map(
            3,3,k+1, 'pol', nlat, slat, centrallon, coastlines=False,
            dlat=10, useLT=True, lonticklabel=lonticklabel)
        ax.contourf(lon1, lat1, rho1, np.linspace(1e-11, 4e-11),
            transform=ccrs.PlateCarree(),cmap='jet', extend='both')
        lon9, lat9, ewind9, nwind9 = g3ca.convert_vector(
            lon1, lat1, ewind1, nwind1, 'pol', projection)
        ax.quiver(
            lon9,lat9,ewind9,nwind9,scale=1500,scale_units='inches',
            regrid_shape=20)
        ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=35, c='k')
        ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=35, c='k')
        plt.title(title)

        plt.figure(3)
        centrallon = g3ca.calculate_centrallon(g1, 'pol', useLT=True)
        ax, projection = gcc.create_map(
            3,3,k+1, 'pol', nlat, slat, centrallon, coastlines=False,
            dlat=10, useLT=True, lonticklabel=lonticklabel)
        ax.contourf(lon0, lat0, rho0, np.linspace(-30, 30, 21),
            transform=ccrs.PlateCarree(),cmap='jet', extend='both')
        lon9, lat9, ewind9, nwind9 = g3ca.convert_vector(
            lon2, lat2, ewind2-ewind1, nwind2-nwind1, 'pol', projection)
        ax.quiver(
            lon9,lat9,ewind9,nwind9,scale=1500,scale_units='inches',
            regrid_shape=20)
        ax.scatter(glons, glats, transform=ccrs.Geodetic(), s=35, c='k')
        ax.scatter(glonn, glatn, transform=ccrs.Geodetic(), s=35, c='k')
        #if k==2:
        #    ut1 = g1['time'].hour+g1['time'].minute/60
        #    onlat = np.linspace(-90,-30, 30)
        #    onlon1 = np.ones(onlat.shape)*(8-ut1)*15
        #    onlon2 = np.ones(onlat.shape)*(20-ut1)*15
        #    ax.scatter(onlon1, onlat, s=5, c='y', transform=ccrs.PlateCarree(),
        #               linewidths=0)
        #    ax.scatter(onlon2, onlat, s=5, c='y', transform=ccrs.PlateCarree(),
        #               linewidths=0)
        #if k==6:
        #    ut1 = g1['time'].hour+g1['time'].minute/60
        #    onlat = np.linspace(-90,-30, 30)
        #    onlon1 = np.ones(onlat.shape)*(4-ut1)*15
        #    onlon2 = np.ones(onlat.shape)*(16-ut1)*15
        #    ax.scatter(onlon1, onlat, s=5, c='y', transform=ccrs.PlateCarree(),
        #               linewidths=0)
        #    ax.scatter(onlon2, onlat, s=5, c='y', transform=ccrs.PlateCarree(),
        #               linewidths=0)
        plt.title(title)

    plt.figure(1)
    plt.subplots_adjust(
        hspace=0.2, wspace=0, left=0.05, right=0.95, top=0.95, bottom=0.05)
    plt.figure(2)
    plt.subplots_adjust(
        hspace=0.2, wspace=0, left=0.05, right=0.95, top=0.95, bottom=0.05)
    plt.figure(3)
    plt.subplots_adjust(
        hspace=0.2, wspace=0, left=0.05, right=0.95, top=0.95, bottom=0.05)
    return

def plot_one_meridian():
    path = 'D:/simulation_results/low_density_cell/run_shrink_igrf_1_c1/UA/data'
    fn01 = path+'/3DALL_t030323_060003.bin'
    fn02 = path+'/3DALL_t030323_180002.bin'
    fn1 = [fn01, fn02]
    path = 'D:/simulation_results/low_density_cell/run_no_shrink_igrf_1_1/UA/data'
    fn01 = path+'/3DALL_t030323_060001.bin'
    fn02 = path+'/3DALL_t030323_180000.bin'
    fn2 = [fn01, fn02]

    nlat, slat = -30, -90
    alt = 300
    plt.figure()
    for k in range(2):
        g1 = gitm.read(fn1[k],varlist=['Rho'])
        g2 = gitm.read(fn2[k],varlist=['Rho'])
        ialt = np.argmin(np.abs(g1['Altitude'][0,0,:]/1000-alt))
        lt1 = 8 if k==0 else 4
        lt2 = 20 if k==0 else 16
        ut = g1['time'].hour + g1['time'].minute/60
        ilon1 = np.argmin(np.abs(g1['dLon'][2:-2,0,0]-(lt1-ut)*15%360))+2
        ilon2 = np.argmin(np.abs(g1['dLon'][2:-2,0,0]-(lt2-ut)*15%360))+2
        ilat = (g1['dLat'][0,:,0]>=slat) & (g1['dLat'][0,:,0]<=nlat)

        lat1 = g1['dLat'][ilon1,ilat,ialt]
        lat2 = g1['dLat'][ilon2,ilat,ialt]

        plt.subplot(2,2,k+1)

        g1rho1 = g1['Rho'][ilon1,ilat,ialt]
        g1rho2 = g1['Rho'][ilon2,ilat,ialt]
        plt.plot(-lat1, g1rho1, 'b')
        plt.plot(180+lat2, g1rho2, 'b')
        plt.plot([-lat1[0],180+lat2[0]],[g1rho1[0],g1rho2[0]],'b')

        g2rho1 = g2['Rho'][ilon1,ilat,ialt]
        g2rho2 = g2['Rho'][ilon2,ilat,ialt]
        plt.plot(-lat1, g2rho1, 'r')
        plt.plot(180+lat2, g2rho2, 'r')
        plt.plot([-lat1[0],180+lat2[0]],[g2rho1[0],g2rho2[0]],'r')
        plt.ylim(1e-11, 4e-11)

        plt.subplot(2,2,k+3)
        drho1 = (100*(g2['Rho']-g1['Rho'])/g1['Rho'])[ilon1,ilat,ialt]
        drho2 = (100*(g2['Rho']-g1['Rho'])/g1['Rho'])[ilon2,ilat,ialt]
        plt.plot(-lat1, drho1, 'k')
        plt.plot(180+lat2, drho2, 'k')
        plt.plot([-lat1[0],180+lat2[0]],[drho1[0],drho2[0]],'k')
        plt.ylim(-40, 40)
    return


if __name__ == '__main__':
    path = 'D:/simulation_results/low_density_cell/run_no_shrink_dipole_1_1/UA/data'
    fn1 = path+'/3DALL_t030322_060000.bin'
    fn2 = path+'/3DALL_t030322_180001.bin'
    path = 'D:/simulation_results/low_density_cell/run_no_shrink_igrf_1_1/UA/data'
    fn1 = path+'/3DALL_t030322_060001.bin'
    fn2 = path+'/3DALL_t030322_180002.bin'

    #  g = gitm.read(fn1)
    #  g['rhoi']=1.66e-27*(
    #      g['O!D2!U+!N']*32 + g['O_4SP_!U+!N']*16 + g['N!U+!N']*14 +
    #      g['N!D2!U+!N']*28 + g['NO!U+!N']*30 + g['He!U+!N']*4 + g['H!U+!N']*1 +
    #      g['O(!U2!NP)!U+!N']*16)
    #  g['vni'] = g['rhoi']/g['Rho']*g['Collision(in)']
    #  g['iondrag'] = g['vni']*np.sqrt(
    #      (g['V!Di!N (east)']-g['V!Dn!N (east)'])**2 +
    #      (g['V!Di!N (north)']-g['V!Dn!N (north)'])**2)

    #  #hc, hq = plot_polar_contour_vector(
    #  #    g,alt=300, zstr='e-', neuion='ion', useLT=True)
    #  #hc[0].set_clim(1e11,1.5e12)
    #  #hc[1].set_clim(1e11,1.5e12)

    #  #hc, hq = plot_polar_contour_vector(
    #  #    g,alt=300, zstr='Rho', neuion='neu', useLT=True)
    #  #hc[0].set_clim(1e-11,4e-11)
    #  #hc[1].set_clim(1e-11,4e-11)

    #  #hc, hq = plot_polar_contour_vector(
    #  #    g,alt=300, zstr='vni', vector=False, useLT=True)
    #  #hc[0].set_clim(5e-5,8e-4)
    #  #hc[1].set_clim(5e-5,8e-4)

    #  hc, hq = plot_polar_contour_vector(
    #      g,alt=300, zstr='iondrag', vector=False, useLT=True)
    #  hc[0].set_clim(0, 0.15)
    #  hc[1].set_clim(0, 0.15)

    #  g = gitm.read(fn2)
    #  g['rhoi']=1.66e-27*(
    #      g['O!D2!U+!N']*32 + g['O_4SP_!U+!N']*16 + g['N!U+!N']*14 +
    #      g['N!D2!U+!N']*28 + g['NO!U+!N']*30 + g['He!U+!N']*4 + g['H!U+!N']*1 +
    #      g['O(!U2!NP)!U+!N']*16)
    #  g['vni'] = g['rhoi']/g['Rho']*g['Collision(in)']
    #  g['iondrag'] = g['vni']*np.sqrt(
    #      (g['V!Di!N (east)']-g['V!Dn!N (east)'])**2 +
    #      (g['V!Di!N (north)']-g['V!Dn!N (north)'])**2)

    #  #hc, hq = plot_polar_contour_vector(
    #  #    g,alt=300, zstr='e-', neuion='ion', useLT=True)
    #  #hc[0].set_clim(1e11,1.5e12)
    #  #hc[1].set_clim(1e11,1.5e12)

    #  #hc, hq = plot_polar_contour_vector(
    #  #    g,alt=300, zstr='Rho', neuion='neu', useLT=True)
    #  #hc[0].set_clim(1e-11,4e-11)
    #  #hc[1].set_clim(1e-11,4e-11)

    #  #hc, hq = plot_polar_contour_vector(
    #  #    g,alt=300, zstr='vni', vector=False, useLT=True)
    #  #hc[0].set_clim(5e-5,8e-4)
    #  #hc[1].set_clim(5e-5,8e-4)

    #  hc, hq = plot_polar_contour_vector(
    #      g,alt=300, zstr='iondrag', vector=False, useLT=True)
    #  hc[0].set_clim(0, 0.15)
    #  hc[1].set_clim(0, 0.15)

    #  plt.show()

    #　plot_rho_wind_ut()
    #　plt.show()

    plot_06ut_18ut()
    plt.show()

    #plot_one_meridian()
    #plt.show()
