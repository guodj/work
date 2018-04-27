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
            useLT=useLT, lonticklabel=(1, 1, 1, 1))
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
            useLT=useLT, lonticklabel=(1, 1, 1, 1))
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

def plot_06ut_18ut(ns='SH',alt=150,tf=[1,0,0,0,0,0]):
    # 0:rho and wind; 1: diff rho and wind; 2:ion drag; 3:vni; 4:ion convection
    # 5:ion density; 6:electron density;
    path = '/home/guod/simulation_output/momentum_analysis/run_shrink_dipole_1_c1/UA/data'
    fn1 = path+'/3DALL_t030322_060002.bin'
    fn2 = path+'/3DALL_t030322_180002.bin'
    g11 = [gitm.read(fn1), gitm.read(fn2)] # shrink, dipole
    #print(g11[0]['Altitude'][0,0,np.argmin(np.abs(g11[0]['Altitude'][0,0,:]-150*1000))])

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

    if ns == 'SH':
        nlat, slat = -50, -90
    elif ns == 'NH':
        nlat, slat = 90, 50
    glats, glons = convert(-90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    glatn, glonn = convert(90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    fign=['(a)','(b)','(c)','(d)']
    # Rho and wind
    if tf[0]:
        plt.figure(1, figsize=(8, 7.55))
        for k in range(4):
            gtemp = g22 if k in [0, 1] else g21
            ax, projection, hc, hq = plot_polar_contour_vector(
                2, 2, k+1, gtemp[k%2], nlat, slat, alt, contour=True, zstr='Rho',
                vector=True, neuion='neu', useLT=True, coastlines=False,
                levels=np.linspace(0.8e-9, 1.6e-9, 21))
            if k in [0, 1]:
                ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
                ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
            ax.scatter(0, 90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            ax.scatter(0, -90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            if k in [0, 1]:
                plt.title('UT: '+gtemp[k%2]['time'].strftime('%H'), y=1.05)
            if k in [0, 2]:
                text = 'IGRF' if k==0 else 'Dipole'
                plt.text(
                    -0.2, 0.5, text, fontsize=14, verticalalignment='center',
                    transform=plt.gca().transAxes, rotation=90)
            plt.text(
                0, 1, fign[k], fontsize=14, transform=plt.gca().transAxes)
        plt.subplots_adjust(bottom=0.13,top=0.93, wspace=0.15, hspace=0.15)
        cax = plt.axes([0.3, 0.08, 0.4, 0.02])
        plt.colorbar(
            hc[0], cax=cax, orientation='horizontal',
            ticks=np.arange(0.8e-9, 1.7e-9, 0.2e-9))
        cax.set_xlabel(r'$\rho (kg/m^3)$')
        plt.quiverkey(hq[0], 0.5,0.12, 500, '500m/s',coordinates='figure')
        plt.savefig(savepath+'2'+ns+'RhoWind.jpeg')
        plt.savefig(savepath+'2'+ns+'RhoWind.eps')

    # ion drag
    if tf[1]:
        plt.figure(1, figsize=(8, 7.55))
        for k in range(4):
            gtemp = g22 if k in [0, 1] else g21
            calc_rhoi(gtemp[0])
            calc_rhoi(gtemp[1])
            gtemp[0]['iondrag'] = \
                gtemp[0]['rhoi']/gtemp[0]['Rho']*gtemp[0]['Collision(in)']*\
                np.sqrt((gtemp[0]['V!Di!N (east)']-gtemp[0]['V!Dn!N (east)'])**2 +\
                        (gtemp[0]['V!Di!N (north)']-gtemp[0]['V!Dn!N (north)'])**2)
            gtemp[1]['iondrag'] = \
                gtemp[1]['rhoi']/gtemp[1]['Rho']*gtemp[1]['Collision(in)']*\
                np.sqrt((gtemp[1]['V!Di!N (east)']-gtemp[1]['V!Dn!N (east)'])**2 +\
                        (gtemp[1]['V!Di!N (north)']-gtemp[1]['V!Dn!N (north)'])**2)
            ax, projection, hc, hq = plot_polar_contour_vector(
                2, 2, k+1, gtemp[k%2], nlat, slat, alt, contour=True, zstr='iondrag',
                vector=False, neuion='neu', useLT=True, coastlines=False,
                levels=np.linspace(0, 0.05, 21))
            if k in [0, 1]:
                ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
                ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
            ax.scatter(0, 90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            ax.scatter(0, -90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            if k in [0, 1]:
                plt.title('UT: '+gtemp[k%2]['time'].strftime('%H'), y=1.05)
            if k in [0, 2]:
                text = 'IGRF' if k==0 else 'Dipole'
                plt.text(
                    -0.2, 0.5, text, fontsize=14, verticalalignment='center',
                    transform=plt.gca().transAxes, rotation=90)
        plt.subplots_adjust(bottom=0.13,top=0.93, wspace=0.15, hspace=0.15)
        cax = plt.axes([0.3, 0.08, 0.4, 0.02])
        plt.colorbar(
            hc[0], cax=cax, orientation='horizontal',
            ticks=np.arange(0, 0.051, 0.01))
        cax.set_xlabel(r'Ion Drag $(m/s^2)$')
        plt.savefig(savepath+'3'+ns+'IonDrag.jpeg')
        plt.savefig(savepath+'3'+ns+'IonDrag.eps')

    # vni
    if tf[2]:
        plt.figure(1, figsize=(8, 7.55))
        for k in range(4):
            gtemp = g22 if k in [0, 1] else g21
            calc_rhoi(gtemp[0])
            calc_rhoi(gtemp[1])
            gtemp[0]['vni'] = \
                gtemp[0]['rhoi']/gtemp[0]['Rho']*gtemp[0]['Collision(in)']
            gtemp[1]['vni'] = \
                gtemp[1]['rhoi']/gtemp[1]['Rho']*gtemp[1]['Collision(in)']
            ax, projection, hc, hq = plot_polar_contour_vector(
                2, 2, k+1, gtemp[k%2], nlat, slat, alt, contour=True, zstr='vni',
                vector=False, neuion='neu', useLT=True, coastlines=False,
                levels=np.linspace(0, 0.0002, 21))
            if k in [0, 1]:
                ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
                ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
            ax.scatter(0, 90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            ax.scatter(0, -90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            if k in [0, 1]:
                plt.title('UT: '+gtemp[k%2]['time'].strftime('%H'), y=1.05)
            if k in [0, 2]:
                text = 'IGRF' if k==0 else 'Dipole'
                plt.text(
                    -0.2, 0.5, text, fontsize=14, verticalalignment='center',
                    transform=plt.gca().transAxes, rotation=90)
        plt.subplots_adjust(bottom=0.13,top=0.93, wspace=0.15, hspace=0.15)
        cax = plt.axes([0.3, 0.08, 0.4, 0.02])
        plt.colorbar(
            hc[0], cax=cax, orientation='horizontal',
            ticks=np.arange(0, 0.00021, 0.0001))
        cax.set_xlabel(r'$\nu_{ni}$ $(s^{-1})$')
        plt.savefig(savepath+'4'+ns+'vni.jpeg')
        plt.savefig(savepath+'4'+ns+'vni.eps')

    # ion convection
    if tf[3]:
        plt.figure(1, figsize=(8, 7.55))
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
            ax.scatter(0, 90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            ax.scatter(0, -90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            if k in [0, 1]:
                plt.title('UT: '+gtemp[k%2]['time'].strftime('%H'), y=1.05)
            if k in [0, 2]:
                text = 'IGRF' if k==0 else 'Dipole'
                plt.text(
                    -0.2, 0.5, text, fontsize=14, verticalalignment='center',
                    transform=plt.gca().transAxes, rotation=90)
        plt.subplots_adjust(bottom=0.13,top=0.93, wspace=0.15, hspace=0.15)
        cax = plt.axes([0.3, 0.08, 0.4, 0.02])
        plt.colorbar(
            hc[0], cax=cax, orientation='horizontal',
            ticks=np.arange(0, 1001, 200))
        cax.set_xlabel(r'Vi (m/s)')
        plt.savefig(savepath+'5'+ns+'IonV.jpeg')
        plt.savefig(savepath+'5'+ns+'IonV.eps')

    # ion density
    if tf[4]:
        plt.figure(1, figsize=(8, 7.55))
        for k in range(4):
            gtemp = g22 if k in [0, 1] else g21
            calc_rhoi(gtemp[0])
            calc_rhoi(gtemp[1])
            ax, projection, hc, hq = plot_polar_contour_vector(
                2, 2, k+1, gtemp[k%2], nlat, slat, alt, contour=True, zstr='rhoi',
                vector=False, neuion='neu', useLT=True, coastlines=False,
                levels=np.linspace(0, 12e-15, 21))
            if k in [0, 1]:
                ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
                ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
            ax.scatter(0, 90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            ax.scatter(0, -90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            if k in [0, 1]:
                plt.title('UT: '+gtemp[k%2]['time'].strftime('%H'), y=1.05)
            if k in [0, 2]:
                text = 'IGRF' if k==0 else 'Dipole'
                plt.text(
                    -0.2, 0.5, text, fontsize=14, verticalalignment='center',
                    transform=plt.gca().transAxes, rotation=90)
        plt.subplots_adjust(bottom=0.13,top=0.93, wspace=0.15, hspace=0.15)
        cax = plt.axes([0.3, 0.08, 0.4, 0.02])
        plt.colorbar(
            hc[0], cax=cax, orientation='horizontal',
            ticks=np.arange(0, 12.1e-15, 3e-15))
        cax.set_xlabel(r'$\rho_i$ $(kg/m^3)$')
        plt.savefig(savepath+'6'+ns+'Rhoi.jpeg')
        plt.savefig(savepath+'6'+ns+'Rhoi.eps')

    #  electron density
    if tf[5]:
        plt.figure(1, figsize=(8, 7.55))
        for k in range(4):
            gtemp = g22 if k in [0, 1] else g21
            calc_rhoi(gtemp[0])
            calc_rhoi(gtemp[1])
            ax, projection, hc, hq = plot_polar_contour_vector(
                2, 2, k+1, gtemp[k%2], nlat, slat, alt, contour=True, zstr='e-',
                vector=False, neuion='neu', useLT=True, coastlines=False,
                levels=np.linspace(0, 0.9e12, 21))
            if k in [0, 1]:
                ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=55, c='k')
                ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=55, c='k')
            ax.scatter(0, 90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            ax.scatter(0, -90, transform=ccrs.PlateCarree(), s=55, c='w',zorder=100)
            if k in [0, 1]:
                plt.title('UT: '+gtemp[k%2]['time'].strftime('%H'), y=1.05)
            if k in [0, 2]:
                text = 'IGRF' if k==0 else 'Dipole'
                plt.text(
                    -0.2, 0.5, text, fontsize=14, verticalalignment='center',
                    transform=plt.gca().transAxes, rotation=90)
        plt.subplots_adjust(bottom=0.13,top=0.93, wspace=0.15, hspace=0.15)
        cax = plt.axes([0.3, 0.08, 0.4, 0.02])
        plt.colorbar(
            hc[0], cax=cax, orientation='horizontal',
            ticks=np.arange(0, 0.31e12, 0.1e12))
        cax.set_xlabel(r'Ne ($m^{-3}$)')
        # plt.savefig(savepath+'7'+ns+'Ne.jpeg')
        # plt.savefig(savepath+'7'+ns+'Ne.eps')
    return


def plot_rho_wind_ut(ns='SH', alt=150):
    path = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_igrf_1_1/data'
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

    if ns=='SH':
        nlat, slat = -50, -90
    elif ns=='NH':
        nlat, slat = 90, 50
    glats, glons = convert(-90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    glatn, glonn = convert(90, 0, 0, date=dt.date(2003,1,1), a2g=True)
    fig1 = plt.figure(1, figsize=(5.28,8.49))
    titf = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)']
    for k in range(8):
        g2 = gitm.read(fn2[k])

        title = 'UT: '+g2['time'].strftime('%H')

        lon2, lat2, rho2 = g3ca.contour_data('Rho', g2, alt=alt)
        ilat2 = (lat2[:,0]>=slat) & (lat2[:,0]<=nlat)
        min_div_mean = \
            np.min(rho2[ilat2, :]) / \
            (np.mean(rho2[ilat2, :]*np.cos(lat2[ilat2, :]/180*np.pi)) /\
            np.mean(np.cos(lat2[ilat2,:]/180*np.pi)))
        print(ns,':  ',100*(min_div_mean-1))

        lon2, lat2, ewind2, nwind2 = g3ca.vector_data(g2,'neu',alt=alt)

        if k ==0:
            lonticklabel = [1, 0, 0, 1]
        elif k==1:
            lonticklabel = [1, 1, 0, 0]
        elif k in [2, 4]:
            lonticklabel = [0, 0, 0, 1]
        elif k in [3, 5]:
            lonticklabel = [0, 1, 0, 0]
        elif k==6:
            lonticklabel = [0, 0, 1, 1]
        elif k==7:
            lonticklabel = [0, 1, 1, 0]
        olat=False if k<7 else True

        # No shrink
        plt.figure(1)
        centrallon = g3ca.calculate_centrallon(g2, 'pol', useLT=True)
        ax, projection = gcc.create_map(
            4,2,k+1, 'pol', nlat, slat, centrallon, coastlines=False,
            dlat=10, useLT=True, lonticklabel=lonticklabel, olat=olat)
        hc1 = ax.contourf(lon2, lat2, rho2, np.linspace(0.8e-9, 1.6e-9),
            transform=ccrs.PlateCarree(),cmap='jet', extend='both')
        lon9, lat9, ewind9, nwind9 = g3ca.convert_vector(
            lon2, lat2, ewind2, nwind2, 'pol', projection)
        hq1 = ax.quiver(
            lon9,lat9,ewind9,nwind9,scale=1500,scale_units='inches',
            regrid_shape=20)
        ax.scatter(glons, glats, transform=ccrs.PlateCarree(), s=35, c='k')
        ax.scatter(glonn, glatn, transform=ccrs.PlateCarree(), s=35, c='k')
        ax.scatter(0, 90, transform=ccrs.PlateCarree(), s=35, c='w', zorder=100)
        ax.scatter(0, -90, transform=ccrs.PlateCarree(), s=35, c='w', zorder=100)
        plt.title(titf[k]+' '+title,x=0,y=0.93)

    plt.figure(1)
    plt.subplots_adjust(
        hspace=0.01, wspace=0.01, left=0.05, right=0.95, top=0.95, bottom=0.12)
    cax = plt.axes([0.25,0.07,0.5,0.02])
    hcb = plt.colorbar(
        hc1,cax=cax,orientation='horizontal',ticks=np.arange(0.8,1.7,0.2)*1e-9)
    hcb.set_label(r'$\rho$ (kg/$m^3$)')
    plt.quiverkey(hq1, 0.5,0.12, 500, '500m/s',coordinates='figure')
    plt.savefig(savepath+'1'+ns+'AllUT.eps')
    plt.savefig(savepath+'1'+ns+'AllUT.jpeg')

    return

if __name__ == '__main__':
    savepath ='/home/guod/Documents/Work/DensityCell/UTVariation/paper/'

    plt.close('all')

    plot_06ut_18ut(ns='SH', alt=150, tf=[1,1,1,1,1,1])
    #plot_rho_wind_ut(ns='NH', alt=150)

    #plot_one_meridian()
    #plt.show()
