#Global imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
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

def plot_vertical_wind(show=True, save=False):
    plt.close('all')
    plt.figure()
    ilt = np.argmin(np.abs(g1a['LT'][2:-2, 0, 0]-whichlt))+2
    vert_wind1 = g1a['V!Dn!N (up)'][ilt, 2:-2, 2:-2]
    vert_wind2 = g2a['V!Dn!N (up)'][ilt, 2:-2, 2:-2]
    vert_wind = np.array(vert_wind2 - vert_wind1).T
    lat = np.array(g1a['dLat'][ilt, 2:-2, 2:-2]).T
    alt = np.array(g1a['Altitude'][ilt, 2:-2, 2:-2]/1000).T
    plt.contourf(lat, alt, vert_wind, levels=np.linspace(-10, 10, 21),
                 cmap='seismic', extend='both')
    plt.xlim(-90, -40)
    plt.xlabel('Latitude')
    plt.ylabel('Altitude')
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'02_vertical_wind_diff_%s.pdf' % tstring)
    return

def plot_rho_diff(show=True, save=False):
    plt.close('all')
    plt.figure()
    ilt = np.argmin(np.abs(g1a['LT'][2:-2, 0, 0]-whichlt))+2
    rho1 = g1a['Rho'][ilt,2:-2,2:-2]
    rho2 = g2a['Rho'][ilt,2:-2,2:-2]
    drho = np.array(100*(rho2-rho1)/rho1).T
    lat = np.array(g1a['dLat'][ilt, 2:-2, 2:-2]).T
    alt = np.array(g1a['Altitude'][ilt, 2:-2, 2:-2]/1000).T
    plt.contourf(lat, alt, drho, levels=np.linspace(-30, 30, 21),
                 cmap='seismic', extend='both')
    plt.xlim(-90, -40)
    plt.xlabel('Latitude')
    plt.ylabel('Altitude')
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'01_rho_diff_%s.pdf' % tstring)
    return

def plot_divrhov_rho_diff(show=True, save=False):
    plt.close('all')
    plt.figure()
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
    div_rhov1 = (cr.calc_div_hozt(lon1, lat1, alt1, rho1*nwind1, rho1*ewind1)\
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
    div_rhov2 = (cr.calc_div_hozt(lon2, lat2, alt2, rho2*nwind2, rho2*ewind2)\
               +cr.calc_div_vert(alt2, rho2*uwind2))/rho2

    ilt = np.argmin(np.abs(g1a['LT'][2:-2, 0, 0]-whichlt))+2
    divrhov1 = div_rhov1[ilt,2:-2,2:-2]
    divrhov2 = div_rhov2[ilt,2:-2,2:-2]
    ddivrhov = np.array(divrhov1-divrhov2).T
    lat = np.array(g1a['dLat'][ilt, 2:-2, 2:-2]).T
    alt = np.array(g1a['Altitude'][ilt, 2:-2, 2:-2]/1000).T
    plt.contourf(lat, alt, ddivrhov, levels=np.linspace(-1, 1, 21)*1e-4,
                 cmap='seismic', extend='both')
    plt.xlim(-90, -40)
    plt.xlabel('Latitude')
    plt.ylabel('Altitude')
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'03_divrhov_rho_diff_%s.pdf' % tstring)
    return

def plot_vert_divv_diff(show=True, save=False):
    plt.close('all')
    plt.figure()

    velr = np.array(g1a['V!Dn!N (up)'])
    divv1 = cr.calc_div_vert(g1a['Altitude'], velr)

    velr = np.array(g2a['V!Dn!N (up)'])
    divv2 = cr.calc_div_vert(g2a['Altitude'], velr)

    ilt = np.argmin(np.abs(g1a['LT'][2:-2, 0, 0]-whichlt))+2
    divv1 = divv1[ilt,2:-2,2:-2]
    divv2 = divv2[ilt,2:-2,2:-2]
    ddivv = np.array(divv1-divv2).T
    lat = np.array(g1a['dLat'][ilt, 2:-2, 2:-2]).T
    alt = np.array(g1a['Altitude'][ilt, 2:-2, 2:-2]/1000).T
    plt.contourf(lat, alt, ddivv, levels=np.linspace(-1, 1, 21)*1e-4,
                 cmap='seismic', extend='both')
    plt.xlim(-90, -40)
    plt.xlabel('Latitude')
    plt.ylabel('Altitude')
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'04_vert_divv_diff%s.pdf' % tstring)
    return

def plot_hozt_divv_diff(show=True, save=False):
    plt.close('all')
    plt.figure()

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

    ilt = np.argmin(np.abs(g1a['LT'][2:-2, 0, 0]-whichlt))+2
    divv1 = divv1[ilt,2:-2,2:-2]
    divv2 = divv2[ilt,2:-2,2:-2]
    ddivv = np.array(divv1-divv2).T
    lat = np.array(g1a['dLat'][ilt, 2:-2, 2:-2]).T
    alt = np.array(g1a['Altitude'][ilt, 2:-2, 2:-2]/1000).T
    plt.contourf(lat, alt, ddivv, levels=np.linspace(-1, 1, 21)*1e-4,
                 cmap='seismic', extend='both')
    plt.xlim(-90, -40)
    plt.xlabel('Latitude')
    plt.ylabel('Altitude')
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'05_hozt_divv_diff%s.pdf' % tstring)
    return

def plot_vert_vgradrho_rho_diff(show=True, save=False):
    plt.close('all')
    plt.figure()

    rho1 = np.array(g1a['Rho'])
    vgradrho1 = \
            g1a['V!Dn!N (up)']\
            *cr.calc_rusanov_alts_ausm(g1a['Altitude'],rho1)/g1a['Rho']

    rho2 = np.array(g2a['Rho'])
    vgradrho2 = \
            g2a['V!Dn!N (up)']\
            *cr.calc_rusanov_alts_ausm(g2a['Altitude'],rho2)/g2a['Rho']

    ilt = np.argmin(np.abs(g1a['LT'][2:-2, 0, 0]-whichlt))+2
    vgradrho1 = vgradrho1[ilt,2:-2,2:-2]
    vgradrho2 = vgradrho2[ilt,2:-2,2:-2]
    dvgradrho = np.array(vgradrho1-vgradrho2).T
    lat = np.array(g1a['dLat'][ilt, 2:-2, 2:-2]).T
    alt = np.array(g1a['Altitude'][ilt, 2:-2, 2:-2]/1000).T
    plt.contourf(lat, alt, dvgradrho, levels=np.linspace(-1, 1, 21)*1e-4,
                 cmap='seismic', extend='both')
    plt.xlim(-90, -40)
    plt.xlabel('Latitude')
    plt.ylabel('Altitude')
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'06_vert_vgradrho_rho_diff_%s.pdf' % tstring)
    return

def plot_hozt_vgradrho_rho_diff(show=True, save=False):
    plt.close('all')
    plt.figure()

    lon1 = np.array(g1a['Longitude'])
    lat1 = np.array(g1a['Latitude'])
    alt1 = np.array(g1a['Altitude'])
    Re = 6371*1000 # Earth radius, unit: m
    RR = Re+alt1
    omega = 2*np.pi/(24*3600)
    rho1 = np.array(g1a['Rho'])
    nwind1 = np.array(g1a['V!Dn!N (north)'])
    ewind1 = np.array(g1a['V!Dn!N (east)']) + omega*RR*np.cos(lat1)
    vgradrho1 = (nwind1*cr.calc_rusanov_lats(lat1,alt1,rho1)\
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
    vgradrho2 = (nwind2*cr.calc_rusanov_lats(lat2,alt2,rho2)\
               +ewind2*cr.calc_rusanov_lons(lon2,lat2,alt2,rho2))/rho2

    ilt = np.argmin(np.abs(g1a['LT'][2:-2, 0, 0]-whichlt))+2
    vgradrho1 = vgradrho1[ilt,2:-2,2:-2]
    vgradrho2 = vgradrho2[ilt,2:-2,2:-2]
    dvgradrho = np.array(vgradrho1-vgradrho2).T
    lat = np.array(g1a['dLat'][ilt, 2:-2, 2:-2]).T
    alt = np.array(g1a['Altitude'][ilt, 2:-2, 2:-2]/1000).T
    plt.contourf(lat, alt, dvgradrho, levels=np.linspace(-1, 1, 21)*1e-4,
                 cmap='seismic', extend='both')
    plt.xlim(-90, -40)
    plt.xlabel('Latitude')
    plt.ylabel('Altitude')
    plt.text(0.5, 0.95, 'Time: '+tstring, fontsize=15,
            horizontalalignment='center', transform=plt.gcf().transFigure)
    if show:
        plt.show()
    if save:
        plt.savefig(spath+'07_hozt_vgradrho_rho_diff_%s.pdf' % tstring)
    return

if __name__=='__main__':
    from PyPDF2 import PdfFileMerger, PdfFileReader
    whichlt=7
    trange = pd.date_range(
            '2003-03-22 00:00:00', '2003-03-22 04:00:00', freq='10min')
    spath = '/home/guod/Documents/work/fig/density_cell/'\
            'why_no_low_density_cell_at_high_latitude/lat_alt/'
    outpdf_vert_wind = PdfFileMerger()
    outpdf_rho_diff = PdfFileMerger()
    outpdf_divrhov_rho_diff = PdfFileMerger()
    outpdf_vert_divv_diff = PdfFileMerger()
    outpdf_hozt_divv_diff = PdfFileMerger()
    outpdf_vert_vgradrho_diff = PdfFileMerger()
    outpdf_hozt_vgradrho_diff = PdfFileMerger()
    for t in trange:
        tstring = t.strftime('%H%M')
        filename = glob.glob(
                '/media/guod/wd2t/simulation_output/momentum_analysis/run_shrink_iondrift_4_c1'\
                '/data/3DALL_t030322_%s*.bin' % tstring)[0]
        g1a = gitm.GitmBin(filename)
        filename = glob.glob(
                '/media/guod/wd2t/simulation_output/momentum_analysis/run_no_shrink_iondrift_4_1'\
                '/data/3DALL_t030322_%s*.bin' % tstring)[0]
        g2a = gitm.GitmBin(filename)
        # plot_vertical_wind(show=False, save=True)
        # outpdf_vert_wind.append(PdfFileReader(open(spath+'02_vertical_wind_diff_%s.pdf' % tstring, 'rb')))
        # plot_rho_diff(show=False, save=True)
        # outpdf_rho_diff.append(PdfFileReader(open(spath+'01_rho_diff_%s.pdf' % tstring, 'rb')))
        # plot_divrhov_rho_diff(show=False, save=True)
        # outpdf_divrhov_rho_diff.append(PdfFileReader(open(spath+'03_divrhov_rho_diff_%s.pdf' % tstring, 'rb')))
        # plot_vert_divv_diff(show=False, save=True)
        # outpdf_vert_divv_diff.append(PdfFileReader(open(spath+'04_vert_divv_diff%s.pdf' % tstring, 'rb')))
        # plot_hozt_divv_diff(show=False, save=True)
        # outpdf_hozt_divv_diff.append(PdfFileReader(open(spath+'05_hozt_divv_diff%s.pdf' % tstring, 'rb')))
        plot_vert_vgradrho_rho_diff(show=False, save=True)
        outpdf_vert_vgradrho_diff.append(PdfFileReader(open(spath+'06_vert_vgradrho_rho_diff_%s.pdf' % tstring, 'rb')))
        # plot_hozt_vgradrho_rho_diff(show=False, save=True)
        # outpdf_hozt_vgradrho_diff.append(PdfFileReader(open(spath+'07_hozt_vgradrho_rho_diff_%s.pdf' % tstring, 'rb')))
    #outpdf_vert_wind.write(spath+'02_vertical_wind.pdf')
    #outpdf_rho_diff.write(spath+'01_rho_diff.pdf')
    #outpdf_rho_diff.write(spath+'03_divrhov_rho_diff.pdf')
    #outpdf_vert_divv_diff.write(spath+'04_vert_divv_diff.pdf')
    #outpdf_hozt_divv_diff.write(spath+'05_hozt_divv_diff.pdf')
    outpdf_vert_vgradrho_diff.write(spath+'06_vert_vgradrho_rho_diff.pdf')
    #outpdf_hozt_vgradrho_diff.write(spath+'07_hozt_vgradrho_rho_diff.pdf')

