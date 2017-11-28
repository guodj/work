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

def plot_vertical_wind(show=True, save=False):
    apex = Apex(date=2003)
    qlat, qlon = apex.convert(90, 0, source='apex', dest='geo', height=400)

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
        plt.savefig(spath+'vertical_wind_%s.pdf' % tstring)
    return

if __name__=='__main__':
    from PyPDF2 import PdfFileMerger, PdfFileReader
    whichlt=7
    trange = pd.date_range(
            '2003-03-22 00:00:00', '2003-03-22 04:00:00', freq='10min')
    spath = '/home/guod/Documents/work/fig/density_cell/'\
            'why_no_low_density_cell_at_high_latitude/snapshot_alt_lat/'
    outpdf_vert_wind = PdfFileMerger()
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
        plot_vertical_wind(show=False, save=True)
        outpdf_vert_wind.append(PdfFileReader(open(spath+'vertical_wind_%s.pdf' % tstring, 'rb')))
    outpdf_vert_wind.write(spath+'vertical_wind.pdf')

