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
import gitm_divergence as gd

def rho_diff(readdata=False,show=True):
    zz = []
    plt.close('all')
    if readdata:
        for ifn, fn in enumerate(fns1):
            # read gitm data
            g1a = gitm.GitmBin(fns1[ifn])
            g2a = gitm.GitmBin(fns2[ifn])
            alt_ind = np.argmin(np.abs(g1a['Altitude'][0,0,2:-2]/1000-alt))+2
            alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
            lt_ind = np.argmin(np.abs(g1a['LT'][2:-2, 0,0]-LT))+2
            rho_diff = ((100*(g2a['Rho']-g1a['Rho'])/g1a['Rho'])[lt_ind,2:-2,alt_ind])\
                       .reshape(-1, 1)
            zz.append(rho_diff)
        zz=np.concatenate(zz,axis=1)
        xtime = (timeidx-stime)/pd.Timedelta('1h')
        ylat = np.array(g1a['Latitude'][0,2:-2,0]/np.pi*180)
        pd.to_pickle((xtime,ylat,zz), '/home/guod/tmp/rho_diff.dat')
    xtime,ylat,zz = pd.read_pickle('/home/guod/tmp/rho_diff.dat')
    xtime,ylat = np.meshgrid(xtime,ylat)
    hct=plt.contourf(xtime,ylat,zz,np.linspace(-30,30,21),cmap='seismic',
                     extend='both')
    hcb=plt.colorbar(hct)
    hcb.set_label(r'$100*\frac{\rho_2-\rho_1}{\rho_1}$ (%)')
    plt.ylim(-90, -40)
    plt.ylabel('Latitude')
    plt.xlabel('Time (hour)')
    if show:
        plt.show()
    plt.savefig(spath+'rho_diff.pdf')
    return

def pressure_diff(readdata=False,show=True):
    zz = []
    plt.close('all')
    if readdata:
        for ifn, fn in enumerate(fns1):
            # read gitm data
            g1a = gitm.GitmBin(fns1[ifn])
            g2a = gitm.GitmBin(fns2[ifn])
            gp.calc_pressure(g1a)
            gp.calc_pressure(g2a)
            alt_ind = np.argmin(np.abs(g1a['Altitude'][0,0,2:-2]/1000-alt))+2
            alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
            lt_ind = np.argmin(np.abs(g1a['LT'][2:-2, 0,0]-LT))+2
            pressure_diff = ((100*(g2a['pressure']-g1a['pressure'])/g1a['pressure'])\
                             [lt_ind,2:-2,alt_ind]).reshape(-1, 1)
            zz.append(pressure_diff)
        zz=np.concatenate(zz,axis=1)
        xtime = (timeidx-stime)/pd.Timedelta('1h')
        ylat = np.array(g1a['Latitude'][0,2:-2,0]/np.pi*180)
        pd.to_pickle((xtime,ylat,zz), '/home/guod/tmp/pressure_diff.dat')
    xtime,ylat,zz = pd.read_pickle('/home/guod/tmp/pressure_diff.dat')
    xtime,ylat = np.meshgrid(xtime,ylat)
    hct=plt.contourf(xtime,ylat,zz,np.linspace(-30,30,21),cmap='seismic',
                     extend='both')
    hcb=plt.colorbar(hct)
    hcb.set_label(r'$100*\frac{P_2-P_1}{P_1}$ (%)')
    plt.ylim(-90, -40)
    plt.ylabel('Latitude')
    plt.xlabel('Time (hour)')
    if show:
        plt.show()
    plt.savefig(spath+'pressure_diff.pdf')
    return

def wind_diff(readdata=False,show=True):
    zz = []
    plt.close('all')
    if readdata:
        for ifn, fn in enumerate(fns1):
            # read gitm data
            g1a = gitm.GitmBin(fns1[ifn])
            g2a = gitm.GitmBin(fns2[ifn])
            alt_ind = np.argmin(np.abs(g1a['Altitude'][0,0,2:-2]/1000-alt))+2
            alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
            lt_ind = np.argmin(np.abs(g1a['LT'][2:-2, 0,0]-LT))+2
            nwind_diff = g2a['']
        pd.to_pickle((xtime,ylat,zz), '/home/guod/tmp/rho_diff.dat')
    xtime,ylat,zz = pd.read_pickle('/home/guod/tmp/rho_diff.dat')
    xtime,ylat = np.meshgrid(xtime,ylat)
    hct=plt.contourf(xtime,ylat,zz,np.linspace(-30,30,21),cmap='seismic',
                     extend='both')
    hcb=plt.colorbar(hct)
    hcb.set_label(r'$100*\frac{\rho_2-\rho_1}{\rho_1}$ (%)')
    plt.ylim(-90, -40)
    plt.ylabel('Latitude')
    plt.xlabel('Time (hour)')
    if show:
        plt.show()
    plt.savefig(spath+'rho_diff.pdf')
    return

def vert_vgradrho_diff(readdata=False,show=True):
    zz = []
    plt.close('all')
    if readdata:
        for ifn, fn in enumerate(fns1):
            # read gitm data
            g1a = gitm.GitmBin(fns1[ifn])
            g2a = gitm.GitmBin(fns2[ifn])
            alt_ind = np.argmin(np.abs(g1a['Altitude'][0,0,2:-2]/1000-alt))+2
            alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
            lt_ind = np.argmin(np.abs(g1a['LT'][2:-2, 0,0]-LT))+2

            # vgradrho1
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

            # vgradrho2
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

            vgradrho_diff = ((vgradrho1-vgradrho2)[lt_ind,2:-2,alt_ind]).reshape(-1, 1)

            zz.append(vgradrho_diff)
        zz=np.concatenate(zz,axis=1)
        xtime = (timeidx-stime)/pd.Timedelta('1h')
        ylat = np.array(g1a['Latitude'][0,2:-2,0]/np.pi*180)
        pd.to_pickle((xtime,ylat,zz), '/home/guod/tmp/vert_vgradrho_diff.dat')
    xtime,ylat,zz = pd.read_pickle('/home/guod/tmp/vert_vgradrho_diff.dat')
    xtime,ylat = np.meshgrid(xtime,ylat)
    hct=plt.contourf(xtime,ylat,zz,np.linspace(-1.2,1.2,21)*1e-15,cmap='seismic',
                     extend='both')
    hcb=plt.colorbar(hct)
    hcb.set_label(r'$\vec{u}\cdot\nabla\rho$ (up)')
    plt.ylim(-90, -40)
    plt.ylabel('Latitude')
    plt.xlabel('Time (hour)')
    if show:
        plt.show()
    plt.savefig(spath+'vert_vgradrho.pdf')
    return

def hozt_vgradrho_diff(readdata=False,show=True):
    zz = []
    plt.close('all')
    if readdata:
        for ifn, fn in enumerate(fns1):
            # read gitm data
            g1a = gitm.GitmBin(fns1[ifn])
            g2a = gitm.GitmBin(fns2[ifn])
            alt_ind = np.argmin(np.abs(g1a['Altitude'][0,0,2:-2]/1000-alt))+2
            alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
            lt_ind = np.argmin(np.abs(g1a['LT'][2:-2, 0,0]-LT))+2

            # vgradrho1
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

            # vgradrho2
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

            vgradrho_diff = ((vgradrho1-vgradrho2)[lt_ind,2:-2,alt_ind]).reshape(-1, 1)

            zz.append(vgradrho_diff)
        zz=np.concatenate(zz,axis=1)
        xtime = (timeidx-stime)/pd.Timedelta('1h')
        ylat = np.array(g1a['Latitude'][0,2:-2,0]/np.pi*180)
        pd.to_pickle((xtime,ylat,zz), '/home/guod/tmp/hozt_vgradrho_diff.dat')
    xtime,ylat,zz = pd.read_pickle('/home/guod/tmp/hozt_vgradrho_diff.dat')
    xtime,ylat = np.meshgrid(xtime,ylat)
    hct=plt.contourf(xtime,ylat,zz,np.linspace(-1.2,1.2,21)*1e-15,cmap='seismic',
                     extend='both')
    hcb=plt.colorbar(hct)
    hcb.set_label(r'$\vec{u}\cdot\nabla\rho$ (horizontal)')
    plt.ylim(-90, -40)
    plt.ylabel('Latitude')
    plt.xlabel('Time (hour)')
    if show:
        plt.show()
    plt.savefig(spath+'hozt_vgradrho.pdf')
    return

def vert_rhodivv_diff(readdata=False,show=True):
    zz = []
    plt.close('all')
    if readdata:
        for ifn, fn in enumerate(fns1):
            # read gitm data
            g1a = gitm.GitmBin(fns1[ifn])
            g2a = gitm.GitmBin(fns2[ifn])
            alt_ind = np.argmin(np.abs(g1a['Altitude'][0,0,2:-2]/1000-alt))+2
            alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
            lt_ind = np.argmin(np.abs(g1a['LT'][2:-2, 0,0]-LT))+2

            # vgradrho1
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

            # vgradrho2
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

            rhodivv_diff = ((rhodivv1-rhodivv2)[lt_ind,2:-2,alt_ind]).reshape(-1, 1)

            zz.append(rhodivv_diff)
        zz=np.concatenate(zz,axis=1)
        xtime = (timeidx-stime)/pd.Timedelta('1h')
        ylat = np.array(g1a['Latitude'][0,2:-2,0]/np.pi*180)
        pd.to_pickle((xtime,ylat,zz), '/home/guod/tmp/vert_rhodivv_diff.dat')
    xtime,ylat,zz = pd.read_pickle('/home/guod/tmp/vert_rhodivv_diff.dat')
    xtime,ylat = np.meshgrid(xtime,ylat)
    hct=plt.contourf(xtime,ylat,zz,np.linspace(-1.2,1.2,21)*1e-15,cmap='seismic',
                     extend='both')
    hcb=plt.colorbar(hct)
    hcb.set_label(r'$\rho\nabla\cdot\vec{u}$ (up)')
    plt.ylim(-90, -40)
    plt.ylabel('Latitude')
    plt.xlabel('Time (hour)')
    if show:
        plt.show()
    plt.savefig(spath+'vert_rhodivv.pdf')
    return

def hozt_rhodivv_diff(readdata=False,show=True):
    zz = []
    plt.close('all')
    if readdata:
        for ifn, fn in enumerate(fns1):
            # read gitm data
            g1a = gitm.GitmBin(fns1[ifn])
            g2a = gitm.GitmBin(fns2[ifn])
            alt_ind = np.argmin(np.abs(g1a['Altitude'][0,0,2:-2]/1000-alt))+2
            alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
            lt_ind = np.argmin(np.abs(g1a['LT'][2:-2, 0,0]-LT))+2

            # vgradrho1
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

            # vgradrho2
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

            rhodivv_diff = ((rhodivv1-rhodivv2)[lt_ind,2:-2,alt_ind]).reshape(-1, 1)

            zz.append(rhodivv_diff)
        zz=np.concatenate(zz,axis=1)
        xtime = (timeidx-stime)/pd.Timedelta('1h')
        ylat = np.array(g1a['Latitude'][0,2:-2,0]/np.pi*180)
        pd.to_pickle((xtime,ylat,zz), '/home/guod/tmp/hozt_rhodivv_diff.dat')
    xtime,ylat,zz = pd.read_pickle('/home/guod/tmp/hozt_rhodivv_diff.dat')
    xtime,ylat = np.meshgrid(xtime,ylat)
    hct=plt.contourf(xtime,ylat,zz,np.linspace(-1.2,1.2,21)*1e-15,cmap='seismic',
                     extend='both')
    hcb=plt.colorbar(hct)
    hcb.set_label(r'$\rho\nabla\cdot\vec{u}$ (horizontal)')
    plt.ylim(-90, -40)
    plt.ylabel('Latitude')
    plt.xlabel('Time (hour)')
    if show:
        plt.show()
    plt.savefig(spath+'hozt_rhodivv.pdf')
    return

def divrhov_diff(readdata=False,show=True):
    zz = []
    plt.close('all')
    if readdata:
        for ifn, fn in enumerate(fns1):
            # read gitm data
            g1a = gitm.GitmBin(fns1[ifn])
            g2a = gitm.GitmBin(fns2[ifn])
            alt_ind = np.argmin(np.abs(g1a['Altitude'][0,0,2:-2]/1000-alt))
            alt_str = '%6.2f' % (g1a['Altitude'][0, 0, alt_ind]/1000)
            lt_ind = np.argmin(np.abs(g1a['LT'][2:-2, 0,0]-LT))+2

            # vgradrho1
            lon1 = np.array(g1a['Longitude'])[2:-2, 2:-2, 2:-2]
            lat1 = np.array(g1a['Latitude'])[2:-2, 2:-2, 2:-2]
            alt1 = np.array(g1a['Altitude'])[2:-2, 2:-2, 2:-2]
            Re = 6371*1000 # Earth radius, unit: m
            RR = Re+alt1
            omega = 2*np.pi/(24*3600)
            rho1 = np.array(g1a['Rho'])[2:-2, 2:-2, 2:-2]
            nwind1 = np.array(g1a['V!Dn!N (north)'])[2:-2, 2:-2, 2:-2]
            ewind1 = np.array(g1a['V!Dn!N (east)'])[2:-2, 2:-2, 2:-2] \
                     + omega*RR*np.cos(lat1)
            uwind1 = np.array(g1a['V!Dn!N (up)'])[2:-2, 2:-2, 2:-2]
            div_rhov1 = (
                    1.0/(RR**2)
                  * np.gradient((RR**2)*rho1*uwind1, axis=2) / np.gradient(alt1, axis=2)
                  + 1.0/(RR*np.cos(lat1))
                  * (np.gradient(rho1*nwind1*np.cos(lat1), axis=1) / np.gradient(lat1, axis=1)
                     + np.gradient(rho1*ewind1, axis=0) / np.gradient(lon1, axis=0)))

            # vgradrho2
            lon2 = np.array(g2a['Longitude'])[2:-2, 2:-2, 2:-2]
            lat2 = np.array(g2a['Latitude'])[2:-2, 2:-2, 2:-2]
            alt2 = np.array(g2a['Altitude'])[2:-2, 2:-2, 2:-2]
            Re = 6371*1000 # Earth radius, unit: m
            RR = Re+alt2
            omega = 2*np.pi/(24*3600)
            rho2 = np.array(g2a['Rho'])[2:-2, 2:-2, 2:-2]
            nwind2 = np.array(g2a['V!Dn!N (north)'])[2:-2, 2:-2, 2:-2]
            ewind2 = np.array(g2a['V!Dn!N (east)'])[2:-2, 2:-2, 2:-2]\
                     + omega*RR*np.cos(lat2)
            uwind2 = np.array(g2a['V!Dn!N (up)'])[2:-2, 2:-2, 2:-2]
            div_rhov2 = (
                    1.0/(RR**2)
                  * np.gradient((RR**2)*rho2*uwind2, axis=2) / np.gradient(alt2, axis=2)
                  + 1.0/(RR*np.cos(lat2))
                  * (np.gradient(rho2*nwind2*np.cos(lat2), axis=1) / np.gradient(lat2, axis=1)
                     + np.gradient(rho2*ewind2, axis=0) / np.gradient(lon2, axis=0)))

            div_rhov_diff = ((div_rhov1-div_rhov2)[lt_ind,2:-2,alt_ind]).reshape(-1, 1)

            zz.append(div_rhov_diff)
        zz=np.concatenate(zz,axis=1)
        xtime = (timeidx-stime)/pd.Timedelta('1h')
        ylat = np.array(g1a['Latitude'][0,2:-2,0]/np.pi*180)
        pd.to_pickle((xtime,ylat,zz), '/home/guod/tmp/divrhov_diff.dat')
    xtime,ylat,zz = pd.read_pickle('/home/guod/tmp/divrhov_diff.dat')
    xtime,ylat = np.meshgrid(xtime,ylat)
    hct=plt.contourf(xtime,ylat,zz,np.linspace(-1.2,1.2,21)*1e-15,cmap='seismic',
                     extend='both')
    hcb=plt.colorbar(hct)
    hcb.set_label(r'$\nabla\cdot(\rho\vec{u})$')
    plt.ylim(-90, -40)
    plt.ylabel('Latitude')
    plt.xlabel('Time (hour)')
    if show:
        plt.show()
    plt.savefig(spath+'divrhov.pdf')
    return


if __name__=='__main__':
    import gc
    stime = pd.Timestamp('2003-03-22 00:00:00')
    etime = pd.Timestamp('2003-03-22 06:00:00')
    timeidx = pd.DatetimeIndex(start=stime, end=etime, freq='5min')
    fp1 = '/home/guod/simulation_output/momentum_analysis/'\
         + 'run_shrink_iondrift_3_continue/data/'
    fns1 = [glob.glob(fp1+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
           for k in timeidx]
    fns11 = [glob.glob(fp1+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
            for k in timeidx]
    fp2 = '/home/guod/simulation_output/momentum_analysis/'\
         + 'run_no_shrink_iondrift_3/data/'
    fns2 = [glob.glob(fp2+'3DALL_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
            for k in timeidx]
    fns22 = [glob.glob(fp2+'3DMOM_t'+k.strftime('%y%m%d_%H%M')+'*.bin')[0]
            for k in timeidx]

    # save path
    spath = '/home/guod/Documents/work/fig/density_cell/' \
           'why_no_low_density_cell_at_high_latitude/lat_time/'

    alt = 400
    LT = 7
    divrhov_diff(readdata=True)
    gc.collect()
