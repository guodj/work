"""
The thermospheric neutral density at 150 km is lower at the pole than the
surrounding area even the ion drift is set to zero. Why is this happening?
"""
import gitm_new as gitm
import gitm_3D_const_alt as g3ca
import gitm_3D_const_lon as g3cl
import gitm_create_coordinate as gcc
import matplotlib.pyplot as plt
import numpy as np
import os

gallfn = '/Users/guod/data/GITMOutput/run_shrink_notides_tmp/UA/data/3DALL_t030323_000004.bin'
gthmfn = '/Users/guod/data/GITMOutput/run_shrink_notides_tmp/UA/data/3DTHM_t030323_000004.bin'
gall = gitm.read(gallfn)
gthm = gitm.read(gthmfn)
savepath='/Users/guod/tmp/'

def multiple_alts_const_alt():
    print('Density and temperature distribution at different heights.')
    for alt in [100, 105, 110, 115, 120, 125, 135, 145]:
        if not os.path.exists(savepath+str(alt)):
            os.mkdir(savepath+str(alt))
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='Rho',levels=None, vector=True,
            neuion='neu', scale=500, useLT=True)
        plt.colorbar(hc)
        plt.savefig(savepath+str(alt)+'/001_Rho.png')
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='Temperature',levels=None, vector=True,
            neuion='neu', scale=500, useLT=True)
        plt.savefig(savepath+str(alt)+'/002_Temperature.png')
    print('Distributions of energy terms at different heights.')
    for alt in [100, 105, 110, 115, 120, 125, 135, 145]:
        if not os.path.exists(savepath+str(alt)):
            os.mkdir(savepath+str(alt))
        fig, hc = g3ca.test(
            gthm, alt=alt, contour=True, zstr='EUV Heating',
            levels=np.linspace(-0.01, 0.01, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/003_EUVHeating.png')
        fig, hc = g3ca.test(
            gthm, alt=alt, contour=True, zstr='Conduction',
            levels=np.linspace(-0.01, 0.01, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/003_Conduction.png')
        fig, hc = g3ca.test(
            gthm, alt=alt, contour=True, zstr='Chemical Heating',
            levels=np.linspace(-0.01, 0.01, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/003_ChemicalHeating.png')
        fig, hc = g3ca.test(
            gthm, alt=alt, contour=True, zstr='Auroral Heating',
            levels=np.linspace(-0.01, 0.01, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/003_AuroralHeating.png')
        fig, hc = g3ca.test(
            gthm, alt=alt, contour=True, zstr='Joule Heating',
            levels=np.linspace(-0.01, 0.01, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/003_JouleHeating.png')
        fig, hc = g3ca.test(
            gthm, alt=alt, contour=True, zstr='NO Cooling',
            levels=np.linspace(-0.01, 0.01, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/003_NOCooling.png')
        fig, hc = g3ca.test(
            gthm, alt=alt, contour=True, zstr='O Cooling',
            levels=np.linspace(-0.01, 0.01, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/003_OCooling.png')
        gthm['total energy'] =\
            gthm['EUV Heating'] + gthm['Conduction'] + \
            gthm['Chemical Heating'] + gthm['Auroral Heating'] + \
            gthm['Joule Heating'] - gthm['NO Cooling'] - gthm['O Cooling']
        fig, hc = g3ca.test(
            gthm, alt=alt, contour=True, zstr='total energy',
            levels=np.linspace(-0.01, 0.01, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/003_TotalEnergy.png')
        pass
    return

def height_profile(g, var='Rho', whichlt=12, altmin=100, altmax=150, log=True,
                   vector=True, levels=None, cmap='jet'):
    ilt = np.argmin(np.abs(g['LT'][2:-2, 0, 0]-whichlt))+2
    ialt = (g['Altitude'][0,0,:]>=altmin*1000) & \
           (g['Altitude'][0,0,:]<=altmax*1000)
    lat = g['dLat'][ilt, :, ialt]
    alt = g['Altitude'][ilt,:,ialt]/1000
    variable = g[var][ilt,:,ialt]
    if log:
        variable = np.log10(variable)
    if vector:
        ue = g['V!Dn!N (north)'][ilt,:,ialt]
        w = g['V!Dn!N (up)'][ilt,:,ialt]

    fig = plt.figure(figsize=[5.6, 8.7])
    hc = plt.contourf(
        lat, alt, variable, 21 if levels is None else levels, cmap=cmap)
    plt.colorbar(hc)
    if vector:
        hq = plt.quiver(lat, alt, ue, w*100, scale=100, scale_units='inches')
    plt.xlabel('Latitude')
    plt.ylabel('Altitude (km)')
    plt.title(var)
    return fig, hc

def height_profile_all():
    for lt in [0, 12]:
        fig, hc = height_profile(
            gall, var='Rho', whichlt=lt, altmin=100, altmax=150, log=True,
            vector=True, levels=None)
        plt.savefig(savepath+'004_height_profile_rho_%02d.png'%lt)
        fig, hc = height_profile(
            gall, var='Temperature', whichlt=lt, altmin=100, altmax=150, log=False,
            vector=True, levels=None)
        plt.savefig(savepath+'004_height_profile_temperature_%02d.png'%lt)
        fig, hc = height_profile(
            gall, var='V!Dn!N (up)', whichlt=lt, altmin=100, altmax=150, log=False,
            vector=False, levels=np.linspace(-0.5,0.5,21), cmap='seismic')
        plt.savefig(savepath+'004_height_profile_vup_%02d.png'%lt)


if __name__ == '__main__':
    plt.close('all')
    multiple_alts_const_alt()
    height_profile_all()


