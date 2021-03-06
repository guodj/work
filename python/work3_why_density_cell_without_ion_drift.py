"""
The thermospheric neutral density at 150 km is lower at the pole than the
surrounding area even the ion drift is set to zero. Why is this happening?
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import gitm_new as gitm
import gitm_3D_const_alt as g3ca
import gitm_pressure as gp
import calc_rusanov as cr

gallfn = '/Users/guod/data/GITMOutput/run_shrink_notides_tmp/UA/data/3DALL_t030323_000004.bin'
gthmfn = '/Users/guod/data/GITMOutput/run_shrink_notides_tmp/UA/data/3DTHM_t030323_000004.bin'
gall = gitm.read(gallfn)
lontp = gall['Longitude']
lattp = gall['Latitude']
alttp = gall['Altitude']
Re = 6371*1000 # Earth radius, unit: m
RR = Re+alttp
omega = 2*np.pi/(24*3600)
rhotp = gall['Rho']
nwindtp = gall['V!Dn!N (north)']
ewindtp = gall['V!Dn!N (east)'] + omega*RR*np.cos(lattp)
uwindtp = gall['V!Dn!N (up)']
gall['div_rhov'] = \
    -(cr.calc_div_hozt(lontp, lattp, alttp, rhotp*nwindtp, rhotp*ewindtp)\
    +cr.calc_div_vert(alttp, rhotp*uwindtp))/rhotp
gall['vert_divv'] = -cr.calc_div_vert(alttp, uwindtp)
gall['hozt_divv'] = -cr.calc_div_hozt(lontp, lattp, alttp, nwindtp, ewindtp)
gall['vert_vgradrho'] = -uwindtp*cr.calc_rusanov_alts_ausm(alttp,rhotp)/rhotp
gall['hozt_vgradrho'] = \
    -(nwindtp*cr.calc_rusanov_lats(lattp,alttp,rhotp)\
    +ewindtp*cr.calc_rusanov_lons(lontp,lattp,alttp,rhotp))/rhotp
gp.calc_pressure(gall)
gthm = gitm.read(gthmfn)
gthm['total energy'] =\
    gthm['EUV Heating'] + gthm['Conduction'] + \
    gthm['Chemical Heating'] + gthm['Auroral Heating'] + \
    gthm['Joule Heating'] - gthm['NO Cooling'] - gthm['O Cooling']
savepath='/Users/guod/tmp/'

def multiple_alts_const_alt():
    # print('Density and temperature distribution at different heights.')
    # for alt in [100, 105, 110, 115, 120, 125, 135, 145]:
    #     if not os.path.exists(savepath+str(alt)):
    #         os.mkdir(savepath+str(alt))
    #     fig, hc = g3ca.test(
    #         gall, alt=alt, contour=True, zstr='Rho',levels=None, vector=True,
    #         neuion='neu', scale=500, useLT=True)
    #     plt.colorbar(hc)
    #     plt.savefig(savepath+str(alt)+'/001_Rho.png')
    #     fig, hc = g3ca.test(
    #         gall, alt=alt, contour=True, zstr='Temperature',levels=None, vector=True,
    #         neuion='neu', scale=500, useLT=True)
    #     plt.colorbar(hc)
    #     plt.savefig(savepath+str(alt)+'/002_Temperature.png')
    #     fig, hc = g3ca.test(
    #         gall, alt=alt, contour=True, zstr='pressure',levels=None, vector=True,
    #         neuion='neu', scale=500, useLT=True)
    #     plt.colorbar(hc)
    #     plt.savefig(savepath+str(alt)+'/003_Pressure.png')
    # print('Distributions of energy terms at different heights.')
    # for alt in [100, 105, 110, 115, 120, 125, 135, 145]:
    #     if not os.path.exists(savepath+str(alt)):
    #         os.mkdir(savepath+str(alt))
    #     fig, hc = g3ca.test(
    #         gthm, alt=alt, contour=True, zstr='EUV Heating',
    #         levels=np.linspace(-0.01, 0.01, 21), vector=False)
    #     plt.savefig(savepath+str(alt)+'/004_EUVHeating.png')
    #     fig, hc = g3ca.test(
    #         gthm, alt=alt, contour=True, zstr='Conduction',
    #         levels=np.linspace(-0.01, 0.01, 21), vector=False)
    #     plt.savefig(savepath+str(alt)+'/004_Conduction.png')
    #     fig, hc = g3ca.test(
    #         gthm, alt=alt, contour=True, zstr='Chemical Heating',
    #         levels=np.linspace(-0.01, 0.01, 21), vector=False)
    #     plt.savefig(savepath+str(alt)+'/004_ChemicalHeating.png')
    #     fig, hc = g3ca.test(
    #         gthm, alt=alt, contour=True, zstr='Auroral Heating',
    #         levels=np.linspace(-0.01, 0.01, 21), vector=False)
    #     plt.savefig(savepath+str(alt)+'/004_AuroralHeating.png')
    #     fig, hc = g3ca.test(
    #         gthm, alt=alt, contour=True, zstr='Joule Heating',
    #         levels=np.linspace(-0.01, 0.01, 21), vector=False)
    #     plt.savefig(savepath+str(alt)+'/004_JouleHeating.png')
    #     fig, hc = g3ca.test(
    #         gthm, alt=alt, contour=True, zstr='NO Cooling',
    #         levels=np.linspace(-0.01, 0.01, 21), vector=False)
    #     plt.savefig(savepath+str(alt)+'/004_NOCooling.png')
    #     fig, hc = g3ca.test(
    #         gthm, alt=alt, contour=True, zstr='O Cooling',
    #         levels=np.linspace(-0.01, 0.01, 21), vector=False)
    #     plt.savefig(savepath+str(alt)+'/004_OCooling.png')
    #     fig, hc = g3ca.test(
    #         gthm, alt=alt, contour=True, zstr='total energy',
    #         levels=np.linspace(-0.01, 0.01, 21), vector=False)
    #     plt.savefig(savepath+str(alt)+'/004_TotalEnergy.png')
    print('Distributions of terms in the continuity equation at different heights.')
    for alt in [100, 105, 110, 115, 120, 125, 135, 145]:
        if not os.path.exists(savepath+str(alt)):
            os.mkdir(savepath+str(alt))
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='div_rhov',
            levels=np.linspace(-2e-5, 2e-5, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/005_Div_rhov.png')
        plt.colorbar(hc)
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='vert_divv',
            levels=np.linspace(-2e-5, 2e-5, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/005_vert_divv.png')
        plt.colorbar(hc)
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='hozt_divv',
            levels=np.linspace(-2e-5, 2e-5, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/005_hozt_divv.png')
        plt.colorbar(hc)
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='vert_vgradrho',
            levels=np.linspace(-2e-5, 2e-5, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/005_vert_vgradrho.png')
        plt.colorbar(hc)
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='hozt_vgradrho',
            levels=np.linspace(-2e-5, 2e-5, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/005_hozt_vgradrho.png')
        plt.colorbar(hc)
    print('Distributions of momentum terms at different heights.')
    for alt in [100, 105, 110, 115, 120, 125, 135, 145]:
        if not os.path.exists(savepath+str(alt)):
            os.mkdir(savepath+str(alt))
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='div_rhov',
            levels=np.linspace(-2e-5, 2e-5, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/005_Div_rhov.png')
        plt.colorbar(hc)
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='vert_divv',
            levels=np.linspace(-2e-5, 2e-5, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/005_vert_divv.png')
        plt.colorbar(hc)
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='hozt_divv',
            levels=np.linspace(-2e-5, 2e-5, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/005_hozt_divv.png')
        plt.colorbar(hc)
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='vert_vgradrho',
            levels=np.linspace(-2e-5, 2e-5, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/005_vert_vgradrho.png')
        plt.colorbar(hc)
        fig, hc = g3ca.test(
            gall, alt=alt, contour=True, zstr='hozt_vgradrho',
            levels=np.linspace(-2e-5, 2e-5, 21), vector=False)
        plt.savefig(savepath+str(alt)+'/005_hozt_vgradrho.png')
        plt.colorbar(hc)
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
        plt.savefig(savepath+'005_height_profile_rho_%02d.png'%lt)
        fig, hc = height_profile(
            gall, var='Temperature', whichlt=lt, altmin=100, altmax=150, log=False,
            vector=True, levels=None)
        plt.savefig(savepath+'005_height_profile_temperature_%02d.png'%lt)
        fig, hc = height_profile(
            gall, var='pressure', whichlt=lt, altmin=100, altmax=150, log=True,
            vector=True, levels=None)
        plt.savefig(savepath+'005_height_profile_pressure_%02d.png'%lt)
        fig, hc = height_profile(
            gall, var='V!Dn!N (up)', whichlt=lt, altmin=100, altmax=150, log=False,
            vector=False, levels=np.linspace(-0.5,0.5,21), cmap='seismic')
        plt.savefig(savepath+'005_height_profile_vup_%02d.png'%lt)


if __name__ == '__main__':
    plt.close('all')
    multiple_alts_const_alt()
    #height_profile_all()
