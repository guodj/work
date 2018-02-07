import matplotlib.pyplot as plt
import gitm
import numpy as np
import seaborn as sns
import calc_rusanov as cr

sns.set('paper', 'whitegrid')
plt.close('all')
plt.figure(figsize=(6.5,7.2))
fn1 = '/home/guod/simulation_output/momentum_analysis/'\
     'run_shrink_iondrift_4_c1/data/3DALL_t030322_060000.bin'
fn2 = '/home/guod/simulation_output/momentum_analysis/'\
     'run_no_shrink_iondrift_4_1/data/3DALL_t030322_060000.bin'
g1 = gitm.GitmBin(fn1)
g2 = gitm.GitmBin(fn2)
g1['diffrho'] = 100*(g2['Rho']-g1['Rho'])/g1['Rho']
color =['b', 'r']
for ig, g in enumerate([g1,g2]):
    lons,lats,alts = (g[k] for k in ['Longitude', 'Latitude', 'Altitude'])
    rho = g['Rho']

    nlat, slat = -30, -90
    for i, alt in enumerate([130,300,600]):
        ax1 = plt.subplot(3,2,i*2+1)
        ax2 = plt.subplot(3,2,i*2+2)
        lt = 6
        ilt = np.argmin(np.abs(g['LT'][2:-2,0,0]-lt))+2
        ialt = np.argmin(np.abs(alts[0,0,2:-2]/1000-alt))+2
        ilat = (lats[0,:,0]/np.pi*180<nlat) &(lats[0,:,0]/np.pi*180>slat)
        ax1.plot(
            180+lats[ilt,ilat,ialt]/np.pi*180,rho[ilt,ilat,ialt], color[ig])
        ax2.plot(
            180+lats[ilt,ilat,ialt]/np.pi*180,
            g1['diffrho'][ilt,ilat,ialt], 'k')
        lat1, rho1, rho11 = \
            180+lats[ilt,2,ialt]/np.pi*180,\
            rho[ilt,2, ialt], g1['diffrho'][ilt,2,ialt]

        lt = 18
        ilt = np.argmin(np.abs(g['LT'][2:-2,0,0]-lt))+2
        ialt = np.argmin(np.abs(g['Altitude'][0,0,2:-2]/1000-alt))+2
        ilat = (lats[0,:,0]/np.pi*180<nlat) &(lats[0,:,0]/np.pi*180>slat)
        ax1.plot(-lats[ilt,ilat,ialt]/np.pi*180,rho[ilt,ilat,ialt], color[ig])
        ax2.plot(
            -lats[ilt,ilat,ialt]/np.pi*180,g1['diffrho'][ilt,ilat,ialt], 'k')
        lat2, rho2, rho22 = \
                -lats[ilt,2,ialt]/np.pi*180,\
                rho[ilt,2, ialt], g1['diffrho'][ilt,2,ialt]

        ax1.plot([lat1,lat2], [rho1,rho2],color[ig])
        ax2.plot([lat1,lat2], [rho11,rho22],'k')

        for iax in [ax1, ax2]:
            plt.sca(iax)
            plt.xlim([-nlat,180+nlat])
            if iax == ax2:
                plt.ylim([-35,5])
            plt.xticks(np.arange(-nlat,181+nlat,30),
                np.concatenate(
                    [np.arange(nlat,-90,-30), np.arange(-90,nlat+1,30)]))
            if i<2:
                plt.xticks(np.arange(-nlat,181+nlat,30),[])
            if i==2:
                ax1.set_xlabel('Latitude')
                ax2.set_xlabel('Latitude')


plt.figure()
lon2 = np.array(g2['Longitude'])
lat2 = np.array(g2['Latitude'])
alt2 = np.array(g2['Altitude'])
rho2 = np.array(g2['Rho'])
gradrho_rho = np.sqrt(cr.calc_rusanov_lats(lat2,alt2,rho2)**2 +
    cr.calc_rusanov_lons(lon2,lat2,alt2,rho2)**2)/rho2
for lt in [0, 3, 6, 9, 12, 15, 18, 21]:
    ilt = np.argmin(np.abs(g2['LT'][2:-2,0,0]-lt))+2
    ilat = np.argmin(np.abs(g2['dLat'][0,2:-2,0]-75))+2
    plt.plot(gradrho_rho[ilt,ilat,2:-2],g2['Altitude'][0,0,2:-2]/1000)
plt.xlabel(r'$\frac{\nabla\rho}{\rho}$ ($m^{-1}$)')
plt.ylabel(r'Altitude (km)')
plt.show()
