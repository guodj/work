import gitm_new as gitm
import matplotlib.pyplot as plt
import numpy as np
def const_lon(g, whichlt=12, altmin=100, altmax=600):
    ilt = np.argmin(np.abs(g['LT'][2:-2, 0, 0]-whichlt))+2
    ialt = (g['Altitude'][0,0,:]>=altmin*1000) & \
           (g['Altitude'][0,0,:]<=altmax*1000)
    lat = g['dLat'][ilt, :, ialt]
    alt = g['Altitude'][ilt,:,ialt]/1000
    temp = g['Temperature'][ilt,:,ialt]
    rho = np.log10(g['Rho'][ilt,:,ialt])
    ue = g['V!Dn!N (north)'][ilt,:,ialt]
    w = g['V!Dn!N (up)'][ilt,:,ialt]

    plt.figure(figsize=[8.6, 8.7])
    plt.subplot(1,2,1)
    plt.contourf(lat, alt, rho, 21, cmap='jet')
    plt.quiver(lat, alt, ue, w*10, scale=500, scale_units='inches')
    plt.xlabel('Latitude')
    plt.ylabel('Altitude (km)')
    plt.title('log10(Density)')

    plt.subplot(1,2,2)
    plt.contourf(lat, alt, temp, 21, cmap='jet')
    plt.quiver(lat, alt, ue, w*10, scale=500, scale_units='inches')
    plt.xlabel('Latitude')
    plt.title('Temperature')
    plt.show()
    return

def const_lon_diff(g1, g2, whichlt=12, scale=100, altmin=100, altmax=600, latmin=-90, latmax=90):
    ilt = np.argmin(np.abs(g1['LT'][2:-2, 0, 0]-whichlt))+2
    ialt = (g1['Altitude'][0,0,:]>=altmin*1000) & \
           (g1['Altitude'][0,0,:]<=altmax*1000)
    lat = g1['dLat'][ilt, 2:-2, ialt]
    alt = g2['Altitude'][ilt,2:-2,ialt]/1000
    temp = g2['Temperature'][ilt,2:-2,ialt] - g1['Temperature'][ilt,2:-2,ialt]
    rho = 100*(g2['Rho'][ilt,2:-2,ialt] - g1['Rho'][ilt,2:-2,ialt])/g1['Rho'][ilt,2:-2,ialt]
    ue = g2['V!Dn!N (north)'][ilt,2:-2,ialt] - g1['V!Dn!N (north)'][ilt,2:-2,ialt]
    w = g2['V!Dn!N (up)'][ilt,2:-2,ialt] - g1['V!Dn!N (up)'][ilt,2:-2,ialt]

    plt.figure(figsize=[8.6, 8.7])
    plt.subplot(1,2,1)
    plt.contourf(lat, alt, rho, levels=np.linspace(-30,30,21), cmap='seismic')
    plt.quiver(lat, alt, ue, w, scale=scale, scale_units='inches')
    plt.xlim(latmin, latmax)
    plt.ylim(altmin, altmax)
    plt.xlabel('Latitude')
    plt.ylabel('Altitude (km)')
    plt.title('Density')

    plt.subplot(1,2,2)
    plt.contourf(lat, alt, temp, levels=np.linspace(-60,60,21), cmap='seismic')
    plt.quiver(lat, alt, ue, w, scale=scale, scale_units='inches')
    plt.xlim(latmin, latmax)
    plt.ylim(altmin, altmax)
    plt.xlabel('Latitude')
    plt.title('Temperature')
    plt.show()
    return

if __name__=='__main__':
    plt.close('all')
    fn ='/home/guod/tmp/3DALL_t030321_000000.bin'
    fn1 = '/home/guod/simulation_output/momentum_analysis/'\
          'run_shrink_iondrift_4_c1_high_resolution/data/3DALL_t030322_000103.bin'
    fn2 = '/home/guod/simulation_output/momentum_analysis/'\
          'run_no_shrink_iondrift_4_1_high_resolution/data/3DALL_t030322_000103.bin'
    fn1 = '/home/guod/simulation_output/momentum_analysis/'\
          'run_shrink_iondrift_4_c1_high_resolution/data/3DALL_t030322_000203.bin'
    fn2 = '/home/guod/simulation_output/momentum_analysis/'\
          'run_no_shrink_iondrift_4_1_high_resolution/data/3DALL_t030322_000203.bin'
    fn1 = '/home/guod/simulation_output/momentum_analysis/'\
          'run_shrink_iondrift_4_c1_high_resolution/data/3DALL_t030322_000302.bin'
    fn2 = '/home/guod/simulation_output/momentum_analysis/'\
          'run_no_shrink_iondrift_4_1_high_resolution/data/3DALL_t030322_000302.bin'
    fn1 ='/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_4_c1/data/3DALL_t030322_000501.bin'
    fn2 ='/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_4_1/data/3DALL_t030322_000501.bin'
    # fn1 ='/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_4_c1/data/3DALL_t030322_001003.bin'
    # fn2 ='/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_4_1/data/3DALL_t030322_001003.bin'
    fn1 ='/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_4_c1/data/3DALL_t030322_003003.bin'
    fn2 ='/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_4_1/data/3DALL_t030322_003002.bin'
    # fn1 ='/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_4_c1/data/3DALL_t030322_010002.bin'
    # fn2 ='/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_4_1/data/3DALL_t030322_010000.bin'
    # fn1 ='/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_4_c1/data/3DALL_t030322_060000.bin'
    # fn2 ='/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_4_1/data/3DALL_t030322_060000.bin'
    fn1 = '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_4_c4/data/3DALL_t030322_180502.bin'
    fn2 = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_4_4/data/3DALL_t030322_180502.bin'
    fn1 = '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_4_c4/data/3DALL_t030322_181003.bin'
    fn2 = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_4_4/data/3DALL_t030322_181000.bin'
    fn1 = '/home/guod/simulation_output/momentum_analysis/run_shrink_iondrift_4_c4/data/3DALL_t030322_183002.bin'
    fn2 = '/home/guod/simulation_output/momentum_analysis/run_no_shrink_iondrift_4_4/data/3DALL_t030322_183002.bin'
    g1 = gitm.read(fn1)
    g2 = gitm.read(fn2)
    #const_lon(g2, 0, 100, 600)
    const_lon_diff(g1, g2, 4, 50, 100, 600, -90, -30)
