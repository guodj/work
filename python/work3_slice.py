import gitm_new as gitm
import numpy as np
from numpy import pi
import gitm_3D_const_alt as g3ca
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy import interpolate
plt.close('all')
path1 = '/home/guod/simulation_output/momentum_analysis/'+\
        'run_shrink_iondrift_4_c1/data/3DALL_t030322_060000.bin'
path2 = '/home/guod/simulation_output/momentum_analysis/'+\
        'run_no_shrink_iondrift_4_1/data/3DALL_t030322_060000.bin'
g1 = gitm.read(path1)
g2 = gitm.read(path2)
lat1, lon1, alt1, rho1 = \
    g1['Latitude'][2:-2, 2:-2, 2:-2], g1['Longitude'][2:-2, 2:-2, 2:-2], \
    g1['Altitude'][2:-2, 2:-2, 2:-2], g1['Rho'][2:-2, 2:-2, 2:-2]
lat2, lon2, alt2, rho2 = \
    g2['Latitude'][2:-2, 2:-2, 2:-2], g2['Longitude'][2:-2, 2:-2, 2:-2], \
    g2['Altitude'][2:-2, 2:-2, 2:-2], g2['Rho'][2:-2, 2:-2, 2:-2]

ilat = lat1[0,:,0]<-np.pi/6
lat1, lon1, alt1, rho1, rho2 = \
    lat1[:,ilat,:], lon1[:,ilat,:], alt1[:,ilat,:], rho1[:,ilat,:], \
    rho2[:,ilat,:]

fig = plt.figure(figsize=(5.41, 8.54))
ax = fig.add_subplot(111, projection='3d')
for alt in [150, 200, 300, 400, 500, 600]:
    ialt = np.argmin(np.abs(alt1[0,0,:]/1000-alt))
    sp = interpolate.RectBivariateSpline(
        lon1[:,0,ialt], lat1[0,:,ialt], rho2[:,:,ialt])

    lat11 = np.arange(-(np.pi/2-0.04), -np.pi/6, pi/3/50)
    lon11 = np.arange(0, 2*np.pi, 2*np.pi/300)
    rho11 = sp(lon11, lat11)

    lat11, lon11 = np.meshgrid(lat11, lon11)
    alt11 = np.ones(lat11.shape)*alt1[0,0,ialt]

    r = lat11+np.pi/2
    theta = lon11
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    z = alt11/1000
    v = rho11
    ax.scatter(x, y, z, c=v, s=1, cmap='jet', alpha=0.5)
plt.subplots_adjust(top=1,bottom=0)
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])
ax.set_zlabel('Altitude (km)', labelpad=12)
# ax.set_xlabel('LT:00', labelpad=0)
# ax.set_ylabel('LT:06', labelpad=0)
ax.set_zlim(100,600)
ax.text(ax.get_xlim()[-1]*1.5, 0, 100, 'LT: 06', zdir='y', horizontalalignment='center')
ax.text(0, ax.get_ylim()[0]*1.3, 100, 'LT: 00', zdir='x', horizontalalignment='center')
plt.show()
