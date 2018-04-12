import gitm_new as gitm
import numpy as np
from numpy import pi
import gitm_3D_const_alt as g3ca
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from aacgmv2 import convert
import datetime as dt
plt.close('all')
path1 = '/home/guod/simulation_output/momentum_analysis/'+\
        'run_shrink_iondrift_4_c1/data/3DALL_t030322_060000.bin'
path2 = '/home/guod/simulation_output/momentum_analysis/'+\
        'run_no_shrink_iondrift_4_1/data/3DALL_t030322_060000.bin'
g1 = gitm.read(path1)
g2 = gitm.read(path2)
fig = plt.figure(figsize=(5.41, 8.54))
ax = fig.add_subplot(111, projection='3d')

levels = [np.linspace(6.1e-10, 17.8e-10,21), np.linspace(1.1e-10, 2.7e-10, 21),
          np.linspace(1.3e-11, 3.9e-11,21), np.linspace(2.6e-12, 9.503e-12,21),
          np.linspace(4.3e-13, 2.71e-12,21), np.linspace(7.4e-14, 8.2e-13,21)]
olat = -40
olatr = np.abs(olat)/180*np.pi

glats, glons = convert(-90,0,0,date=dt.date(2003,3,22), a2g=True)
glats, glons = glats/180*pi, glons/180*pi

artalt = [10, 200, 380, 540, 680, 820]
realalt = [150, 200, 300, 400, 500, 600]
for ik, alt in enumerate([150, 200, 300, 400, 500, 600]):
    lon0, lat0, rho1 = g3ca.contour_data('Rho', g1, alt=alt)
    lon0, lat0, rho2 = g3ca.contour_data('Rho', g2, alt=alt)
    lat0, lon0 = lat0/180*pi, lon0/180*pi
    ilat = lat0[:, 0]<=olat/180*np.pi
    lat0, lon0, rho1, rho2 = \
        lat0[ilat,:], lon0[ilat,:], rho1[ilat,:], rho2[ilat,:]
    r = lat0+pi/2
    theta = lon0
    x = r*np.cos(theta)
    x = -x # for sorthern hemisphere
    y = r*np.sin(theta)
    v = rho2
    print (np.min(v), np.max(v))
    ax.contourf(
        x, y, v, zdir='z', offset=artalt[ik], levels=levels[ik],
        cmap='jet', zorder=100)
    ax.scatter(
        -(glats+pi/2)*np.cos(glons), (glats+pi/2)*np.sin(glons),
        artalt[ik], s=50, c='k', zorder=200)
    for rtick in np.arange(10,50,10)/180*np.pi:
        thetatick = np.arange(0, 2*pi, 0.01)
        xtick = rtick*np.cos(thetatick)
        ytick = rtick*np.sin(thetatick)
        hp = ax.plot(
            xtick,ytick, artalt[ik], color='0.7',linestyle=':',
            zorder = 300, linewidth=0.8)
plt.subplots_adjust(top=1,bottom=0)
ax.set_xlim(-olatr, olatr)
ax.set_ylim(-olatr, olatr)
ax.set_zlim(0,800)
ax.axis('off')

for ik, alt in enumerate(artalt):
    ialt = np.argmin(np.abs(g1['Altitude'][0,0,:]/1000-realalt[ik]))
    ax.plot([-olatr, -olatr+0.1], [-olatr, -olatr], alt, color='k')
    ax.text(-olatr-0.1, -olatr, alt,
            '%.0f km'%(g1['Altitude'][0,0,ialt]/1000),
            verticalalignment='center')
ax.plot([-olatr, -olatr], [-olatr, -olatr], [artalt[0], artalt[-1]], color='k')

ax.text(-olatr-0.3, 0, artalt[-1], '06', horizontalalignment='center',
        zdir='x', zorder=1000)
ax.text(0, -olatr-0.3, artalt[-1], '00', horizontalalignment='center',
        zdir='x', zorder=1000)
ax.text(olatr+0.3, 0, artalt[-1], '18', horizontalalignment='center',
        zdir='x', zorder=1000)
ax.text(0, olatr+0.3, artalt[-1], '12', horizontalalignment='center',
        zdir='x', zorder=1000)
plt.show()
