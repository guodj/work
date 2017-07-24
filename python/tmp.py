import gitm
import gitm_gradient as gg
import gitm_divergence as gd
import gitm_3D_const_alt as g3ca
import gitm_create_coordinate as gcc
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt

Re = 6371*1000
RR = Re+400000
g = gitm.GitmBin('/home/guod/big/raid4/guod/run_imfby/run2c/data/'
                 '3DALL_t100323_013000.bin', varlist=[
                     'Rho', 'V!Dn!N (north)', 'V!Dn!N (east)', 'V!Dn!N (up)'])
gg.calc_gradient(g, 'Rho', 'gradn', component='north')
gg.calc_gradient(g, 'Rho', 'grade', component='east')
gg.calc_gradient(g, 'Rho', 'gradr', component='radial')
gd.calc_divergence(g, neuion='neutral', name='divergence')

lon0, lat0, ewind0 = g3ca.contour_data('V!Dn!N (east)', g, alt=400)
lon0, lat0, nwind0 = g3ca.contour_data('V!Dn!N (north)', g, alt=400)
lon0, lat0, uwind0 = g3ca.contour_data('V!Dn!N (up)', g, alt=400)
ewind0 = ewind0 + ((2*np.pi)/(24*3600))*RR+np.cos(lat0*np.pi/180)
lon0, lat0, gradn0 = g3ca.contour_data('gradn', g, alt=400)
lon0, lat0, grade0 = g3ca.contour_data('grade', g, alt=400)
lon0, lat0, gradr0 = g3ca.contour_data('gradr', g, alt=400)

lon0, lat0, rho0 = g3ca.contour_data('Rho', g, alt=400)
lon0, lat0, div0 = g3ca.contour_data('divergence', g, alt=400)

hadvect = ewind0 * grade0 + nwind0 * gradn0
vadvect = uwind0*gradr0

eee = rho0*div0

plt.close('all')
plt.figure(figsize=(9, 9))
ax, projection = gcc.create_map(
        2, 2, 1, 'polar', nlat=90, slat=50,
        centrallon=g3ca.calculate_centrallon(g, 'polar', useLT=True),
        coastlines=False, dlat=30, useLT=True)
ax.contourf(lon0, lat0, -(hadvect+vadvect+eee), transform=ccrs.PlateCarree(),
            levels = np.linspace(-4e-15, 4e-15, 21), cmap='seismic')

ax, projection = gcc.create_map(
        2, 2, 2, 'polar', nlat=90, slat=50,
        centrallon=g3ca.calculate_centrallon(g, 'polar', useLT=True),
        coastlines=False, dlat=30, useLT=True)
ax.contourf(lon0, lat0, -hadvect, transform=ccrs.PlateCarree(),
            levels = np.linspace(-4e-15, 4e-15, 21), cmap='seismic')

ax, projection = gcc.create_map(
        2, 2, 3, 'polar', nlat=90, slat=50,
        centrallon=g3ca.calculate_centrallon(g, 'polar', useLT=True),
        coastlines=False, dlat=30, useLT=True)
ax.contourf(lon0, lat0, -(vadvect+eee), transform=ccrs.PlateCarree(),
            levels = np.linspace(-4e-15, 4e-15, 21), cmap='seismic')

# ax, projection = gcc.create_map(
#         2, 2, 4, 'polar', nlat=90, slat=50,
#         centrallon=g3ca.calculate_centrallon(g, 'polar', useLT=True),
#         coastlines=False, dlat=30, useLT=True)
# ax.contourf(lon0, lat0, -eee, transform=ccrs.PlateCarree(),
#             levels = np.linspace(-4e-15, 4e-15, 21), cmap='seismic')
plt.show()
