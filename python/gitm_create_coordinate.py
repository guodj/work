import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import gitm
import matplotlib.path as mpath
from matplotlib.offsetbox import AnchoredText
import sys

def create_polar(nrow, ncolumn, n, nlat=90, slat=0, centrallon=0,
                 coastlines=True, lonticklabel=(1, 1, 1, 1),
                 dlat=10, useLT=True, **kw):
    '''
    lonticklabel : ticklabel at `top`, `right`, `bottom`, `left`
    centrallon : central longitude at 6 o'clock direction
    '''
    # In the northern hemisphere, central longitude is at 6 o'clock direction
    # In the southern hemisphere, central longitude is at 12 o'clock direction
    if nlat*slat<0:
        sys.exit('Input of `nlat` or `slat` is wrong!!!')
    if nlat == 90:
        centrallat, boundinglat = nlat, slat
    else:
        centrallat, boundinglat = slat, nlat
    if centrallat == -90:
        centrallon = centrallon+180
    projection=ccrs.AzimuthalEquidistant(
            central_latitude=centrallat, central_longitude=centrallon)
    ax = plt.subplot(nrow, ncolumn, n, projection=projection, **kw)

    if coastlines:
        ax.coastlines(lw=0.5, linestyle='--', resolution='110m')

    # longitude grid and label
    lonticklabel_onoff = {'top':lonticklabel[0], 'right':lonticklabel[1],
                          'bottom':lonticklabel[2], 'left':lonticklabel[3]}
    lonticklabel_pos = {'top': 8, 'right': 6, 'bottom': 9, 'left': 7}
    bbox_to_anchor = {'top': [0.5, 1], 'right': [1, 0.5],
                      'bottom':[0.5, 0], 'left': [0, 0.5]}
    for k in range(4):
        ax.plot([centrallon+90*k, centrallon+90*k], [centrallat, boundinglat],
                 '--', color='gray', transform=ccrs.PlateCarree(), lw=0.5)
    if centrallat==90:
        lontickp = ['bottom', 'right', 'top', 'left']
    else:
        lontickp = ['bottom', 'left', 'top', 'right']
    for k in range(4):
        lonbegin = centrallon if centrallat==90 else centrallon-180
        if lonticklabel_onoff[lontickp[k]]==1:
            if useLT:
                lontick_str = '{:02.0f}'.format(90*k/15)
            else:
                lontick_t = (lonbegin+90*k)%360
                lontick_str = '{:.0f}'.format(
                        lontick_t if lontick_t<=180 else lontick_t-360)
            at = AnchoredText(lontick_str, loc=lonticklabel_pos[lontickp[k]],
                              frameon=False,
                              bbox_to_anchor=bbox_to_anchor[lontickp[k]],
                              prop=dict(weight='normal'), pad=0, borderpad=0.12,
                              bbox_transform=ax.transAxes)
            ax.add_artist(at)

    # latitude grid and label
    for k in np.arange(
            min(boundinglat, centrallat), max(boundinglat, centrallat), dlat):
        lon = np.linspace(0, 360, 1000)
        lat = np.ones(lon.shape)*k
        ax.plot(lon, lat, '--', color='gray', lw=0.5,
                transform=ccrs.PlateCarree())
    at = AnchoredText('{:.0f}\u00B0'.format(boundinglat),
                      loc=2, frameon=False,
                      bbox_to_anchor=[0.5+0.5*np.sin(np.pi/4),
                                      0.5-0.5*np.cos(np.pi/4)],
                      pad=0, borderpad=0.12, prop=dict(weight='bold'),
                      bbox_transform=ax.transAxes)
    ax.add_artist(at)

    # boundary
    x0, y0 = projection.transform_point(
            centrallon, centrallat, ccrs.PlateCarree())
    x1, y1 = projection.transform_point(
            centrallon, boundinglat, ccrs.PlateCarree())
    xy = np.sqrt(x1**2+y1**2)
    ax.set_extent([x0-xy, x0+xy, y0-xy, y0+xy], crs=projection)
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    return ax, projection

def create_rectangular(
        nrow, ncolumn, n, slat=-90, nlat=90, centrallon=0,
        coastlines=True, dlat=30, dlon=90, useLT=True, aspect=1, **kw):

    projection=ccrs.PlateCarree(central_longitude=centrallon)
    ax = plt.subplot(nrow, ncolumn, n, projection=projection, **kw)

    if coastlines:
        ax.coastlines(lw=0.5, linestyle='--', resolution='110m')

    ax.set_aspect(aspect)

    lonticks = np.arange(centrallon-180.0, centrallon+181.0, dlon)
    lonticks[0] = lonticks[0]+0.000001
    lonticks[-1] = lonticks[-1]-0.000001
    latticks = np.arange(-90, 91, dlat)
    ax.set_yticks(latticks, crs=ccrs.PlateCarree())
    ax.set_yticklabels(['{:.0f}\u00B0'.format(k) for k in latticks])
    ax.set_xticks(lonticks, crs=ccrs.PlateCarree())
    if not useLT:
        ax.set_xticklabels('{:.0f}\u00B0'.format(k) for k in lonticks)
    else:
        xticklabel = [
                '{:02.0f}'.format(k)
                for k in (12-centrallon/15+lonticks/180*12)%24]
        ax.set_xticklabels(xticklabel)
    ax.set_extent([-180.01, 180.01, slat, nlat], crs=ccrs.PlateCarree())
    return ax, projection


def create_map(nrow, ncolumn, n, plot_type='polar', nlat=90, slat=-90,
               centrallon=0, coastlines=True, dlat=30, dlon=90, useLT=True,
               aspect=1, lonticklabel=(1, 1, 1, 1)):
    if 'po' in plot_type:
        ax, projection = create_polar(
                nrow, ncolumn, n, nlat=nlat, slat=slat, centrallon=centrallon,
                coastlines=coastlines, lonticklabel=lonticklabel, dlat=dlat,
                useLT=useLT)
    else:
        ax, projection = create_rectangular(
                nrow, ncolumn, n, nlat=nlat, slat=slat, centrallon=centrallon,
                coastlines=coastlines, dlat=dlat, dlon=dlon, useLT=useLT,
                aspect=aspect)
    return ax, projection
if __name__ == '__main__':
    plt.close('all')
    plt.figure(figsize=(10, 10))
    gitmfn='/home/guod/WD4T/gitm/run_imfby/run1c/data/3DALL_t100323_040500.bin'
    g = gitm.GitmBin(gitmfn, varlist=['Rho'])
    gt = g['time']
    centrallon = (12-(gt.hour+gt.minute/60+gt.second/3600))*15
    ax, projection = create_map(
            1, 1, 1, 'polar', nlat=90, slat=0, centrallon=centrallon,
            useLT=True, dlat=30)
    lon = g['Longitude'][2:-2, 0, 0]/np.pi*180
    lat = g['Latitude'][0, 2:-2, 0]/np.pi*180
    rho = g['Rho'][2:-2, 2:-2, 36]
    ax.contourf(np.array(lon), np.array(lat), np.array(rho).T,
                transform=ccrs.PlateCarree())
    plt.show()

