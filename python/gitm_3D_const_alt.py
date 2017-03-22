#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
# By Dongjie Guo, 2016-12-06 10:56, UM
# Comments: Routine to make contour or scatter plot of GITM output on a polar or
#           rectangular coordinates as a function of (magnetic) latitude and
#           (magnetic) local time/longitude.
#
# Include: _contour_data() - obtain 2D latitude, LT/Lon and z data for a
#                  contour or scatter plot. Fixed altitude, specified latitudes.
#          _contour_data_mag() - similar with _contour_data(), but for
#                  geomagnetic coordinates
#          _contour_plot() - contour or scatter plot on a polar or rectangular
#                  coordinates. Fixed altitude, specified latitudes.
#          contour_single() - contour or scatter plot on a polar or rectangular
#                  coordinates. Fixed altitude, specified latitudes. For 1 Gitm
#                  output file.
#          contour_diff() - contour or scatter plot on a polar or rectangular
#                  coordinates. Fixed altitude, specified latitudes. For
#                  difference of 2 Gitm output files.
#          _vector_data() - obtain 2D latitude, local time and neutral or ion
#                  velocity for vector plot. Fixed altitude, specified latitudes
#          _vector_plot() - vector plot of neutral or ion velocity on a polar
#                  or rectangular coordinates
#          vector_single() - vector plot of one gitm 3D file
#          vector_diff() - vector plot of the difference of 2 gitm 3D files.
#-------------------------------------------------------------------------------

import gitm
import numpy as np
import matplotlib.pyplot as plt


def _contour_data(
        zkey, gdata, alt=400, nlat=90, slat=-90, ialt=None, useLT=True):
    '''
    Get Latitude, Local Time (or Longitude) and z data for plot
    Input: zkey       = key for z variable (ie 'Vertical TEC')
           gdata      = gitm bin structure
           alt        = altitude (default 400 km)
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)
           ialt       = another way to assign the altitude. if None will use
                        alt, else use ialt
           useLT      = Whether you use LT or Longitude (default LT)

    Output: Lat, LT(Lon), zdata
    '''
    # Find altitude index
    if ialt == None:
        altitude = gdata['Altitude'][0, 0, :]
        ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m
    # Confine latitudes and longitudes
    latitude = gdata['dLat'][0, :, 0]
    ilat = (latitude >= slat) & (latitude <= nlat)
    longitude = gdata['dLon'][:, 0, 0]
    ilon = (longitude>=0) & (longitude<=360)
    # Find the data
    lat0, lon0, lt0, zdata0 = (
            gdata[k][ilon, ...][:, ilat, :][..., ialt].copy()
            for k in ['dLat', 'dLon', 'LT', zkey])
    lon0[lon0>180]=lon0[lon0>180]-360
    lonlt0 = lt0 if useLT else lon0
    # Sort LT or longitude
    ilt = np.argmin(lonlt0[:,0])
    lat0, lonlt0, zdata0 = (
            np.roll(k, -ilt, 0) for k in [lat0, lonlt0, zdata0])
    # Extend Lon/LT (for contour)
    maxlonlt = 24 if useLT else 360
    lat0, lonlt0, zdata0 = (
            np.vstack([k[-1, :][None, :], k]) for k in [lat0, lonlt0, zdata0])
    lat0, lonlt0, zdata0 = (
            np.vstack([k, k[1, :][None, :]]) for k in [lat0, lonlt0, zdata0])
    lonlt0[0, :] = lonlt0[0, :]-maxlonlt
    lonlt0[-1, :] = lonlt0[-1, :]+maxlonlt
    return lat0, lonlt0, zdata0


def _contour_data_mag(
        zkey, gdata, alt=400, nlat=90, slat=-90, ialt=None, useLT=True):
    '''
    Get magnetic Latitude, Local Time and z data for plot
    Input: zkey        = key for z variable (ie 'Vertical TEC')
           gdata       = gitm bin structure
           alt         = altitude (default 400 km)
           nlat        = northern magnetic latitude limit
                         (degrees North, default 90)
           slat        = southern magneticlatitude limit
                         (degrees North, defalut -90)
           ialt        = another way to assign the altitude. if None will use
                         alt, else use ialt
           useLT       = Whether you use LT or Longitude (default LT)

    Output: mLat, mLT, zdata
    '''
    from apexpy import Apex
    from scipy.spatial import cKDTree
    # Find altitude index
    if ialt == None:
        altitude = gdata['Altitude'][0, 0, :]
        ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m
    # Confine latitudes and longitudes
    latitude = gdata['dLat'][0, :, 0]
    ilat = np.argwhere((latitude >= -90) & (latitude <= 90))
    ilatmin, ilatmax = ilat.min(), ilat.max()+1
    longitude = gdata['dLon'][:, 0, 0]
    ilon = np.argwhere((longitude>=0) & (longitude<=360))
    ilonmin, ilonmax = ilon.min(), ilon.max()+1
    lat0 = gdata['dLat'][ilonmin:ilonmax, ilatmin:ilatmax, ialt]
    lon0 = gdata['dLon'][ilonmin:ilonmax, ilatmin:ilatmax, ialt]
    lonlt0 = gdata['LT'][ilonmin:ilonmax, ilatmin:ilatmax, ialt]
    zdata0 = gdata[zkey][ilonmin:ilonmax, ilatmin:ilatmax, ialt]
    # Convert lat, lt to mlat, mlt
    yearf = gdata['time'].year+gdata['time'].month/12
    apx = Apex(date=yearf)
    mlat0, mlt0 = apx.convert(
            lat0, lon0, source='geo', dest='mlt',
            height=gdata['Altitude'][0, 0, ialt]/1000,
            datetime=gdata['time'])
    mlat0, mlon0 = apx.convert(
            lat0, lon0, source='geo', dest='qd',
            height=gdata['Altitude'][0, 0, ialt]/1000,
            datetime=gdata['time'])
    mlonlt0 = mlt0 if useLT else mlon0
    mlat0, mlonlt0 = mlat0.reshape(-1, 1), mlonlt0.reshape(-1, 1)
    zdata0 = zdata0.reshape(-1, 1)
    # Convert spherical coordinate to cartesion coordinate
    maxlonlt = 24 if useLT else 360
    x0 = np.cos(mlat0*np.pi/180)*np.sin(mlonlt0*2*np.pi/maxlonlt)
    y0 = np.cos(mlat0*np.pi/180)*np.cos(mlonlt0*2*np.pi/maxlonlt)
    z0 = np.sin(mlat0*np.pi/180)
    # Create cKDtree
    tree = cKDTree(np.concatenate([x0, y0, z0], axis=1))
    # Create magnetic coordinate grids
    mlat1 = np.arange(slat, nlat+1, 1)
    mlonlt1 = np.arange(0, maxlonlt+0.000001, maxlonlt/180)
    mlat1, mlonlt1 = np.meshgrid(mlat1, mlonlt1)
    # Convert spherical coordinate to cartesion coordinate
    x1 = np.cos(mlat1*np.pi/180)*np.sin(mlonlt1*2*np.pi/maxlonlt)
    y1 = np.cos(mlat1*np.pi/180)*np.cos(mlonlt1*2*np.pi/maxlonlt)
    z1 = np.sin(mlat1*np.pi/180)
    # Get results
    x1, y1, z1 = x1.reshape(-1, 1), y1.reshape(-1, 1), z1.reshape(-1, 1)
    d, izdata1 = tree.query(np.concatenate([x1, y1, z1], axis=1), k=1)
    #w = 1.0/d**2
    #zdata1 = np.sum(
    #        w*(zdata0[izdata1].reshape(w.shape)), axis=1)/np.sum(w, axis=1)
    #zdata1 = zdata1.reshape(mlat1.shape)
    zdata1 = zdata0[izdata1].reshape(mlat1.shape)
    return  mlat1, mlonlt1, zdata1


def _contour_plot(
        ax, plot_type, lat0, lonlt0, zdata0, nlat=90, slat=-90, dlat=10,
        dlonlt=6, lonlt00='S', zmax=None, zmin=None, nzlevels=20, zcolor=None,
        data_type="contour", fill=True, useLT=True, log10=False, *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: ax         = ax to draw the plot on
           plot_type  = key to determine plot type (rectangular, polar)
           lat0       = 2D latitude data
           lonlt0     = 2D local time data
           zdata0     = 2D data to be drawn
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)
           dlat       = increment of latitude ticks (default 10)
           dlonlt     = increment of longitude/local time ticks (default 6)
           lonlt00    = 00 local direction (default 'S')
           zmax       = Maximum z range (default None)
           zmin       = Minimum z range (default None)
           nzlevels   = split [zmin, zmax] into nzlevels sets (default=20)
           zcolor     = Color map for the z variable.  If none, will be chosen
                        based on the z range (default=None)
           data_type  = scatter or contour (default=scatter)
           fill       = whether to use contour fill (default=True)
           useLT      = Whether you use LT or Longitude (default LT)
           log10      = whether the log10(zdata0) is used

    Output: ax = handle of the axis
            h  = handle of contourf or scatter plot
    '''
    if log10:
        zdata0 = np.log10(zdata0)
    maxlonlt = 24 if useLT else 360
    if (zmin is None) | (zmax is None):
        zmin, zmax = np.min(zdata0), np.max(zdata0)
    # set color map
    if zcolor is None:
        if zmin < 0.0 and zmax > 0.0:
            zcolor = "seismic"
            # Change zrange to be centered, making white = 0.0
            if abs(zmin) > zmax:
                zmax = abs(zmin)
            else:
                zmin = -zmax
        else:
            zcolor = "Spectral_r"

    if 'pol' in plot_type.lower(): # polar
        if slat*nlat < 0:
            print('Error: For polar plot, nlat and slat should '
                  'have the same sign')
            return
        csign = 1 if nlat > 0 else -1
        # Convert lat and lt to r and theta
        r0 = 90-csign*lat0
        theta0 = lonlt0*2*np.pi/maxlonlt
        if csign == -1: # south pole
            theta0 = 2*np.pi-theta0
        if 'cont' in data_type.lower(): # contour
            if fill:
                hcont = ax.contourf(
                        np.array(theta0), np.array(r0), np.array(zdata0),
                        levels=np.linspace(zmin, zmax, nzlevels+1),
                        cmap=zcolor, extend='both', *args, **kwargs)
            if not fill:
                if 'colors' in kwargs:
                    zcolor=None
                hcont = ax.contour(
                        np.array(theta0), np.array(r0), np.array(zdata0),
                        levels=np.linspace(zmin, zmax, nzlevels+1),
                        cmap=zcolor, *args, **kwargs)
        if 'sca' in data_type.lower(): # scatter
            hcont = ax.scatter(
                    np.array(theta0), np.array(r0),
                    c=np.array(zdata0), vmin=zmin, vmax=zmax,
                    cmap=zcolor, lw=0, *args, **kwargs)
        # Set polar coordinates
        rticks = np.arange(dlat, r0.max()+dlat/2, dlat)
        rlabels = ['{:02.0f}$^\circ$'.format(k) for k in csign*(90-rticks)]
        thetaticks = np.arange(0, 360, dlonlt*360/maxlonlt)
        thetafmt = '{:02.0f}' if useLT else '{:.0f}'
        thetaticklabel = thetaticks*maxlonlt/360
        if csign==-1:
            thetaticklabel = ((360-thetaticks)*maxlonlt/360)%maxlonlt
        if not useLT:
            thetaticklabel[thetaticklabel>180] -=360
        thetalabels = [thetafmt.format(k) for k in thetaticklabel]
        ax.set_rgrids(rticks, rlabels)
        ax.set_thetagrids(thetaticks, thetalabels)
        ax.set_theta_zero_location(lonlt00)
        ax.set_rlim(0, min(90-slat, 90+nlat))
        return ax, hcont
    if 'rec' in plot_type.lower(): # rectangular
        if 'cont' in data_type.lower(): # contour
            if fill:
                hcont = ax.contourf(
                        np.array(lonlt0), np.array(lat0), np.array(zdata0),
                        levels=np.linspace(zmin, zmax, nzlevels+1),
                        cmap=zcolor, extend='both', *args, **kwargs)
            if not fill:
                if 'colors' in kwargs:
                    zcolor=None
                hcont = ax.contour(
                        np.array(lonlt0), np.array(lat0), np.array(zdata0),
                        levels=np.linspace(zmin, zmax, nzlevels+1),
                        cmap=zcolor, *args, **kwargs)
        if 'sca' in data_type.lower(): # scatter
            hcont = ax.scatter(np.array(lonlt0), np.array(lat0),
                               c=np.array(zdata0), vmin=zmin, vmax=zmax,
                               cmap=zcolor, lw=0, *args, **kwargs)
        if useLT:
            xticks = np.arange(0, maxlonlt+dlonlt/2, dlonlt)
        else:
            xticks = np.arange(-180, 180+dlonlt/2, dlonlt)
        xticklabels = ['{:.0f}'.format(k) for k in xticks]
        if useLT:
            xticklabels = ['{:02.0f}'.format(k) for k in xticks]
        yticks = np.arange(slat, nlat+0.00001, dlat)
        yticklabels = ['{:02.0f}'.format(k) for k in yticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        if useLT:
            ax.set_xlim(0,24)
        else:
            ax.set_xlim(-180,180)
        ax.set_ylim(slat, nlat)
        return ax, hcont


def contour_single(
        ax, zkey, plot_type, gdata, alt=400, title=False,
        nlat=90, slat=-90, dlat=10, dlonlt=6, mag=False,
        lonlt00='S', zmax=None, zmin=None, nzlevels=20, zcolor=None,
        data_type="contour", fill=True, ialt=None, useLT=True, log10=False,
        *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range at a fixed altitude for one 3D file.
    Input: ax         = ax to draw the plot on
           zkey       = key for z variable (ie 'Vertical TEC')
           plot_type  = key to determine plot type (rectangular, polar)
           gdata      = gitm bin structure
           alt        = altitude (km, default 400 km)
           title      = whether to use default title (default False)
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)
           dlat       = increment of latitude ticks (default 10)
           dlonlt     = increment of local time ticks (default 6)
           mag        = whether use magnetic coordinates (default False)
           lonlt00    = 00 local direction (default 'S')
           zmax       = Maximum z range (default None)
           zmin       = Minimum z range (default None)
           nzlevels   = split [zmin, zmax] into nzlevels sets (default=20)
           zcolor     = Color map for the z variable.  If none, will be chosen
                        based on the z range (default=None)
           data_type  = scatter or contour (default=contour)
           fill       = whether to use contour fill (default=True)
           ialt       = another way to assign the altitude. if None will use
                        alt, else use ialt
           useLT      = Whether you use LT or Longitude (default LT)
           log10      = whether the log10 of data is used

    Output: ax = handle of the axis
            h  = handle of contourf or scatter plot
    '''
    from apexpy import Apex
    if not mag:
        lat0, lonlt0, zdata0 = _contour_data(
                zkey=zkey, gdata=gdata, alt=alt, nlat=nlat, slat=slat,
                ialt=ialt, useLT=useLT)
    if mag:
        lat0, lonlt0, zdata0 = _contour_data_mag(
                zkey=zkey, gdata=gdata, alt=alt, nlat=nlat, slat=slat,
                ialt=ialt, useLT=useLT)
    ax, hc = _contour_plot(
             ax=ax, plot_type=plot_type, lat0=lat0, lonlt0=lonlt0, zdata0=zdata0,
             nlat=nlat, slat=slat, dlat=dlat, dlonlt=dlonlt, lonlt00=lonlt00,
             zmax=zmax, zmin=zmin, nzlevels=nzlevels, zcolor=zcolor,
             data_type=data_type, fill=fill, useLT=useLT, log10=log10,
             *args, **kwargs)
    if title:
        # Find altitude index
        if ialt == None:
            altitude = gdata['Altitude'][0, 0, :]
            ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m
        altalt = altitude[ialt]
        ax.set_title(gdata['time'].strftime('%y-%m-%d %H:%M:%S')+
                     ' @ '+'%5.1f'%(altalt/1000)+' km')
    return ax, hc


def contour_diff(
        ax, zkey, plot_type, gdata1, gdata2, alt=400, diff_type='relative',
        title=False,  nlat=90, slat=-90, dlat=10, dlonlt=6, mag=False,
        lonlt00='S', zmax=None, zmin=None, nzlevels=20,zcolor=None,
        data_type="contour", fill=True, ialt=None, useLT=True, log10=False,
        *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range at a fixed altitude for one 3D file.
    Input: ax         = ax to draw the plot on
           zkey       = key for z variable (ie 'Vertical TEC')
           plot_type  = key to determine plot type (rectangular, polar)
           gdata1     = gitm bin structure 1
           gdata2     = gitm bin structure 2
           alt        = altitude (km, default 400 km)
           diff_type  = relative or absolute density difference
                        (default relative)
           title      = whether to use default title (default False)
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)
           dlat       = increment of latitude ticks (default 10)
           dlonlt     = increment of local time ticks (default 6)
           mag        = whether use magnetic coordinates (default False)
           lonlt00    = 00 local direction (default 'S')
           zmax       = Maximum z range (default None)
           zmin       = Minimum z range (default None)
           nzlevels   = split [zmin, zmax] into nzlevels sets (default=20)
           zcolor     = Color map for the z variable.  If none, will be chosen
                        based on the z range (default=None)
           data_type  = scatter or contour (default=scatter)
           fill       = whether to use contour fill (default=True)
           ialt       = another way to assign the altitude. if None will use
                        alt, else use ialt
           useLT      = Whether you use LT or Longitude (default LT)
           log10      = whether the log10 of data is used

    Output: ax = handle of the axis
            h  = handle of contourf or scatter plot
    '''
    if not mag:
        lat1, lonlt1, zdata1 = _contour_data(
                zkey=zkey, gdata=gdata1, alt=alt, nlat=nlat, slat=slat,
                ialt=ialt, useLT=useLT)
        lat2, lonlt2, zdata2 = _contour_data(
                zkey=zkey, gdata=gdata2, alt=alt, nlat=nlat, slat=slat,
                ialt=ialt, useLT=useLT)
    if mag:
        lat1, lonlt1, zdata1 = _contour_data_mag(
                zkey=zkey, gdata=gdata1, alt=alt, nlat=nlat, slat=slat,
                ialt=ialt, useLT=useLT)
        lat2, lonlt2, zdata2 = _contour_data_mag(
                zkey=zkey, gdata=gdata2, alt=alt, nlat=nlat, slat=slat,
                ialt=ialt, useLT=useLT)
    # Make sure that lat1, lonlt1 are the same with lat2, lonlt2
    print('Contour delta lat max: {:5.2f}, delta lt max: {:5.2f}'.format(
            np.abs(lat1-lat2).max(), (np.abs(lonlt1-lonlt2).max())))
    if (np.abs(lat1-lat2).max() > 0.5) | (np.abs(lonlt1-lonlt2).max() > 0.1):
        print('The latitudes or local times obtained from the '
              '2 files are different.')
        return
    if 'rel' in diff_type.lower():
        zdata = 100*(zdata2-zdata1)/zdata1
    if 'abs' in diff_type.lower():
        zdata = zdata2-zdata1
    ax, hc =_contour_plot(
            ax=ax, plot_type=plot_type, lat0=lat1, lonlt0=lonlt1, zdata0=zdata,
            nlat=nlat, slat=slat, dlat=dlat, dlonlt=dlonlt, lonlt00=lonlt00,
            zmax=zmax, zmin=zmin, nzlevels=nzlevels, zcolor=zcolor,
            data_type=data_type, fill=fill, useLT=useLT, log10=log10, *args, **kwargs)
    if title:
        # Find altitude index
        if ialt == None:
            altitude = gdata1['Altitude'][0, 0, :]
            ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m
        altalt = altitude[ialt]
        ax.set_title(
                gdata2['time'].strftime('%m-%d %H:%M')+' - '+
                gdata1['time'].strftime('%m-%d %H:%M')+' @ '+
                '%5.1f'%(altalt/1000)+' km')
    return ax, hc


def _vector_data(
        species, gdata, alt=400, nlat=90, slat=-90, ialt=None, useLT=True):
    '''
    Obtain data for _vector_plot.
    Input: species    = 'neutral' or 'ion'
           gdata      = gitm bin structure
           alt        = altitude (default 400 km)
           nlat       = northern latitude limit (degrees North, default 90)
           slat       = southern latitude limit (degrees North, defalut -90)
           ialt       = another way to assign the altitude. if None will use
                        alt, else use ialt
           useLT      = Whether you use LT or Longitude (default LT)

    Output: lat, lt, nwind, ewind
    '''
    # Find altitude index
    if ialt == None:
        altitude = gdata['Altitude'][0, 0, :]
        ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m
    # Confine latitudes and longitudes
    latitude = gdata['dLat'][0, :, 0]
    splat = int(len(latitude)/36)
    ilat = np.argwhere((latitude >= slat) & (latitude <= nlat))
    ilatmin, ilatmax = ilat.min(), ilat.max()+1
    longitude = gdata['dLon'][:, 0, 0]
    splon = int(len(longitude)/36)
    ilon = np.argwhere((longitude>=0) & (longitude<=360))
    ilonmin, ilonmax = ilon.min(), ilon.max()+1
    # Find the data
    lat0 = gdata['dLat'][ilonmin:ilonmax:splon, ilatmin:ilatmax:splat, ialt]
    lonlt0 = gdata['LT'][ilonmin:ilonmax:splon, ilatmin:ilatmax:splat, ialt]
    if not useLT:
        lonlt0 = gdata['dLon'][ilonmin:ilonmax:splon,
                               ilatmin:ilatmax:splat, ialt]
        lonlt0[lonlt0>180]=lonlt0[lonlt0>180]-360
    if 'neu' in species.lower():
        nwind = gdata['V!Dn!N (north)']
        ewind = gdata['V!Dn!N (east)']
    elif 'ion' in species.lower():
        nwind = gdata['V!Di!N (north)']
        ewind = gdata['V!Di!N (east)']
    nwind = nwind[ilonmin:ilonmax:splon, ilatmin:ilatmax:splat, ialt]
    ewind = ewind[ilonmin:ilonmax:splon, ilatmin:ilatmax:splat, ialt]
    # Sort LT
    if useLT:
        ilt = np.argmin(lonlt0[:,0]) # LT range [0, 24]
        lat0 = np.roll(lat0, -ilt, 0)
        lonlt0 = np.roll(lonlt0, -ilt, 0)
        nwind = np.roll(nwind, -ilt, 0)
        ewind = np.roll(ewind, -ilt, 0)
    return lat0, lonlt0, nwind, ewind


def _vector_plot(
        ax, lat, lonlt, nwind, ewind, plot_type, nlat=90, slat=-90,
        dlat=10, dlonlt=6, lonlt00='S', scale=None, scale_units='inches',
        useLT=True, *args, **kwargs):
    '''
    Creates a rectangular or polar map projection vector plot for a specified
    latitude range.
    Input: ax          = ax to draw the plot on
           lat         = 2D latitudes
           lonlt       = 2D local time
           nwind       = northward wind
           ewind       = eastward wind
           plot_type   = key to determine plot type (rectangular, polar)
           nlat        = northern latitude limit (degrees North, default 90)
           slat        = southern latitude limit (degrees North, defalut -90)
           dlat        = increment of latitude ticks (default 10)
           dlonlt      = increment of local time ticks (default 6)
           lonlt00     = 00 local time/longitude direction (default 'S')
           scale       = the same as in plt.quiver function
           scale_units = the same as in plt.quiver function
           useLT       = Whether you use LT or Longitude (default LT)

    Output: ax = handle of the axis
            h  = handle of contourf or scatter plot
    '''
    #------------------------------------------------------------
    # For now, 00 local time must be southward, i.e. lonlt00='S' is required.
    #------------------------------------------------------------

    maxlonlt = 24 if useLT else 360
    if 'pol' in plot_type.lower(): # polar
        if slat*nlat < 0:
            print('Error: For polar plot, nlat and slat should '
                  'have the same sign')
            return
        csign = 1 if nlat > 0 else -1
        theta0 = lonlt*2*np.pi/maxlonlt
        if csign == -1:
            theta0 = 2*np.pi - theta0
        r0 = 90-csign*lat
        xwind = csign*(ewind*np.cos(theta0)-nwind*np.sin(theta0))
        ywind = csign*(ewind*np.sin(theta0)+nwind*np.cos(theta0))
        hq =ax.quiver(
                theta0, r0, xwind, ywind,
                scale=scale, scale_units=scale_units, *args, **kwargs)
        # Set polar coordinates
        rticks = np.arange(dlat, r0.max()+dlat/2, dlat)
        rlabels = ['{:02.0f}$^\circ$'.format(k) for k in csign*(90-rticks)]
        thetaticks = np.arange(0, 360, dlonlt*360/maxlonlt)
        thetafmt = '{:02.0f}' if useLT else '{:.0f}'
        thetaticklabel = thetaticks*maxlonlt/360
        if csign==-1:
            thetaticklabel = ((360-thetaticks)*maxlonlt/360)%maxlonlt
        if not useLT:
            thetaticklabel[thetaticklabel>180] -=360
        thetalabels = [thetafmt.format(k) for k in thetaticklabel]
        ax.set_rgrids(rticks, rlabels)
        ax.set_thetagrids(thetaticks, thetalabels)
        ax.set_theta_zero_location(lonlt00)
        ax.set_rlim(0, min(90-slat, 90+nlat))
        return ax, hq
    if 'rec' in plot_type.lower(): # rectangular
        hq = ax.quiver(
                lonlt, lat, ewind, nwind, scale=scale, scale_units=scale_units,
                *args, **kwargs)
        if useLT:
            xticks = np.arange(0, maxlonlt+dlonlt/2, dlonlt)
        else:
            xticks = np.arange(-180, 180+dlonlt/2, dlonlt)
        xticklabels = ['{:.0f}'.format(k) for k in xticks]
        if useLT:
            xticklabels = ['{:02.0f}'.format(k) for k in xticks]
        yticks = np.arange(slat, nlat+0.00001, dlat)
        yticklabels = ['{:02.0f}'.format(k) for k in yticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        if useLT:
            ax.set_xlim(0,24)
        else:
            ax.set_xlim(-180,180)
        ax.set_ylim(slat, nlat)
        return ax, hq


def vector_single(
        ax, gdata, species, plot_type, alt=400, title=False,
        nlat=90, slat=-90, dlat=10, dlonlt=6, lonlt00='S', scale=None,
        scale_units='inches', ialt=None, useLT=True, *args, **kwargs):
    '''
    Vector plot for one gitm 3D file
    Input: ax          = ax to draw the plot on
           gdata       = gitm bin structure
           species     = 'neutral' or 'ion'
           plot_type   = key to determine plot type (rectangular, polar)
           alt         = altitude (default 400 km)
           title       = whether to add a default title (default False
           nlat        = northern latitude limit (degrees North, default 90)
           slat        = southern latitude limit (degrees North, defalut -90)
           dlat        = increment of latitude ticks (default 10)
           dlonlt      = increment of local time ticks (default 6)
           lonlt00     = 00 local direction (default 'S')
           scale       = the same as in plt.quiver function
           scale_units = the same as in plt.quiver function
           ialt        = another way to assign the altitude. if None will use
                         alt, else use ialt
           useLT       = Whether you use LT or Longitude (default LT)

    Output: ax = handle of the axis
            hq = handle of the vector plot
    '''
    lat, lonlt, nwind, ewind = _vector_data(
            species=species, gdata=gdata, alt=alt, nlat=nlat, slat=slat,
            ialt=ialt, useLT=useLT)
    ax, hq = _vector_plot(
            ax=ax, lat=lat, lonlt=lonlt, nwind=nwind, ewind=ewind,
            plot_type=plot_type, nlat=nlat, slat=slat, dlat=dlat, dlonlt=dlonlt,
            lonlt00=lonlt00, scale=scale, scale_units=scale_units, useLT=useLT,
            *args, **kwargs)
    if title:
        # Find altitude index
        if ialt == None:
            altitude = gdata['Altitude'][0, 0, :]
            ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m
        altalt = altitude[ialt]
        ax.set_title(gdata['time'].strftime('%y-%m-%d %H:%M:%S')+
                     ' @ '+'%5.1f'%(altalt/1000)+' km')
    return ax, hq


def vector_diff(
        ax, g1, g2, species, plot_type, alt=400, title=False,
        nlat=90, slat=-90, dlat=10, dlonlt=6, lonlt00='S',
        scale=None, scale_units='inches', ialt=None, useLT=True,
        *args, **kwargs):
    '''
    Vector plot for the difference of two gitm 3D files
    Input: ax          = ax to draw the plot on
           g1          = the first gitm 3D file
           g2          = the second gitm 3D file
           species     = 'neutral' or 'ion'
           plot_type   = key to determine plot type (rectangular, polar)
           alt         = altitude (default 400 km)
           title       = whether to add a default title (default False
           nlat        = northern latitude limit (degrees North, default 90)
           slat        = southern latitude limit (degrees North, defalut -90)
           dlat        = increment of latitude ticks (default 10)
           dlonlt      = increment of local time ticks (default 6)
           lonlt00     = 00 local direction (default 'S')
           scale       = the same as in plt.quiver function
           scale_units = the same as in plt.quiver function
           ialt        = another way to assign the altitude. if None will use
                         alt, else use ialt
           useLT       = Whether you use LT or Longitude (default LT)

    Output: ax = handle of the axis
            hq = handle of the vector plot
    '''
    lat1, lonlt1, nwind1, ewind1 = _vector_data(
            species=species, gdata=g1, alt=alt, nlat=nlat, slat=slat,
            ialt=ialt, useLT=useLT)
    lat2, lonlt2, nwind2, ewind2 = _vector_data(
            species=species, gdata=g2, alt=alt, nlat=nlat, slat=slat,
            ialt=ialt, useLT=useLT)
    print('Vector delta lat max: {:5.2f}, delta lonlt max: {:5.2f}'.format(
            np.abs(lat1-lat2).max(), (np.abs(lonlt1-lonlt2).max())))
    if (np.abs(lat2-lat1).max() > 0.5) | (np.abs(lonlt1-lonlt2).max() > 0.1):
        print('The latitudes or local times obtained from the '
              '2 files are different.')
        return
    ax, hq = _vector_plot(
            ax=ax, lat=lat1, lonlt=lonlt1, nwind=nwind2-nwind1,
            ewind=ewind2-ewind1, plot_type=plot_type, nlat=nlat,
            slat=slat, dlat=dlat, dlonlt=dlonlt, lonlt00=lonlt00,
            scale=scale, scale_units=scale_units, useLT=useLT,
            *args, **kwargs)
    if title:
        # Find altitude index
        if ialt == None:
            altitude = g1['Altitude'][0, 0, :]
            ialt = np.argmin(abs(altitude-alt*1000)) # in GITM, unit of alt is m
        altalt = altitude[ialt]
        ax.set_title(g2['time'].strftime('%m-%d %H:%M')+' - '+
                     g1['time'].strftime('%m-%d %H:%M')+' @ '+
                     '%5.1f'%(altalt/1000)+' km')
    return ax, hq


#END
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    # test
    import gitm
    import pandas as pd
    path = '/home/guod/WD2T/run_imfby/'
    g = gitm.GitmBin(
            path+'run2/data/3DALL_t100323_000000.bin',
            varlist=['Rho', 'V!Dn!N (east)', 'V!Dn!N (north)'])
    # Tested parameters
    polar = False
    nlat, slat = 90, -90
    useLT, dlonlt = False , 90
    #-----------------------
    pr = 'pol' if polar else 'rec'
    ax = plt.subplot(polar=polar)
    ax, hc = contour_single(
            ax, 'Rho', pr, g, alt=400,
            nlat=nlat, slat=slat,
            dlonlt=dlonlt, useLT=useLT)
    ax, hc = vector_single(
            ax, g, 'neu', pr, alt=400,
            nlat=nlat, slat=slat,
            dlonlt=dlonlt, useLT=useLT)
    plt.show()
