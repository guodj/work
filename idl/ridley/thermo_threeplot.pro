
pro thermo_threeplot, psfile, value, time, lons, lats, mini, maxi, $
                      title, colortitle, maxrange, $
                      vn = vn, ve = ve

  nLevels = 31
  if (mini*maxi lt 0.0) then begin
     maxi = max([maxi,abs(mini)])
     mini = -maxi
  endif
  levels = findgen(nlevels)*(maxi-mini)/(nlevels-1) + mini

  nLons = n_elements(value(*,0))
  nLats = n_elements(value(0,*))
  newvalue = reform(value(1:nLons-2,1:nLats-2))

  newlat = lats(1:nLons-2,1:nLats-2)/!dtor
  newlon = lons(1:nLons-2,1:nLats-2)/!dtor
  nLo  = nLons-2
  nLa  = nLats-2

  newvalue(0,*)       = (newvalue(1,*)+newvalue(nLo-2,*))/2.0
  newvalue(nLo-1,*) = (newvalue(1,*)+newvalue(nLo-2,*))/2.0
  newvalue(*,0)       = mean(newvalue(*,1))
  newvalue(*,nLa-1) = mean(newvalue(*,nLa-2))

  newlon(0,*)       = 0.0
  newlon(nLo-1,*) = 360.0
  newlat(*,0) = -90.0
  newlat(*,nLa-1) =  90.0

  setdevice, psfile, 'p', 5 

  if (mini lt 0) then makect,'mid' $
  else makect,'bristow'

  ppp = 4
  space = 0.01
  pos_space, ppp, space, sizes

  plotdumb

;--------------------------------------------------------------
; Figure out where on the page we should be
;--------------------------------------------------------------

  get_position, ppp, space, sizes, 0, pos

;--------------------------------------------------------------
; Figure out where we are on the page, and whether we need to
; labels or not for the MLT grid
;--------------------------------------------------------------

  no00 = 1
  no06 = 1
  no12 = 0
  no18 = 0

  c_r_to_a, itime, time
  c_a_to_s, itime, stime

  ut = itime(3) + itime(4)/60.0 + itime(5)/3600.0
  utrot = ut * 15.0
  
  lon = reform(newlon(*,0))+utrot
  contour_circle, newvalue, lon, reform(newlat(0,*)), $
                  mini = mini, maxi = maxi, $
                  nLevels = nLevels, $
                  no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                  pos = pos, $
                  maxrange = maxrange, /nolines

  if (n_elements(vn) gt 0 and n_elements(ve) gt 0) then begin

     if (n_elements(factor) eq 0) then factor = max(sqrt(ve^2+vn^2))/10.0

     for iLon = 2, nLons-3,2 do for iLat=2,nLats-3,2 do begin
        r = 90.0 - lats(iLon,iLat)/!dtor 
        if (r lt maxrange and r gt 0) then begin
           lo = lons(iLon,iLat) + utrot*!dtor - !pi/2 
           x = r * cos(lo)
           y = r * sin(lo)
           ux = - vn(iLon,iLat) * cos(lo) - ve(iLon,iLat) * sin(lo)
           uy = - vn(iLon,iLat) * sin(lo) + ve(iLon,iLat) * cos(lo)
           ux = ux/factor
           uy = uy/factor
           oplot,[x,x+ux],[y,y+uy], color = 0, thick = 2.0

           u = sqrt(ux^2+uy^2)
           if (u gt 0) then begin
              t = asin(uy/u)
              if (ux lt 0) then t = !pi-t
              t1 = t+!pi/12
              t2 = t-!pi/12
              ux1 = 0.6 * u * cos(t1)
              uy1 = 0.6 * u * sin(t1)
              ux2 = 0.6 * u * cos(t2)
              uy2 = 0.6 * u * sin(t2)
              oplot,[x+ux, x+ux1],[y+uy,y+uy1], color = 0, thick = 2.0
              oplot,[x+ux, x+ux2],[y+uy,y+uy2], color = 0, thick = 2.0
           endif

        endif
     endfor

  endif


  xc = (pos(2)+pos(0))/2
  xr = (pos(2)-pos(0))/2 * 1.01
  yc = (pos(3)+pos(1))/2
  yr = (pos(3)-pos(1))/2 * 1.01
  xp = xc - xr*sin(!pi/4)
  yp = yc + yr*sin(!pi/4)
  xyouts, xp, yp, 'North', $
          /norm, charsize = 0.9, align = 0.5, orient = 45

  xyouts, pos(2)+space/2, pos(3), title, align = 0.5, /norm

;************************************************************************

  get_position, ppp, space, sizes, 1, pos

;--------------------------------------------------------------
; Figure out where we are on the page, and whether we need to
; labels or not for the MLT grid
;--------------------------------------------------------------

  no00 = 1
  no06 = 0
  no12 = 0
  no18 = 1

  contour_circle, newvalue, lon, -reform(newlat(0,*)), $
                  mini = mini, maxi = maxi, $
                  nLevels = nLevels, $
                  no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                  pos = pos, $
                  maxrange = maxrange, /nolines

  if (n_elements(vn) gt 0 and n_elements(ve) gt 0) then begin

     print, 'Plotting Winds (south)!'

     for iLon = 2, nLons-3,2 do for iLat=2,nLats-3,2 do begin
        r = 90.0 + lats(0,iLat)/!dtor 
        if (r lt maxrange and r gt 0) then begin
           lo = lons(iLon,iLat) + utrot*!dtor - !pi/2 
           x = r * cos(lo)
           y = r * sin(lo)
           ux =   vn(iLon,iLat) * cos(lo) - ve(iLon,iLat) * sin(lo)
           uy =   vn(iLon,iLat) * sin(lo) + ve(iLon,iLat) * cos(lo)
           ux = ux/factor
           uy = uy/factor
           oplot,[x,x+ux],[y,y+uy], color = 0, thick = 2.0

           u = sqrt(ux^2+uy^2)
           if (u gt 0) then begin
              t = asin(uy/u)
              if (ux lt 0) then t = !pi-t
              t1 = t+!pi/12
              t2 = t-!pi/12
              ux1 = 0.6 * u * cos(t1)
              uy1 = 0.6 * u * sin(t1)
              ux2 = 0.6 * u * cos(t2)
              uy2 = 0.6 * u * sin(t2)
              oplot,[x+ux, x+ux1],[y+uy,y+uy1], color = 0, thick = 2.0
              oplot,[x+ux, x+ux2],[y+uy,y+uy2], color = 0, thick = 2.0
           endif

        endif
     endfor

  endif

  xc = (pos(2)+pos(0))/2
  xr = (pos(2)-pos(0))/2 * 1.01
  yc = (pos(3)+pos(1))/2
  yr = (pos(3)-pos(1))/2 * 1.01
  xp = xc + xr*sin(!pi/4)
  yp = yc + yr*sin(!pi/4)
  xyouts, xp, yp, 'South', $
          /norm, charsize = 0.9, align = 0.5, orient = -45

;************************************************************************

;--------------------------------------------------------------
; Figure out where on the page we should be
;--------------------------------------------------------------

  get_position, ppp, space, sizes, 2, pos1
  get_position, ppp, space, sizes, 3, pos2

  pos = pos1
  pos(2) = pos2(2)-0.04

  !p.position = pos

  utime = itime(3)*3600.0 + $
          itime(4)*60.0 + $
          itime(5)
  utime = utime(0)
  p0lon = utime/3600.0 * 360.0 / 24.0

  map_set, 0.0, 180.0-p0lon, /noerase

  newvalue_limited = newvalue

  l = where(newvalue_limited lt levels(1),c)
  if (c gt 0) then newvalue_limited(l) = levels(1)
  l = where(newvalue_limited gt levels(nLevels-2),c)
  if (c gt 0) then newvalue_limited(l) = levels(nLevels-2)

  contour, newvalue_limited, newlon, newlat, $
           /follow, /cell_fill, /over, $
           levels = levels

  map_continents
  map_grid, lats = findgen(19)*10-90, glinethick=3

  ;maxs = 'Max : '+string(max(newvalue),format="(f5.1)")+' m/s'

;  xp = pos(2)
;  yp = pos(1)-yr/10.0
;  xyouts, xp, yp, maxs, charsize = 0.9, align = 1.0, /norm

  xp = pos(0)
  yp = pos(3)+yr/20.0
  xyouts, xp, yp, strmid(stime,0,9), $
          /norm, charsize = 0.9, align = 0.0
  
  xp = pos(2)
  yp = pos(3)+yr/20.0
  xyouts, xp, yp, strmid(stime,10,5)+' UT', $
          /norm, charsize = 0.9, align = 1.0

  ctpos = [pos(2)+0.005, pos(1), pos(2)+0.02, pos(3)]
  range = [mini,maxi]
  ncolors = 255
  plotct, ncolors, ctpos, range, colortitle, /right

  closedevice
  
end
