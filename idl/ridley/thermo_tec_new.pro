
wildcard = ask('files to plot tec (e.g., 3DALL*.bin)','3DALL*.bin')

filelist = findfile(wildcard)
file = filelist(0)

gitm_read_header, file, Vars, Time, nLats, nLons, nAlts, version

;c_r_to_a, itime, time
;print, itime
;print, nLons, nLats, nAlts, version

;display, Vars
IsDone = 0
iVar = [0,1,2,33]

nFiles = n_elements(filelist)

;nFiles = 20

makect,'mid'

totaltec    = fltarr(nFiles)
TECnorth    = fltarr(nFiles)
TECnorthDay = fltarr(nFiles)

TECtest1 = fltarr(nFiles)
TECtest2 = fltarr(nFiles)

TECsouth    = fltarr(nFiles)
TECsouthDay = fltarr(nFiles)
time = dblarr(nFiles)

for iFile = 0, nFiles-1 do begin

   file = filelist(iFile)
   print, file
   gitm_read_bin_1var,file, data, time1, nVars, Vars, version, VarsToGet = iVar

   lats = reform(data(1,*,*,0))
   lons = reform(data(0,*,*,0))

   c_r_to_a, itime, time1
   ut = float(itime(3)) + float(itime(4))/60.0 + float(itime(5))/3600.0

   localtime = (lons/!dtor/15.0 + ut) mod 24.0

   if (iFile eq 0) then begin

      nLons = n_elements(lons(*,0))
      nLats = n_elements(lats(0,*))
      dlon = lons
      dlon(0:nlons-2,*) = lons(1:nLons-1,*) - lons(0:nLons-2,*)
      dlon(nlons-1,*) = dlon(0,*)

      dlat = lats
      dlat(*,1:nlats-2) = (lats(*,2:nLats-1) - lats(*,0:nLats-3))/2.0
      dlat(*,0) = lats(*,1) - lats(*,0)
      dlat(*,nLats-1) = lats(*,nLats-1) - lats(*,nLats-2)
      area = 6372000.0*6372000.0 * cos(lats) * dlon * dlat

      nAlts = n_elements(data(0,0,0,*))
      dAlt = reform(data(2,*,*,*))

      tec = fltarr(nFiles, nLons, nLats)

   endif

   time(iFile) = time1
   for iAlt=0,nAlts-2 do begin
      dAlt(*,*,iAlt) = reform(data(2,*,*,iAlt+1))-reform(data(2,*,*,iAlt))
      tec(iFile,*,*) = tec(iFile,*,*) + dAlt(*,*,iAlt)*data(3,*,*,iAlt)
   endfor

   totaltec(iFile) = total(tec(iFile,2:nLons-3,2:nLats-3)* $
                           area(2:nLons-3,2:nLats-3))/ $
                     1.0e16/total(area(2:nLons-3,2:nLats-3))

   l = where(localtime ge 10.0 and localtime le 14.0 and $
             lats/!dtor ge 30.0 and lats/!dtor lt 80.0, c)
   tecsub = reform(tec(iFile,*,*))
   TECnorthDay(iFile) = total(tecsub(l)*area(l))/1.0e16/total(area(l))

   TECtest1(iFile) = total(lons(l)/!dtor*area(l))/total(area(l))
   TECtest2(iFile) = total(lats(l)/!dtor*area(l))/total(area(l))

   l = where(lats/!dtor ge 50.0,c)
   TECnorth(iFile) = total(tecsub(l)*area(l))/1.0e16/total(area(l))

   l = where(localtime ge 10.0 and localtime le 14.0 and $
             lats/!dtor le -30.0 and lats/!dtor gt -80.0)
   TECsouthDay(iFile) = total(tecsub(l)*area(l))/1.0e16/total(area(l))

   l = where(lats/!dtor le -50.0)
   TECsouth(iFile) = total(tecsub(l)*area(l))/1.0e16/total(area(l))

endfor

etime = max(time)
stime = etime - 24.0*3600.0*3.0

tmptime = time - stime

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

setdevice, 'tec_hemisphere.ps','p',5, 0.9

ppp = 2
space = 0.04
pos_space, ppp, space, sizes, ny = ppp

plotnum = 0

l = where(time ge stime and time le etime)

range = mm([tecnorthday(l)*1.2, tecsouthday(l)*1.2, 20.0])
range(0) = 0.0

get_position, ppp, space, sizes, 0, pos, /rect
pos(0) = pos(0)+0.075

plot, tmptime, tecnorthday, $
      xtickname = xtickname,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn,   $
      pos = pos, /noerase, $
      yrange = range, ystyle = 1, $
      thick = 3, xrange = [btr,etr], xstyle = 1,$
      ytitle = 'Northern Hemisphere TEC'

oplot, tmptime, tecnorth, linestyle = 2, thick = 3

dy = max(range)
dx = max(tmptime)/10.0

oplot, [dx*5,dx*6], [0.9,0.9]*dy, thick = 3
xyouts, dx*6.1, 0.9*dy, '10-14 LT Only (30-80 deg)'
oplot, [dx*5,dx*6], [0.825,0.825]*dy, thick = 3, linestyle = 2
xyouts, dx*6.1, 0.825*dy, 'Whole Polar Region (>50 deg)'

get_position, ppp, space, sizes, 1, pos, /rect
pos(0) = pos(0)+0.075

plot, tmptime, tecsouthday, $
      xtickname = xtickname,			$
      xtitle = xtitle,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn,   $
      pos = pos, /noerase, $
      yrange = range, ystyle = 1, $
      thick = 3, xrange = [btr,etr], xstyle = 1, $
      ytitle = 'Southern Hemisphere TEC'

oplot, tmptime, tecsouth, linestyle = 2, thick = 3

dy = max(range)
dx = max(tmptime)/10.0

oplot, [dx*5,dx*6], [0.9,0.9]*dy, thick = 3
xyouts, dx*6.1, 0.9*dy, '10-14 LT Only (30-80 deg)'
oplot, [dx*5,dx*6], [0.825,0.825]*dy, thick = 3, linestyle = 2
xyouts, dx*6.1, 0.825*dy, 'Whole Polar Region (<-50 deg)'

spawn,'pwd',pwd
xyouts, 1.05,0.0,pwd,/norm,charsize=0.8, orient=90

closedevice

end

;;--------------------------------------------------------------
;; Set up plot sizes for the circles
;;--------------------------------------------------------------
;
;setdevice, 'tec.ps', 'p', 5
;
;ppp = 12
;
;space = 0.005
;pos_space, ppp, space, sizes
;total_plots = sizes.nby*sizes.nbx
;
;;--------------------------------------------------------------
;; determine how much slop we have for the color table.
;;--------------------------------------------------------------
;
;get_position, ppp, space, sizes, 0, pos
;x_left = pos(0)
;get_position, ppp, space, sizes, sizes.nbx-1, pos
;x_right = pos(2)
;slop = 1.0 - (x_right - x_left)
;
;;--------------------------------------------------------------
;; If there is too little room for the color table, we need
;; to make a little more room.
;;--------------------------------------------------------------
;
;if (slop gt 0.02) then begin
;  if (x_right lt 0.98) then shift_left = 0.0 				$
;  else shift_left = 0.02 - (1.0-x_right)
;  ct_left  = x_right+0.01
;  ct_right = ct_left+0.02
;  d_shift_left = 0.0
;endif else begin
;  space = space*2.0
;  pos_space, ppp, space, sizes
;  get_position, ppp, space, sizes, 0, pos
;  x_left = pos(0)
;  get_position, ppp, space, sizes, sizes.nbx-1, pos
;  x_right = pos(2)
;  shift_left = x_left
;  d_shift_left = space*0.5
;  ct_right = 1.00
;  ct_left  = ct_right - shift_left - d_shift_left*(sizes.nbx-1)
;endelse

tec = tec/1.0e16

maxrange = 30.0

lo = reform(lons(*,0))/!dtor
la = reform(lats(0,*))/!dtor

l = where(90-la lt maxrange)
mini = 0.0
maxi = float(fix(max(tec))+1.0)

for iFile = 0, nFiles-1 do begin

   setdevice, 'tec_'+tostr(iFile,3)+'.ps', 'p', 5 

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

   c_r_to_a, itime, time(iFile)
   c_a_to_s, itime, stime

   ut = itime(3) + itime(4)/60.0 + itime(5)/3600.0
   utrot = ut * 15.0

   nLevels = 31
   lon = lo+utrot
   contour_circle, reform(tec(iFile,*,*)), lon, la, $
                   mini = mini, maxi = maxi, $
                   nLevels = nLevels, $
                   no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                   pos = pos, $
                   maxrange = maxrange, title = title 

   xc = (pos(2)+pos(0))/2
   xr = (pos(2)-pos(0))/2 * 1.01
   yc = (pos(3)+pos(1))/2
   yr = (pos(3)-pos(1))/2 * 1.01
   xp = xc - xr*sin(!pi/4)
   yp = yc + yr*sin(!pi/4)
   xyouts, xp, yp, 'North', $
           /norm, charsize = 0.9, align = 0.5, orient = 45

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

   nLevels = 31
   lon = lo+utrot
   contour_circle, reform(tec(iFile,*,*)), lon, -la, $
                   mini = mini, maxi = maxi, $
                   nLevels = nLevels, $
                   no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                   pos = pos, $
                   maxrange = maxrange, title = title 

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
   ;pos([1,3]) = pos([1,3]) - 0.01

   nLevels = 31
   lon = lo+utrot

   !p.position = pos

   utime = itime(3)*3600.0 + $
           itime(4)*60.0 + $
           itime(5)
   utime = utime(0)
   p0lon = utime/3600.0 * 360.0 / 24.0

   map_set, 0.0, 180.0-p0lon, /noerase

   newrat = reform(tec(iFile,1:nLons-2,1:nLats-2))

   newlat = lats(1:nLons-2,1:nLats-2)/!dtor
   newlon = lons(1:nLons-2,1:nLats-2)/!dtor
   nLo  = nLons-2
   nLa  = nLats-2

   newrat(0,*)       = (newrat(1,*)+newrat(nLo-2,*))/2.0
   newrat(nLo-1,*) = (newrat(1,*)+newrat(nLo-2,*))/2.0
   newrat(*,0)       = mean(newrat(*,1))
   newrat(*,nLa-1) = mean(newrat(*,nLa-2))

   newlon(0,*)       = 0.0
   newlon(nLo-1,*) = 360.0
   newlat(*,0) = -90.0
   newlat(*,nLa-1) =  90.0

   levels = findgen(nlevels)*(maxi-mini)/(nlevels-1) + mini

   l = where(newrat lt levels(1),c)
   if (c gt 0) then newrat(l) = levels(1)

   contour, newrat, newlon, newlat, $
            /follow, /cell_fill, /over, $
            levels = levels

   map_continents
   map_grid, lats = findgen(19)*10-90, glinethick=3

   maxs = 'Max : '+string(max(newrat),format="(f5.1)")+'(10!U16!N/m!U2!N)'

   xp = pos(2)
   yp = pos(1)-yr/10.0
   xyouts, xp, yp, maxs, charsize = 0.9, align = 1.0, /norm

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
   units = 'TEC (10!U16!N/m!U2!N)'
   ncolors = 255
   plotct, ncolors, ctpos, range, units, /right

;   contour, tec, data(0,*,*,40)/!dtor, data(1,*,*,40)/!dtor, $
;            levels = findgen(31),/fill
;   
;   contour, tec, data(0,*,*,40)/!dtor, data(1,*,*,40)/!dtor, $
;            levels = findgen(11)*3,/follow, /noerase

  closedevice
   
endfor


end
 
