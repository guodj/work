
if (n_elements(filelist) eq 0) then filelist = findfile("-t *.bin")
filelist = ask('filename to plot',filelist(0))
filelist = findfile(filelist)

nfiles = n_elements(filelist)

if (nFiles eq 1 and strlen(filelist(0)) eq 0) then begin
   print, 'can not find file!'
   stop
endif

if (n_elements(psfile) eq 0) then psfile = filelist(0)+'.ps'
psfile = ask('ps file name',psfile)

file = filelist(0)
VarsToGet=[0,1,2]
gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                        VarsToGet = VarsToGet

alts = reform(data(2,0,0,*)) / 1000.0
lats = reform(data(1,*,*,0))/!dtor
lons = reform(data(0,*,*,0))/!dtor
lon1d = reform(lons(*,0))
lat1d = reform(lats(0,*))
nAlts = n_elements(alts)
nLons = n_elements(lons(*,0))
nLats = n_elements(lats(0,*))

display, lon1d
if (n_elements(iLon) eq 0) then iLon = (nLons*3)/4
iLon = fix(ask('lon to plot',tostr(iLon)))

AllValues = fltarr(nFiles, nLons, nLats)
AllTimes = dblarr(nFiles)
VarsToGet=[33]

for iFile = 0, nFiles-1 do begin

   filename = filelist(iFile)

   print, 'Reading file ',filename

   file = filelist(iFile)
   gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                       VarsToGet = VarsToGet
   AllTimes(iFile) = time

   tec = fltarr(nLons,nLats)

   gitm_read_bin_1var, file, He, timeHe, nVars, Vars, version, $
                       VarsToGet = [31,32]

   for iAlt = 2,nAlts-3 do begin
      tec = tec + (data(0,*,*,iAlt)-He(0,*,*,iAlt)-He(1,*,*,iAlt)) * $
            (alts(iAlt)-alts(iAlt-1))*1000.0
   endfor
   AllValues(iFile,*,*) = tec/1e16

endfor


determine_min_max, AllValues, mini, maxi
nLevels = 31
levels = findgen(nlevels)*(maxi-mini)/(nlevels-1) + mini
maxrange = 90.0

setdevice, 'tec_all_90.ps','p',5

makect,'mid'

ppp = nFiles*4
space = 0.01
pos_space, ppp, space, sizes, nx = 4

plotdumb

for iT = 0, nFiles-1 do begin

   off = iT*4
   get_position, ppp, space, sizes, off, pos

;--------------------------------------------------------------
; Figure out where we are on the page, and whether we need to
; labels or not for the MLT grid
;--------------------------------------------------------------

   no00 = 1
   no06 = 1
   no12 = 1
   no18 = 0

   if (iT eq 0) then n12 = 0
   if (iT eq nFiles-1) then n00 = 0

   c_r_to_a, itime, alltimes(iT)
   c_a_to_s, itime, stime

   ut = itime(3) + itime(4)/60.0 + itime(5)/3600.0
   utrot = ut * 15.0

   print, utrot
  
   lon = lon1d+utrot
   newvalue = reform(allvalues(iT,*,*))
   contour_circle, newvalue, lon, lat1d, $
                   mini = mini, maxi = maxi, $
                   nLevels = nLevels, $
                   no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                   pos = pos, $
                   maxrange = maxrange, /nolines, /nomaxmin

   xc = (pos(2)+pos(0))/2
   xr = (pos(2)-pos(0))/2 * 1.01
   yc = (pos(3)+pos(1))/2
   yr = (pos(3)-pos(1))/2 * 1.01
   xp = xc - xr*sin(!pi/4)
   yp = yc + yr*sin(!pi/4)
   xyouts, xp, yp, 'North', $
           /norm, charsize = 0.8, align = 0.5, orient = 45

   xyouts, pos(0)-0.02, yc, strmid(stime,10,5)+' UT', $
           /norm, charsize = 0.8, align = 0.5, orient = 90

;--------------------------------------------------------------
; Figure out where on the page we should be
;--------------------------------------------------------------

  get_position, ppp, space, sizes, off+1, pos1
  get_position, ppp, space, sizes, off+2, pos2

  pos = pos1
  pos(2) = pos2(2)-0.04

  !p.position = pos

  utime = itime(3)*3600.0 + $
          itime(4)*60.0 + $
          itime(5)
  utime = utime(0)
  p0lon = utime/3600.0 * 360.0 / 24.0

  if (iT eq 0) then title = strmid(stime,0,9) else title = ' '

  center = (180.0-p0lon+360.0) mod 360
  lon0 = center - 45.0
  lon1 = lon0 + 180.0

  map_set, 0.0, 180.0-p0lon, /noerase, title = title, $
           limit = [-90.0,lon0,90.0,lon1]

  newvalue_limited = newvalue

  l = where(newvalue_limited lt levels(1),c)
  if (c gt 0) then newvalue_limited(l) = levels(1)
  l = where(newvalue_limited gt levels(nLevels-2),c)
  if (c gt 0) then newvalue_limited(l) = levels(nLevels-2)

  l = where(lons gt 0 and lons lt 360 and $
            lats gt -90 and lats lt 90)

  nvl = newvalue_limited(l)
  lo = lons(l)
  la = lats(l)

  contour, nvl, lo, la, $
           /follow, /cell_fill, /over, /irr, $
           levels = levels

  oplot, lon1d(iLon)*[1.0,1.0], [-90,90], thick = 5, color = 0

  map_continents
  map_grid, lats = findgen(19)*10-90, glinethick=3

  ;maxs = 'Max : '+string(max(newvalue),format="(f5.1)")+' m/s'

;  xp = pos(2)
;  yp = pos(1)-yr/10.0
;  xyouts, xp, yp, maxs, charsize = 0.9, align = 1.0, /norm

;  xp = pos(0)
;  yp = pos(3)+yr/20.0
;  xyouts, xp, yp, strmid(stime,0,9), $
;          /norm, charsize = 0.9, align = 0.0
;  
;  xp = pos(2)
;  yp = pos(3)+yr/20.0
;  xyouts, xp, yp, strmid(stime,10,5)+' UT', $
;          /norm, charsize = 0.9, align = 1.0


  get_position, ppp, space, sizes, off+3, pos
  pos(0) = pos(0)-0.04
  pos(2) = pos(2)-0.04

;--------------------------------------------------------------
; Figure out where we are on the page, and whether we need to
; labels or not for the MLT grid
;--------------------------------------------------------------

  no00 = 1
  no06 = 1
  no12 = 1
  no18 = 1

  if (iT eq 0) then n12 = 0
  if (iT eq nFiles-1) then n00 = 0

  contour_circle, newvalue, lon, -lat1d, $
                  mini = mini, maxi = maxi, $
                  nLevels = nLevels, $
                  no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                  pos = pos, $
                  maxrange = maxrange, /nolines, /nomaxmin


   xc = (pos(2)+pos(0))/2
   xr = (pos(2)-pos(0))/2 * 1.01
   yc = (pos(3)+pos(1))/2
   yr = (pos(3)-pos(1))/2 * 1.01
   xp = xc - xr*sin(!pi/4)
   yp = yc + yr*sin(!pi/4)
;   xyouts, xp, yp, strmid(stime,10,5)+' UT', $
;           /norm, charsize = 0.8, align = 0.5, orient = 45

   xp = xc + xr*sin(!pi/4)
   yp = yc + yr*sin(!pi/4)
   xyouts, xp, yp, 'South', $
           /norm, charsize = 0.8, align = 0.5, orient = -45

   colortitle = 'TEC (10!U16!Nm!U-2!N)'
  ctpos = [pos(2)+0.005, pos(1), pos(2)+0.02, pos(3)]
  range = [mini,maxi]
  ncolors = 255
  plotct, ncolors, ctpos, range, colortitle, /right
   

endfor

closedevice

setdevice, psfile, 'p', 5

makect,'bristow'

l = where(lat1d ge -60.0 and lat1d le 60.0)

pos = [0.05, 0.3, 0.95, 0.8]

plot, lat1d(l), allvalues(0,iLon,l), $
      yrange = [0,maxi], xstyle = 1, ystyle = 1, $
      pos = pos, $
      xtitle = 'Latitudes (deg)', $
      ytitle = 'TEC', thick = 4

for it = 1, nFiles-1 do begin
   color = 250.0/nFiles * iT
   oplot, lat1d(l), allvalues(iT,iLon,l), color = color, thick = 5
endfor

colortitle = 'UT (hours)'
mini = (min(alltimes) mod (3600.0*24.0)) / 3600.0
maxi = (max(alltimes) mod (3600.0*24.0)) / 3600.0

ctpos = [pos(2)+0.005, pos(1), pos(2)+0.02, pos(3)]
range = [mini,maxi]
ncolors = 255
plotct, ncolors, ctpos, range, colortitle, /right


closedevice

end
