
filelist1 = findfile('3DALL_t121124_0[23456789]*bin')
filelist2 = findfile('3DALL_t121124_1[012345]*bin')

filelist = [filelist1,filelist2]
nFiles = n_elements(filelist)

lon = 360.0-150.0
minlat = 55.0
maxlat = 75.0
alt = 240.0

VarsToGet = [0,1,2,16,17,18,36,37,38]

for iFile = 0, nFiles-1 do begin

   file = filelist(iFile)
   print, 'Reading ',file
   gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                       VarsToGet = VarsToGet

   if (iFile eq 0) then begin

      lon1d = data(0,*,0,0)/!dtor
      Llon = where(lon1d gt lon)
      iLon = Llon(0)-1
      xLon = (lon-lon1d(iLon))/(lon1d(iLon+1)-lon1d(iLon))

      lat1d = data(1,0,*,0)/!dtor
      iLat = where(lat1d ge minlat and lat1d le maxlat, nLats)

      alt1d = data(2,0,0,*)/1000.0
      iAlt = where(alt1d ge alt)
      iAlt = iAlt(0)

      neutral_north = fltarr(nLats, nFiles)
      neutral_east  = fltarr(nLats, nFiles)
      neutral_up    = fltarr(nLats, nFiles)

      ion_north = fltarr(nLats, nFiles)
      ion_east  = fltarr(nLats, nFiles)
      ion_up    = fltarr(nLats, nFiles)

      time_all = dblarr(nFiles)

   endif

   time_all(iFile) = time

   neutral_east(*,iFile)  = data(3,iLon,iLat,iAlt)
   neutral_north(*,iFile) = data(4,iLon,iLat,iAlt)
   neutral_up(*,iFile)    = data(5,iLon,iLat,iAlt)

   ion_east(*,iFile)  = data(6,iLon,iLat,iAlt)
   ion_north(*,iFile) = data(7,iLon,iLat,iAlt)
   ion_up(*,iFile)    = data(8,iLon,iLat,iAlt)

endfor

stime = min(time_all)
etime = max(time_all)

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

ion_max     = float(fix(max(abs(ion_north))/100.0)+1)*100.0
neutral_max = float(fix(max(abs(neutral_north))/100.0)+1)*100.0

neutral_max = 500.0
ion_max = 2000.0

mini = -neutral_max
maxi = neutral_max*nLats 

setdevice, 'gitm_winds_ns.ps','l',5

makect,'mid'

plot, [btr, etr], [mini, maxi], /nodata, $
      ystyle = 5, xstyle = 1,	  $
      xtickname = xtickname,			$
      xtitle = xtitle,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn, $
      pos = [0.1, 0.1, 0.9, 0.95]

for i=0,nLats-1 do begin

   offset = i*neutral_max

   oplot, [btr, etr], offset+[0,0], linestyle = 1
   oplot, time_all-stime, offset+neutral_north(i,*), thick = 4

   for iTime = 0, nFiles-1 do $
      oplot, time_all(iTime)+[0.0,0.0]-stime, $
             offset+[neutral_north(i,iTime),0.0], $
             thick = 2

   dx = etr/100.0
   xyouts, -dx, offset, string(lat1d(iLat(i)), format='(f4.1)'), alignment = 1.0

endfor

dx = etr/100.0

plots, [etr+dx,etr+dx], [nLats, nLats-1]*neutral_max, thick = 2
xyouts, etr+dx, (nLats-1)*neutral_max-neutral_max/10.0, $
        'Neutrals-> '+tostr(neutral_max)+' (m/s)', orient = 270

mini = -ion_max
maxi = ion_max*nLats

plot, [btr, etr], [mini, maxi], /nodata, $
      ystyle = 5, xstyle = 5, /noerase, $
      pos = [0.1, 0.1, 0.9, 0.95]

for i=0,nLats-1 do begin

   offset = i*ion_max

   oplot, time_all-stime, offset+ion_north(i,*), color = 200, thick = 3

endfor

plots, [etr+dx,etr+dx], [3, 2]*ion_max, thick = 2, color = 200
xyouts, etr+dx, 2*ion_max-ion_max/10.0, $
        'Ions-> '+tostr(ion_max)+' (m/s)', orient = 270

closedevice

; ----------------------------------------------------------------------

setdevice, 'gitm_winds_ew.ps','l',5

makect,'mid'

ion_max     = float(fix(max(abs(ion_east))/100.0)+1)*100.0
neutral_max = float(fix(max(abs(neutral_east))/100.0)+1)*100.0

neutral_max = 500.0
ion_max = 2000.0

mini = -neutral_max
maxi = neutral_max*nLats 

plot, [btr, etr], [mini, maxi], /nodata, $
      ystyle = 5, xstyle = 1,	  $
      xtickname = xtickname,			$
      xtitle = xtitle,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn, $
      pos = [0.1, 0.1, 0.9, 0.95]

for i=0,nLats-1 do begin

   offset = i*neutral_max

   oplot, [btr, etr], offset+[0,0], linestyle = 1
   oplot, time_all-stime, offset+neutral_east(i,*), thick = 4

   for iTime = 0, nFiles-1 do $
      oplot, time_all(iTime)+[0.0,0.0]-stime, $
             offset+[neutral_east(i,iTime),0.0], $
             thick = 2

   dx = etr/100.0
   xyouts, -dx, offset, string(lat1d(iLat(i)), format='(f4.1)'), alignment = 1.0

endfor

dx = etr/100.0

plots, [etr+dx,etr+dx], [nLats, nLats-1]*neutral_max, thick = 2
xyouts, etr+dx, (nLats-1)*neutral_max-neutral_max/10.0, $
        'Neutrals-> '+tostr(neutral_max)+' (m/s)', orient = 270

mini = -ion_max
maxi = ion_max*nLats

plot, [btr, etr], [mini, maxi], /nodata, $
      ystyle = 5, xstyle = 5, /noerase, $
      pos = [0.1, 0.1, 0.9, 0.95]

for i=0,nLats-1 do begin

   offset = i*ion_max

   oplot, time_all-stime, offset+ion_east(i,*), color = 200, thick = 3

endfor

plots, [etr+dx,etr+dx], [3, 2]*ion_max, thick = 2, color = 200
xyouts, etr+dx, 2*ion_max-ion_max/10.0, $
        'Ions-> '+tostr(ion_max)+' (m/s)', orient = 270

closedevice

end


