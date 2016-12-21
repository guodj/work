
filelist = findfile('3DION*.bin')
nFiles = n_elements(filelist)

Is3DAll = 0
if (nFiles eq 1) then begin
   Is3DAll = 1
   filelist = findfile('3DALL_t000321_1[012345]*.bin')
   nFiles = n_elements(filelist)
endif
   
vd = fltarr(nFiles)
ut = fltarr(nFiles)
local = fltarr(nFiles)

ref_lon = 283.13

iVar = [0,1,2,33]

apexfile = getenv('IDL_EXTRAS')+'apex2002.dat'
print, 'reading file : ',apexfile
readapex, apexfile, apex

for iFile = 0, nFiles-1 do begin

   file = filelist(iFile)
   print, "Reading File : ", file

   gitm_read_bin_1var,file, data, time, nVars, Vars, version, VarsToGet = iVar
   c_r_to_a, itime, time

;   read_thermosphere_file, file, nvars, nalts, nlats, nlons,vars,data, $
;                           nBLKlat, nBLKlon, nBLK, iTime, Version

   ut(iFile) = float(iTime(3))+float(iTime(4))/60.0 + float(itime(5))/3600.0
   local(iFile) = (ut(iFile) + ref_lon/15.0) mod 24.0
   if (iFile gt 0 and ut(iFile) eq 0) then ut(iFile) = 24.0
   sMin  = chopr('0'+tostr(itime(4)),2)
   sHour = chopr('0'+tostr(itime(3)),2)
   sDay  = chopr('0'+tostr(itime(2)),2)

   psfile = 'nmf2_'+sDay+sHour+sMin+'.ps'

   nLons = n_elements(data(0,*,0,0))
   nLats = n_elements(data(0,0,*,0))
   nAlts = n_elements(data(0,0,0,*))

   if (not Is3DAll) then begin

      print, "Not a 3DALL file!"

      lons = reform(data(0,*,0,0)/!dtor)
      d    = abs(lons - ref_lon)
      iLon = where(d eq min(d))

      if (iFile eq 0) then begin
         alldrifts = fltarr(nFiles,nLons-4)
         alllocal  = fltarr(nFiles,nLons-4)
         allmlats  = fltarr(nFiles,nLons-4)
      endif
      alllocal(iFile,*) = (ut(iFile) + lons(2:nLons-3)/15.0) mod 24.0

      iMagLat = 22
      iE = 12
      iVd = 17
      iVe = 16
      iVu = 17
      Alts  = reform(data(2,0,0,*)/1000.0)

      nmf2 = fltarr(nLats)
      hmf2 = fltarr(nLats)
      mlat = fltarr(nLats)
      vin  = fltarr(nLats)
      viu = fltarr(nLats)

      for i=0,nLats-1 do begin
         nmf2(i) = max(data(iE,iLon,i,*))
         loc     = where(data(iE,iLon,i,*) eq max(data(iE,iLon,i,*)))
         hmf2(i) = data(2,iLon,i,loc(0))
         vin(i) = data(iVe,iLon,i,loc(0))
         viu(i) = data(iVu,iLon,i,loc(0))
         mlat(i) = data(iMagLat,iLon,i,0)
      endfor

      lmin = where(abs(mlat) eq min(abs(mlat)))

      loc = where(data(iE,iLon,lmin(0),*) eq max(data(iE,iLon,lmin(0),*)))
      vd(iFile) = data(iVd,iLon,lmin(0),loc(0))

      l = where(abs(mlat) lt 30, c)
      l = [l(0)-1,l,l(c-1)+1]
      c = c + 2

      outfile = 'nmf2_'+sDay+sHour+'.ascii'
      openw,1,outfile
      printf, 1, c
      for i=0,c-1 do begin
         ii = l(i)
         printf,1,i,mlat(ii), nmf2(ii)
      endfor
      close,1

      for i=2,nLons-3 do begin
         mlatslice = reform(data(iMagLat,i,*,0))
         dp = abs(mlatslice-1.0)
         lminp = where(dp eq min(dp))
         dm = abs(mlatslice+1.0)
         lminm = where(dm eq min(dm))

         x = mlatslice(lminp(0)) / (mlatslice(lminp(0))-mlatslice(lminm(0)))

         locp = where(data(iE,i,lminp(0),*) eq max(data(iE,i,lminp(0),*)))
         locm = where(data(iE,i,lminm(0),*) eq max(data(iE,i,lminm(0),*)))
         alldrifts(iFile, i-2) = $
             (1.0-x) * data(iVd,i,lminp(0),locp(0)) + $
             (    x) * data(iVd,i,lminm(0),locm(0))
         allmlats(iFile,i-2) = $
            (1.0-x) * mlatslice(lminp(0))+ x*mlatslice(lminm(0))
      endfor

      x = alllocal(iFile,*)
      y = alldrifts(iFile,*)
      sort_a,x,y
      alllocal(iFile,*) = x
      alldrifts(iFile,*) = y

   endif else begin

      iE = 3
;      iE = 33
      Lons  = reform(data(0,*,*,0)/!dtor)
      Lats  = reform(data(1,*,*,0)/!dtor)
      Alts  = reform(data(2,*,*,*)/1000.0)
      e     = reform(data(iE,*,*,*))

      nmf2 = fltarr(nLons,nLats)
      hmf2 = fltarr(nLons,nLats)
      viu = fltarr(nLons,nLats)

      for iLon = 0,nLons-1 do for iLat=0,nLats-1 do begin
         l = where(e(iLon,iLat,*) eq max(e(iLon,iLat,*)))
         nmf2(iLon,iLat) = e(iLon,iLat,l(0))/1.0e12
         hmf2(iLon,iLat) = Alts(iLon,iLat,l(0))
      endfor

   endelse

;   setdevice, psfile, 'p', 5

   if (not Is3DAll) then begin

      plot, mlat(l), nmf2(l), $
            pos = [0.1, 0.6, 0.9, 1.0], $
            xtitle = 'Magnetic Latitude', ytitle = 'NMF2 (/m3)', $
            xrange = [-30.0,30.0], xstyle = 1, thick = 4, $
            title = 'Local Time : '+tostr(local(iFile))

      plot, mlat(l), vin(l), $
            pos = [0.1, 0.1, 0.9, 0.5], $
            xtitle = 'Magnetic Latitude', ytitle = 'Vi (m/s)', $
            xrange = [-30.0,30.0], xstyle = 1, thick = 4, /noerase, $
            yrange = [-100,100]

      oplot, mlat(l), viu(l), linestyle = 2, thick = 4

   endif else begin

;      makect, 'mid'

      maxi = max(nmf2)
      mini = 0.0
      levels = findgen(31)/30 * (maxi-mini) + mini
      l = where(nmf2 lt levels(1),c)
      if (c gt 0) then nmf2(l) = levels(1)
      l = where(nmf2 gt levels(29),c)
      if (c gt 0) then nmf2(l) = levels(29)

      colortitle = 'NmF2 (1e12 /m3)'
      title = colortitle
      maxrange = 40.0
      thermo_threeplot, psfile, nmf2, time, $
                        lons*!dtor, lats*!dtor, mini, maxi, $
                        title, colortitle, maxrange, apex=apex
;, $
;                        vn = vn, ve = ve

;      contour, nmf2, lons, lats, levels = levels, $
;               pos = [0.1, 0.6, 0.9, 1.0], $
;               xtitle = 'Longitude (deg)', ytitle = 'Latitude (deg)', $
;               xrange = [0.0,360.0], xstyle = 1, $
;               yrange = [-90.0,90.0], ystyle = 1, $
;               /fill, xtickv = findgen(9)*45.0, xticks = 9, xminor = 3, $
;               ytickv = findgen(5)*45.0-90.0, yticks = 5, yminor = 3
;
;      its = itime
;      c_a_to_r, its, time
;      its(3:5) = 0
;      c_a_to_r, its, basetime
;      hour = (time/3600.0 mod 24.0) + fix((time-basetime)/(24.0*3600.0))*24.0
;      localtime = (Lons/15.0 + hour) mod 24.0
;    
;      angle = 23.0 * !dtor * $
;              sin((jday(its(0),its(1),its(2)) - $
;                   jday(its(0),3,21))*2*!pi/365.0)
;      sza =  acos(sin(angle)*sin(Lats*!dtor) + $
;                  cos(angle)*cos(Lats*!dtor) * $ 
;                  cos(!pi*(LocalTime-12.0)/12.0))/!dtor
;
;      contour, sza, lons, lats, levels = [45,90], thick = 4, /over, $
;               c_linestyle = 1
;
;      if (n_elements(apex) eq 0) then begin
;         apexfile = getenv('IDL_EXTRAS')+'apex2002.dat'
;         print, 'reading file : ',apexfile
;         readapex, apexfile, apex
;      endif
;      contour, apex.alats, apex.glons, apex.glats, $
;               levels = [0.0], thick = 4, $
;               xrange = [0.0,360.0], xstyle = 1, $
;               yrange = [-90.0,90.0], ystyle = 1, /over, $
;               c_linestyle = 2
;
;      ctpos = [0.91, 0.6, 0.93, 1.0]
;      title = 'NmF2 (1e12 /m3)'
;      plotct, 255, ctpos, mm(levels)/1.0e12, title, /right
;
;      ; HMF2
;
;      maxi = max(alts)
;      mini = min(alts)
;      levels = findgen(31)/30 * (maxi-mini) + mini
;      l = where(hmf2 lt levels(1),c)
;      if (c gt 0) then hmf2(l) = levels(1)
;      l = where(hmf2 gt levels(29),c)
;      if (c gt 0) then hmf2(l) = levels(29)
;
;      contour, hmf2, lons, lats, levels = levels, $
;               pos = [0.1, 0.1, 0.9, 0.5], $
;               xtitle = 'Longitude (deg)', ytitle = 'Latitude (deg)', $
;               xrange = [0.0,360.0], xstyle = 1, $
;               yrange = [-90.0,90.0], ystyle = 1, $
;               /fill, xtickv = findgen(9)*45.0, xticks = 9, xminor = 3, $
;               ytickv = findgen(5)*45.0-90.0, yticks = 5, yminor = 3, /noerase
;
;      contour, sza, lons, lats, levels = [45,90], thick = 4, /over, $
;               c_linestyle = 1
;
;      contour, apex.alats, apex.glons, apex.glats, $
;               levels = [0.0], thick = 4, $
;               xrange = [0.0,360.0], xstyle = 1, $
;               yrange = [-90.0,90.0], ystyle = 1, /over, $
;               c_linestyle = 2
;
;      ctpos = [0.91, 0.1, 0.93, 0.5]
;      title = 'HmF2 (km)'
;      plotct, 255, ctpos, mm(levels), title, /right

   endelse

   closedevice

   if (not Is3DAll) then begin

      outfile = 'hmf2_'+sDay+sHour+'.ascii'
      psfile = 'hmf2_'+sDay+sHour+'.ps'
      openw,1,outfile
      printf, 1, c
      for i=0,c-1 do begin
         ii = l(i)
         printf,1,i, mlat(ii), hmf2(ii)
      endfor
      close,1

      setdevice, psfile, 'p', 4
      plot, mlat(l), hmf2(l)/1000.0, $
            pos = [0.1, 0.3, 0.9, 0.7], $
            xtitle = 'Magnetic Latitude', ytitle = 'hmf2 (km)', $
            xrange = [-30.0,30.0], xstyle = 1, thick = 4, $
            title = 'Local Time : '+tostr(local(iFile))
      closedevice

      outfile = 'eden_'+sDay+sHour+'.ascii'
      psfile = 'eden_'+sDay+sHour+'.ps'
      openw,1,outfile
      printf, 1, nAlts-4
      for i=2,nAlts-3 do begin
         printf,1,i, data(iE,iLon,lmin,i)
      endfor
      close,1

      setdevice, psfile, 'p', 4
      plot, data(iE,iLon,lmin,*), Alts, $
            pos = [0.3, 0.2, 0.8, 0.9], $
            ytitle = 'Altitude (km)', xtitle = 'Electron Density (/m3)', $
            thick = 4, yrange = [100.0,1000.0], ystyle = 1, $
            title = 'Local Time : '+tostr(local(iFile))
      closedevice

   endif

endfor

if (not Is3DAll) then begin

   outfile = 'verticaldrift.ascii'
   psfile = 'verticaldrift.ps'
   openw,1,outfile
   printf, 1, nFiles
   for i=0,nFiles-1 do begin
      printf,1,i, vd(i)
   endfor
   close,1

   setdevice, psfile, 'p', 4
   plot, local, vd, psym = 4, $
         pos = [0.1, 0.3, 0.9, 0.7], $
         ytitle = 'Vertical Drift (m/s)', xtitle = 'Local Time (hours)', $
         thick = 4, yrange = [-30.0,30.0], ystyle = 1, $
         xrange = [0.0,24.0], xstyle =1

   oplot, [0,25],[0.0,0.0], linestyle = 1
   closedevice
   
   psfile ='alldrifts.ps'
   setdevice, psfile, 'p', 4
   
   plot, alllocal(0,*), alldrifts(0,*), $
         pos = [0.1, 0.3, 0.9, 0.7], $
         ytitle = 'Vertical Drift (m/s)', xtitle = 'Local Time (hours)', $
         thick = 4, yrange = [-50.0,50.0], ystyle = 1, $
         xrange = [0.0,24.0], xstyle =1, /nodata
   makect,'bry'

   iStart = 1
   if (nFiles ge 48) then iStart = 24

   for i=iStart,nFiles-1 do begin
      c = 250.0/24.0 * ut(i)
      print, c
      oplot, alllocal(i,*), alldrifts(i,*), thick = 4, color = c
   endfor

   oplot, [0,25],[0.0,0.0], linestyle = 1
   closedevice
   
endif


end
