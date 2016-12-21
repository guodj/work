
if (n_elements(filelist) eq 0) then begin
   filelist = findfile("-t /*.save")
   if (strlen(filelist(0)) eq 0) then filelist = findfile("-t *.bin")
endif

filelist = ask('filename to plot',filelist(0))
filelist = findfile(filelist)

nfiles = n_elements(filelist)

if (nFiles eq 1 and strlen(filelist(0)) eq 0) then begin
   print, 'can not find file!'
   stop
endif

psfile = filelist(0)+'.ps'
psfile = ask('ps file name',psfile)

gitm_read_header, filelist(0), GitmTimes, nVars, Vars, $
                  nLons, nLats, nAlts, version

display,vars

if (n_elements(iVar) eq 0) then iVar = '15' else iVar = tostr(iVar)
iVar = fix(ask('which var to plot',iVar))
variable = vars(iVar)

if (n_elements(iVar2) eq 0) then iVar2 = '-1' else iVar2 = tostr(iVar2)
iVar2 = fix(ask('which var to plot with line contours (-1 for none)',iVar2))
if (iVar2 ge 0 and iVar2 ne iVar) then variable2 = vars(iVar2)

if (n_elements(plotlog) eq 0) then plotlog='n' $
else if (plotlog) then plotlog='y' else plotlog='n' 
plotlog = ask('whether you want log or not (y/n)',plotlog)
if (strpos(plotlog,'y') eq 0) then plotlog = 1 else plotlog = 0

if plotlog then variable = 'log('+variable+')'

if (n_elements(IsZonalAverage) eq 0) then IsZonalAverage='0' $
else IsZonalAverage = tostr(IsZonalAverage)
IsZonalAverage = $
   fix(ask('whether you want a Zonal Average (0=n,1=y)',$
           IsZonalAverage))
if (not IsZonalAverage) then begin

   VarsToGet = [0]
   gitm_read_bin_1var, filelist(0), data, time1, nVars, Vars, version, $
                       VarsToGet = VarsToGet
   longitudes = reform(data(0,*,0,0)/!dtor)
   nLons = n_elements(longitudes)

   if (n_elements(IsLocalTime) eq 0) then IsLocalTime = 0
   IsLocalTime = fix(ask('whether you want local time (1) or longitude (0)',$
                         tostr(IsLocalTime)))
   if (IsLocalTime) then begin
      if (n_elements(LocalTime) eq 0) then LocalTime = 12.0
      LocalTime = float(ask('local time to plot',string(LocalTime)))
   endif else begin
      for i=0,nlons-1 do print, tostr(i)+'. '+string(longitudes(i))
      iLon = fix(ask('which longitude to plot','0'))
   endelse

endif else iLon = 0

mini = float(ask('minimum for color (0.0 for automatic)','0.0'))
maxi = float(ask('maximum for color (0.0 for automatic)','0.0'))
if (mini eq maxi) then iAuto = 1 else iAuto = 0

iAuto2 = 1
if (iVar2 ge 0 and iVar2 ne iVar and iAuto eq 0) then begin
   mini2 = float(ask('minimum for lines (0.0 for automatic)','0.0'))
   maxi2 = float(ask('maximum for lines (0.0 for automatic)','0.0'))
   if (mini2 eq maxi2) then iAuto2 = 1 else iAuto2 = 0
endif else begin
   if (iAuto) then iAuto2 = 1 else begin
      if (iVar2 eq iVar) then begin 
         iAuto2 = 0
         mini2 = mini
         maxi2 = maxi
      endif
   endelse
endelse

if (n_elements(minlat) eq 0) then minlat = -90.0
if (n_elements(maxlat) eq 0) then maxlat =  90.0
minlat = float(ask('minimum latitude to plot',string(minlat)))
maxlat = float(ask('maximum latitude to plot',string(maxlat)))

VarsToGet = [1,2,iVar]
if (iVar2 ge 0 and iVar2 ne iVar) then VarsToGet = [VarsToGet,iVar2]

for iFile = 0, nFiles-1 do begin

   file = filelist(iFile)
   print, file
   gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                       VarsToGet = VarsToGet
   c_r_to_a, itime, time
   c_a_to_s, itime, stime
   ut = float(itime(3)) + itime(4)/60.0 + itime(5)/3600.0

   if (iFile eq 0) then begin
      lats = reform(data(0,0,*,*))/!dtor
      alts = reform(data(1,0,*,*))/1000.0
      nLats = n_elements(lats(*,0))
      nAlts = n_elements(alts(0,*))
      minalt = alts(0,2)
      maxalt = alts(0,nAlts-3)
   endif

   iV = n_elements(VarsToGet)-1
   iV2 = -1
   if (iVar2 ge 0) then begin
      iV2 = iV
      if (iVar2 ne iVar) then begin
         iV--
;         iTemp = iV
;         iV = iV2
;         iV2 = iTemp
      endif
   endif

   print, iV, iV2, iVar, iVar2
   print, VarsToGet
;   stop

   if (IsZonalAverage) then begin
      datatoplot  = fltarr(nLats,nAlts)
      datatoplot2 = fltarr(nLats,nAlts)
      for iLat = 0,nLats-1 do for iAlt = 0,nAlts-1 do begin
         datatoplot(iLat,iAlt) = mean(data(iV,2:nLons-3,iLat,iAlt))
         if (iV2 ge 0) then $
            datatoplot2(iLat,iAlt) = mean(data(iV2,2:nLons-3,iLat,iAlt))
      endfor
      location = 'Zonal Mean'
   endif else begin
      if (IsLocalTime) then begin
         lts = (ut + longitudes/15.0) mod 24.0
         diff = abs(lts - LocalTime)
         iLon = where(diff eq min(diff))
         iLon = iLon(0)
      endif
      datatoplot=reform(data(iV,iLon,*,*))
      if (iV2 ge 0) then $
         datatoplot2=reform(data(iV2,iLon,*,*))
      longitude = longitudes(iLon)
      location = string(longitude,format='(f5.1)')+' deg Longitude'
   endelse

   if (nFiles gt 1) then begin
      p = strpos(psfile,'.ps')
      if (p gt -1) then psfile = strmid(psfile,0,p)
      psfile_final = psfile+'_'+chopr('000'+tostr(iFile),4)+'.ps'
   endif else begin
      psfile_final = psfile
   endelse

   setdevice, psfile_final, 'p', 5

   if (plotlog) then begin
      datatoplot = alog10(datatoplot)
      if (iVar2 eq iVar) then datatoplot2 = alog10(datatoplot2)
   endif

   yMin = 0.3
   yMax = 0.8

   if (iAuto) then begin
      determine_min_max, datatoplot, mini, maxi
      if (maxi gt 0 and mini lt 0) then begin
         maxi = max(abs([mini,maxi]))*1.02
         mini = -maxi
      endif else begin
         if (maxi gt 0) then maxi = maxi*1.02 else maxi = maxi*0.98
         if (mini gt 0) then mini = mini*0.98 else mini = mini*1.02
      endelse
   endif
   levels = findgen(61)*(maxi-mini)/60 + mini

   if (iAuto2 and iVar2 ge 0) then begin
      determine_min_max, datatoplot2, mini2, maxi2
      if (maxi2 gt 0 and mini2 lt 0) then begin
         maxi2 = max(abs([mini2,maxi2]))*1.02
         mini2 = -maxi2
      endif else begin
         if (maxi2 gt 0) then maxi2 = maxi2*1.02 else maxi2 = maxi2*0.98
         if (mini2 gt 0) then mini2 = mini2*0.98 else mini2 = mini2*1.02
      endelse
   endif
   if (iV2 ge 0) then levels2 = findgen(11)*(maxi2-mini2)/10 + mini2

   if (maxi gt 0 and mini lt 0) then makect, 'mid' else makect,'all'

   l = where(datatoplot gt levels(59),c)
   if (c gt 0) then datatoplot(l) = levels(59)
   l = where(datatoplot lt levels(1),c)
   if (c gt 0) then datatoplot(l) = levels(1)

   contour, datatoplot, lats, alts, levels = levels, /fill, $
            xstyle = 1, ystyle = 1, $
            xtitle = 'Latitude (deg)', ytitle = 'Altitude (km)', $
            title = variable+' at '+location, $
            xrange = [minlat,maxlat], $
            yrange = [minalt,maxalt], $
            pos = [0.05, yMin, 0.93, yMax]

   if (iV2 ge 0) then begin
      contour, datatoplot2, lats, alts, levels = levels2, $
               xstyle = 1, ystyle = 1, $
               xrange = [minlat,maxlat], $
               yrange = [minalt,maxalt], $
               pos = [0.05, yMin, 0.93, yMax], /noerase, /follow, $
               c_linestyle = (levels2 lt 0.0)
   endif

   xyouts, 0.05, yMax+0.005, strmid(stime,0,9), /norm
   xyouts, 0.93, yMax+0.005, strmid(stime,10,5)+' UT', /norm, align=1.0

   title = variable
   pos = [0.94, yMin, 0.96, yMax]
   plotct,254,pos,mm(levels),variable,/right

   VarsToGet = [iVar]
   if (iVar2 ge 0 and iVar2 ne iVar) then VarsToGet = [VarsToGet,iVar2]

   closedevice

endfor

end

