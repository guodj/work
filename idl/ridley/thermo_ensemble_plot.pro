
if (n_elements(files) eq 0) then files = ''
files = ask('files in the ensemble to plot',files)

;filelist = findfile('Flares0[01234]/data/3DALL_t110215_030000.bin')
;filelist = findfile('Flares0[56789]/data/3DALL_t110215_030000.bin')
;filelist = findfile('Flares1[01234]/data/3DALL_t110215_040000.bin')
;filelist = findfile('Flares1[56789]/data/3DALL_t110215_030000.bin')
filelist = findfile(files)

gitm_read_bin, filelist, alldata, time, nVars, Vars, version

alt = reform(alldata(0,2,*,*,*)) / 1000.0
lat = reform(alldata(0,1,*,*,*)) / !dtor
lon = reform(alldata(0,0,*,*,*)) / !dtor

nLons = n_elements(lon(*,0,0))
nLats = n_elements(lon(0,*,0))
nAlts = n_elements(lon(0,0,*))

for i=0,nvars-1 do print, tostr(i)+'. '+vars(i)
if (n_elements(iVar) eq 0) then iVar = '3' else iVar = tostr(iVar)
iVar = fix(ask('which var to plot',iVar))

if (n_elements(psfile) eq 0) then psfile = 'tmp.ps'
psfile = ask('psfile name',psfile)

for i=0,nalts-1 do print, tostr(i)+'. '+string(alt(2,2,i))
if (n_elements(iAlt) eq 0) then iAlt='0' else iAlt=tostr(iAlt)
iAlt = fix(ask('which altitude to plot',iAlt))

if (n_elements(IsPolar) eq 0) then IsPolar='1' else IsPolar = tostr(IsPolar)
IsPolar = fix(ask('polar (1) or non-polar (0)',IsPolar))

if (IsPolar) then begin
   if (n_elements(IsNorth) eq 0) then IsNorth='1'
   IsNorth = fix(ask('North (1) or South (0)',IsNorth))
   MinLat  = abs(float(ask('minimum latitude to plot','50.0')))
endif

if (n_elements(smini) eq 0) then smini = '0.0'
if (n_elements(smaxi) eq 0) then smaxi = '0.0'
smini = ask('minimum (0.0 for automatic)',smini)
smaxi = ask('maximum (0.0 for automatic)',smaxi)

if (n_elements(sminis) eq 0) then sminis = '0.0'
if (n_elements(smaxis) eq 0) then smaxis = '0.0'
sminis = ask('minimum (0.0 for automatic)',sminis)
smaxis = ask('maximum (0.0 for automatic)',smaxis)

if (n_elements(plotVector) eq 0) then plotvector='y' $
else if (plotvector) then plotvector='y' else plotvector='n'
plotvector=ask('whether you want vectors or not (y/n)',plotvector)
if strpos(plotvector,'y') eq 0 then plotvector=1 $
else plotvector = 0

if (plotvector) then begin

   PlotNeutrals = fix(ask('plot neutral winds (1) or ions (0)','1'))

   print,'-1  : automatic selection'
   factors = [1.0, 5.0, 10.0, 20.0, 25.0, $
              50.0, 75.0, 100.0, 150.0, 200.0,300.0]
   nfacs = n_elements(factors)
   for i=0,nfacs-1 do print, tostr(i+1)+'. '+string(factors(i)*10.0)
   vector_factor = fix(ask('velocity factor','-1'))
endif else vector_factor = 0

meandata = fltarr(nLons, nLats)
stddata  = fltarr(nLons, nLats)

for iLon = 0, nLons-1 do for iLat = 0, nLats-1 do begin

   meandata(iLon,iLat) = mean(alldata(*,iVar,iLon,iLat,iAlt))
   stddata(iLon,iLat)  = stddev(alldata(*,iVar,iLon,iLat,iAlt))/meandata(iLon,iLat)*100.0

endfor

Lon1D = reform(lon(*,0,iAlt))
Lat1D = reform(lat(0,*,iAlt))

c_r_to_a, itime, time(0)
c_a_to_s, itime, sDate

ut = float(itime(3)) + float(itime(4))/60.0 + float(itime(5))/3600.0
utrot = ut * 15.0

setdevice, psfile, 'p', 5

makect,'mid'

mini = float(smini)
maxi = float(smaxi)
if (mini eq 0.0 and maxi eq 0.0) then begin
   mini = min(stddata)
   maxi = max(stddata)
endif

if (IsPolar) then begin

   ppp = 2
   space = 0.01
   pos_space, ppp, space, sizes

   get_position, ppp, space, sizes, 1, pos

   if (not IsNorth) then Lat1D=-Lat1D

   minis = float(sminis)
   maxis = float(smaxis)
   if (minis eq 0.0 and maxis eq 0.0) then begin
      minis = min(stddata)
      maxis = max(stddata)
   endif
   no12=1
   MaxRange = 90.0-MinLat
   contour_circle, stddata, Lon1D+utrot, Lat1D, $
                   no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                   pos = pos, $
                   maxrange = MaxRange, $
                   colorbar = vars(iVar), $
                   mini = minis, maxi = maxis

   get_position, ppp, space, sizes, 0, pos

   mini = float(smini)
   maxi = float(smaxi)
   if (mini eq 0.0 and maxi eq 0.0) then begin
      mini = min(meandata)
      maxi = max(meandata)
   endif
   no12=0
   no00=1
   contour_circle, meandata, Lon1D+utrot, Lat1D, $
                   no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                   pos = pos, $
                   maxrange = MaxRange, $
                   colorbar = vars(iVar), $
                   mini = mini, maxi = maxi

endif else begin

   ppp = 3
   space = 0.01
   pos_space, ppp, space, sizes, nx = 1

   get_position, ppp, space, sizes, 1, pos
   dx = pos(2)-pos(0)
   pos(0) = pos(0)-dx/2
   pos(2) = pos(2)+dx/2

   lon2d = reform(lon(*,*,iAlt))
   lat2d = reform(lat(*,*,iAlt))

   loc = where(abs(lat2d) le 90.0 and $
               lon2d ge 0.0 and $
               lon2d lt 360.0, count)

   lon2d = (lon2d(loc) + utrot) mod 360.0
   lat2d = lat2d(loc)

   nLevels = 31

   minis = float(sminis)
   maxis = float(smaxis)
   if (minis eq 0.0 and maxis eq 0.0) then begin
      minis = min(stddata)
      maxis = max(stddata)
   endif

   contour, stddata(loc), lon2d, lat2d, $
            /noerase, pos = pos, /fill, $
            nlevels = 31, $
            xstyle = 5, ystyle = 5, /irr, $
            levels = findgen(nlevels)*(maxis-minis)/(nlevels-1)+minis, $
            c_colors = findgen(nlevels)*250/(nlevels-1) + 3, $
            yrange = [-90,90], xrange = [0,360]

   !p.position = pos
   t = (180.0-utrot+360.0) mod 360.0
   map_set, 0.0, t, /noerase
   map_continents, color = 0
   !p.position = -1

   ctpos = [pos(2)+0.01, pos(1), pos(2)+0.03, pos(3)]
   minmax = [minis,maxis]
   title = 'stddev('+Vars(iVar)+') (%)'
   plotct, 255, ctpos, minmax, title, /right

   get_position, ppp, space, sizes, 0, pos
   dx = pos(2)-pos(0)
   pos(0) = pos(0)-dx/2
   pos(2) = pos(2)+dx/2

   mini = float(smini)
   maxi = float(smaxi)
   if (mini eq 0.0 and maxi eq 0.0) then begin
      mini = min(meandata)
      maxi = max(meandata)
   endif

   contour, meandata(loc), lon2d, lat2d, $
            /noerase, pos = pos, /fill, $
            nlevels = 31, $
            xstyle = 5, ystyle = 5, /irr, $
            levels = findgen(nlevels)*(maxi-mini)/(nlevels-1)+mini, $
            c_colors = findgen(nlevels)*250/(nlevels-1) + 3, $
            yrange = [-90,90], xrange = [0,360]

   !p.position = pos
   map_set, 0.0, t, /noerase
   map_continents, color = 0
   !p.position = -1

   ctpos = [pos(2)+0.01, pos(1), pos(2)+0.03, pos(3)]
   minmax = [mini,maxi]
   title = Vars(iVar)
   plotct, 255, ctpos, minmax, title, /right

   xyouts, pos(0), pos(3)+0.01, sDate, /norm
   xyouts, pos(2), pos(3)+0.01, tostr(alt(0,0,iAlt))+' km', /norm, alignment=1.0

endelse

closedevice

end

