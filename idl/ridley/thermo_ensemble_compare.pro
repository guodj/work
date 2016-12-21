
if (n_elements(files) eq 0) then files = ''
files = ask('files in the ensemble to plot',files)
filelist = findfile(files)

if (n_elements(basefiles) eq 0) then files = ''
basefiles = ask('basefiles in the ensemble to plot',basefiles)
basefilelist = findfile(basefiles)

;filelist = findfile('Flares0[01234]/data/3DALL_t110215_030000.bin')
;filelist = findfile('Flares0[56789]/data/3DALL_t110215_030000.bin')
;filelist = findfile('Flares1[01234]/data/3DALL_t110215_040000.bin')
;filelist = findfile('Flares1[56789]/data/3DALL_t110215_030000.bin')

gitm_read_bin, filelist, alldata, time, nVars, Vars, version
gitm_read_bin, basefilelist, basedata, basetime, nVars, Vars, version

alt = reform(alldata(0,2,*,*,*)) / 1000.0
lat = reform(alldata(0,1,*,*,*)) / !dtor
lon = reform(alldata(0,0,*,*,*)) / !dtor

nLons = n_elements(lon(*,0,0))
nLats = n_elements(lon(0,*,0))
nAlts = n_elements(lon(0,0,*))
nFiles = n_elements(time)

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

if (n_elements(smini2) eq 0) then smini2 = '0.0'
if (n_elements(smaxi2) eq 0) then smaxi2 = '0.0'
smini2 = ask('minimum percent diff (0.0 for automatic)',smini2)
smaxi2 = ask('maximum percent diff (0.0 for automatic)',smaxi2)

if (n_elements(smini3) eq 0) then smini3 = '0.0'
if (n_elements(smaxi3) eq 0) then smaxi3 = '0.0'
smini3 = ask('minimum std dev % (0.0 for automatic)',smini3)
smaxi3 = ask('maximum std dev % (0.0 for automatic)',smaxi3)

if (n_elements(smini4) eq 0) then smini4 = '0.0'
if (n_elements(smaxi4) eq 0) then smaxi4 = '0.0'
smini4 = ask('minimum diff / std dev (0.0 for automatic)',smini4)
smaxi4 = ask('maximum diff / std dev (0.0 for automatic)',smaxi4)

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

stddata      = fltarr(nLons, nLats)
meanbasedata = fltarr(nLons, nLats)
stddiffdata  = fltarr(nLons, nLats)
meandiffdata = fltarr(nLons, nLats)

for iLon = 0, nLons-1 do for iLat = 0, nLats-1 do begin

   meanbasedata(iLon,iLat) = mean(basedata(*,iVar,iLon,iLat,iAlt))
   d = alldata(*,iVar,iLon,iLat,iAlt) - basedata(*,iVar,iLon,iLat,iAlt)
   
   meandiffdata(iLon,iLat)  = mean(d)/meanbasedata(iLon,iLat)*100.0
   stddiffdata(iLon,iLat)   = stddev(d)/meanbasedata(iLon,iLat)*100.0
   stddata(iLon,iLat)       = mean(d)/stddev(basedata(*,iVar,iLon,iLat,iAlt))

endfor

Lon1D = reform(lon(*,0,iAlt))
Lat1D = reform(lat(0,*,iAlt))

c_r_to_a, itime, time(0)
c_a_to_s, itime, sDate

ut = float(itime(3)) + float(itime(4))/60.0 + float(itime(5))/3600.0
utrot = ut * 15.0

setdevice, psfile, 'l', 5

makect,'mid'

if (IsPolar) then begin

   ppp = 2
   space = 0.01
   pos_space, ppp, space, sizes

   get_position, ppp, space, sizes, 1, pos

   if (not IsNorth) then Lat1D=-Lat1D

   mini2 = float(smini2)
   maxi2 = float(smaxi2)
   if (mini2 eq 0.0 and maxi2 eq 0.0) then begin
      mini2 = min(stddata)
      maxi2 = max(stddata)
   endif
   no12=1
   MaxRange = 90.0-MinLat
   contour_circle, stddata, Lon1D+utrot, Lat1D, $
                   no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                   pos = pos, $
                   maxrange = MaxRange, $
                   colorbar = vars(iVar), $
                   mini = mini2, maxi = maxi2

   get_position, ppp, space, sizes, 0, pos

   mini = float(smini)
   maxi = float(smaxi)
   if (mini eq 0.0 and maxi eq 0.0) then begin
      mini = min(meanbasedata)
      maxi = max(meanbasedata)
   endif
   no12=0
   no00=1
   contour_circle, meanbasedata, Lon1D+utrot, Lat1D, $
                   no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                   pos = pos, $
                   maxrange = MaxRange, $
                   colorbar = vars(iVar), $
                   mini = mini, maxi = maxi

endif else begin

   ppp = 6
   space = 0.01
   pos_space, ppp, space, sizes, nx = 2

   lon2d = reform(lon(*,*,iAlt))
   lat2d = reform(lat(*,*,iAlt))

   loc = where(abs(lat2d) le 90.0 and $
               lon2d ge 0.0 and $
               lon2d lt 360.0, count)

   lon2d = (lon2d(loc) + utrot) mod 360.0
   lat2d = lat2d(loc)

   nLevels = 31

   get_position, ppp, space, sizes, 0, pos
   dx = pos(2)-pos(0)
   pos(0) = 0.0
   pos(2) = pos(0)+dx*2

   mini = float(smini)
   maxi = float(smaxi)
   if (mini eq 0.0 and maxi eq 0.0) then begin
      mini = min(meanbasedata)
      maxi = max(meanbasedata)
   endif

   contour, meanbasedata(loc), lon2d, lat2d, $
            /noerase, pos = pos, /fill, $
            nlevels = 31, $
            xstyle = 5, ystyle = 5, /irr, $
            levels = findgen(nlevels)*(maxi-mini)/(nlevels-1)+mini, $
            c_colors = findgen(nlevels)*250/(nlevels-1) + 3, $
            yrange = [-90,90], xrange = [0,360]

   !p.position = pos
   t = (180.0-utrot+360.0) mod 360.0
   map_set, 0.0, t, /noerase
   map_continents, color = 0
   !p.position = -1

   ctpos = [pos(0), pos(1)-0.03, pos(2), pos(1)-0.01]
   minmax = [mini,maxi]
   title = Vars(iVar)
   plotct, 255, ctpos, minmax, title, /bottom

   xyouts, pos(0), pos(3)+0.01, sDate, /norm
   xyouts, pos(2), pos(3)+0.01, tostr(alt(0,0,iAlt))+' km', /norm, alignment=1.0

   get_position, ppp, space, sizes, 1, pos
   dx = pos(2)-pos(0)
   pos(0) = 0.5
   pos(2) = pos(0)+dx*2

   mini2 = float(smini2)
   maxi2 = float(smaxi2)
   if (mini2 eq 0.0 and maxi2 eq 0.0) then begin
      mini2 = min(meandiffdata)
      maxi2 = max(meandiffdata)
   endif

   contour, meandiffdata(loc), lon2d, lat2d, $
            /noerase, pos = pos, /fill, $
            nlevels = 31, $
            xstyle = 5, ystyle = 5, /irr, $
            levels = findgen(nlevels)*(maxi2-mini2)/(nlevels-1)+mini2, $
            c_colors = findgen(nlevels)*250/(nlevels-1) + 3, $
            yrange = [-90,90], xrange = [0,360]

   !p.position = pos
   map_set, 0.0, t, /noerase
   map_continents, color = 0
   !p.position = -1

   ctpos = [pos(0), pos(1)-0.03, pos(2), pos(1)-0.01]
   minmax = [mini2,maxi2]
   title = 'Mean Difference ('+Vars(iVar)+') (%)'
   plotct, 255, ctpos, minmax, title, /bottom

   get_position, ppp, space, sizes, 4, pos
   dx = pos(2)-pos(0)
   pos(0) = 0.0
   pos(2) = pos(0)+dx*2
   pos([1,3])  = pos([1,3]) + 0.25

   mini3 = float(smini3)
   maxi3 = float(smaxi3)
   if (mini3 eq 0.0 and maxi3 eq 0.0) then begin
      mini3 = min(stddiffdata)
      maxi3 = max(stddiffdata)
   endif

   contour, stddiffdata(loc), lon2d, lat2d, $
            /noerase, pos = pos, /fill, $
            nlevels = 31, $
            xstyle = 5, ystyle = 5, /irr, $
            levels = findgen(nlevels)*(maxi3-mini3)/(nlevels-1)+mini3, $
            c_colors = findgen(nlevels)*250/(nlevels-1) + 3, $
            yrange = [-90,90], xrange = [0,360]

   !p.position = pos
   t = (180.0-utrot+360.0) mod 360.0
   map_set, 0.0, t, /noerase
   map_continents, color = 0
   !p.position = -1

   ctpos = [pos(0), pos(1)-0.03, pos(2), pos(1)-0.01]
   minmax = [mini3,maxi3]
   title = 'stddev of diff ('+Vars(iVar)+') (%)'
   plotct, 255, ctpos, minmax, title, /bottom

   get_position, ppp, space, sizes, 5, pos
   dx = pos(2)-pos(0)
   pos(0) = 0.5
   pos(2) = pos(0)+dx*2
   pos([1,3])  = pos([1,3]) + 0.25

   mini4 = float(smini4)
   maxi4 = float(smaxi4)
   if (mini4 eq 0.0 and maxi4 eq 0.0) then begin
      mini4 = min(stddata)
      maxi4 = max(stddata)
   endif

   contour, stddata(loc), lon2d, lat2d, $
            /noerase, pos = pos, /fill, $
            nlevels = 31, $
            xstyle = 5, ystyle = 5, /irr, $
            levels = findgen(nlevels)*(maxi4-mini4)/(nlevels-1)+mini4, $
            c_colors = findgen(nlevels)*250/(nlevels-1) + 3, $
            yrange = [-90,90], xrange = [0,360]

   !p.position = pos
   t = (180.0-utrot+360.0) mod 360.0
   map_set, 0.0, t, /noerase
   map_continents, color = 0
   !p.position = -1

   ctpos = [pos(0), pos(1)-0.03, pos(2), pos(1)-0.01]
   minmax = [mini4,maxi4]
   title = 'diff / std(base) ('+Vars(iVar)+') (ratio)'
   plotct, 255, ctpos, minmax, title, /bottom

endelse

closedevice

end

