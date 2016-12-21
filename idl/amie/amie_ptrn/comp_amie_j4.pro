
;--------------------------------------------------------------
; Get Inputs from the user
;--------------------------------------------------------------

if (n_elements(amie_file) eq 0) then begin

  initial_guess = findfile('-t b*')
  initial_guess = initial_guess(0)
  if strlen(initial_guess) eq 0 then initial_guess='b970101'

endif else initial_guess = amie_file

amie_file = ask('AMIE binary file name',initial_guess)
psfile = ask('ps file',amie_file+'_j4.ps')

mr = 40.0

read_amie_binary, amie_file, AMIEData, lats, mlts, AMIETime, fields, 	$
                  imf, ae, dst, hp, cpcp

nmlts = n_elements(mlts)
nlats = n_elements(lats)

; Figure out which hemisphere we are in

Hem = 1
if lats(0) lt 0 then Hem = -1

; For most current AMIE runs, Auroral Energy Flux is 6, while
; Mean energy is 5.

;ef_ = 14
;me_ = 13
ef_ = 6
me_ = 5
ha_ = 4

dmlt = mlts(1) - mlts(0)
dlat = lats(0) - lats(1)

; Determine x and y coordinates of AMIE plots

lat2d = fltarr(nmlts,nlats)
lon2d = fltarr(nmlts,nlats)

for i=0,nlats-1 do lon2d(*,i) = mlts*!pi/12.0 - !pi/2.0
for i=0,nmlts-1 do lat2d(i,*) = lats

x = (90.0-lat2d)*cos(lon2d)
y = (90.0-lat2d)*sin(lon2d)

LatLoc = where(90.0 - lats le mr)

x = x(*,LatLoc)
y = y(*,LatLoc)

c_r_to_a, itime, AMIETime(0)
J4Date = tostr(itime(0))+ $
         chopr('0'+tostr(itime(1)),2) + $
         chopr('0'+tostr(itime(2)),2)
;;chopr('0'+tostr(itime(0)),2) + $
yymm   = chopr('0'+tostr(itime(0)),2) + $
         chopr('0'+tostr(itime(1)),2)

filelist = findfile('/swrdata/raidzone1-5/dmsp/'+yymm+'/j4/19'+J4Date+'/*.ssj4')
nfiles = n_elements(filelist)

if (nfiles eq 1) then begin
  filelist = findfile('/swrdata/raidzone1-5/AMIE/'+yymm+'/Data/dmsp/J4/19'+J4Date+'/*.ssj4')
  nfiles = n_elements(filelist)
endif

setdevice, psfile, 'p', 4, 0.95

; Metrics is an array which is going to compare some values for all of
; the passes.  The 13 is for 13 metrics. The 2 is for 2 cuts through the
; aurora.
nMetrics = 13
Metrics = fltarr(nMetrics,2,nFiles)
MetricsTimes = dblarr(nFiles)

DummyData = fltarr(6)

pnBase = 8

for iFile=0,nfiles-1 do begin

  j4file = filelist(iFile)
  print, 'J4File:'  + J4file

  read_j4_data, j4file, J4Data, J4Time, GeoPos, MagPos, nPts

  ; I am tired of dealing with skimmer passes.  Lets only
  ; deal with passes which cut through the main oval
  ; in two places.  These passes are typically the ones
  ; which go about 80 degree latitude.

  loc = where(MagPos(0,*)*Hem gt 80.0, count)
  if (count eq 0) then nPts = 1

  if nPts eq 0 then nPts = 1

  ; for J4 data, Auroral Energy Flux is 1, while
  ; Mean energy is 2.

  CompIter  = intarr(nPts)
  InterLong = fltarr(nPts)
  InterLat  = fltarr(nPts)

  AMIEEFlux = fltarr(nPts)
  AMIEAveE  = fltarr(nPts)
  AMIEHall   = fltarr(nPts)

  J4EFlux = fltarr(nPts)
  J4AAveE  = fltarr(nPts)
  J4Hall   = fltarr(nPts)

  nAve = 10

  for i=1,nPts-1 do $
    if (J4Data(1,i) gt J4Data(1,i-1)*5.0 and J4Data(1,i) gt 8.0) then $
      J4Data(1,i) = J4Data(1,i-1)

  for i=0,nPts-1 do begin

    ; Find Closest Time in AMIE Data

    d = abs(AMIETime - J4Time(i))
    loc = where(d eq min(d))
    CompIter(i) = loc(0)

    ; Find Interpolation points in AMIE grid

    InterLong(i) = MagPos(2,i) / dmlt
    InterLat(i)  = (90-MagPos(0,i)) / dlat

    AMIEEFlux(i) = interpolate(AMIEData(CompIter(i),ef_,*,*), $
			InterLong(i), InterLat(i))
    AMIEAveE(i)  = interpolate(AMIEData(CompIter(i),me_,*,*), $
			InterLong(i), InterLat(i))
    AMIEHall(i)  = interpolate(AMIEData(CompIter(i),ha_,*,*), $
			InterLong(i), InterLat(i))

    ped = 40.0 * j4data(2,i) / (16.0 + j4data(2,i)^2) * sqrt(J4Data(1,i))
    j4hall(i) = 0.45 * (j4data(2,i)^0.85) * ped

    ; Lets calculate running averages to compare to
    if (i gt nAve/2) and (i lt nPts-nAve/2-1) then begin
      for j=i-nAve/2,i+nAve/2-1 do begin
        J4EFlux(i) = J4EFlux(i) + J4Data(1,j)/nAve
        J4AAveE(i)  = J4AAveE(i) + J4Data(2,j)/nAve
      endfor
    endif

  endfor

  loc = where(MagPos(0,*)*Hem gt 90-mr, nTruePts)

  if (nTruePts gt 0) then begin
    loczero = where(J4Data(1,loc) eq 0.0)
    nBadPts = n_elements(loczero)
    if (float(nBadPts)/float(nTruePts) gt 0.5) then begin
      nTruePts = 0
    endif
  endif

  if (nTruePts gt 0 and npts gt 1) then begin

    stime = min(J4Time(loc))
    etime = max(J4Time(loc))

    time_axis, stime, etime, s_time_range, e_time_range,        	$
	xtickname, xtitle, xtickvalue, xminor, xtickn

    if pnBase eq 8 then pnBase = 0 else pnBase = 8
    if pnBase eq 0 then begin
      plotdumb
      ppp = 16
      space = 0.03
      pos_space, ppp, space, sizes, ny = 4
      fac = 1.0
    endif else fac = -1.0

    pn = pnBase + 0
    get_position, ppp, space, sizes, pn, pos
    pos([1,3]) = pos([1,3]) + fac*space

    iter = mean(CompIter(loc))

    mini = 0.0
    maxi = max(AMIEData(iter,ef_,*,*))
    range = maxi
    dc      = 10.0^fix(alog10(range/20.0))
    factor  = 0.0
    while (range gt dc*20.0*factor) do factor=factor+0.1
    dc = factor*dc
    levels = findgen(21)*dc

    readct, ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
    ;readct, ncolors, getenv("IDL_EXTRAS")+"white_red.ct"
    clevels = (ncolors-1)*findgen(21)/20.0 + 1

    contour, reform(AMIEData(iter,ef_,*,LatLoc)),x,y,/follow, 		$
	xstyle = 5, ystyle = 5,						$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase, /cell_fill, c_color = clevels

  ;  xyouts, pos(0)-0.01,pos(3)-0.01, amie_time_1(n),charsize=0.9, /norm

    J4x = (90.0 - MagPos(0,loc)) * cos(MagPos(2,loc)*!pi/12.0 - !pi/2.0)
    J4y = (90.0 - MagPos(0,loc)) * sin(MagPos(2,loc)*!pi/12.0 - !pi/2.0)

    plotmlt, mr, /no06, /no00
    oplot, J4x, J4y
    xyouts, J4x(0), J4y(0), 'S'

    pn = pnBase + 1
    get_position, ppp, space, sizes, pn, pos1

    pn = pnBase + 3
    get_position, ppp, space, sizes, pn, pos2

    pos(0) = pos1(0)
    pos(1) = pos1(1)
    pos(2) = pos2(2)
    pos(3) = pos1(3)
    pos([1,3]) = pos([1,3]) + fac*space

    plot, J4Time(loc)-stime,J4Data(1,loc), pos = pos,	 		$
	xstyle = 1, /noerase,		$
	xtickname = strarr(10)+' ', xtickv=xtickvalue, 			$
	xticks = xtickn, xminor = xminor, xtitle = ' ',			$
	xrange = [s_time_range, e_time_range], $
	ytickname = strarr(20)+' ', yrange=[0.0,min([10.0,max(J4Data(1,loc))])]
    axis, yax=1, charsize=0.9, yrange=[0.0,min([10.0,max(J4Data(1,loc))])]

    xyouts, pos(0)-0.005, (pos(1)+pos(3))/2.0, 'Elec. E. (ergs/cm2/s)', $
  	/norm, alignment = 0.5, orient = 90, charsize = 0.75

    oplot, J4Time(loc)-stime, AMIEEFlux(loc), linestyle = 2
    oplot, J4Time(loc)-stime, J4EFlux(loc), linestyle = 1

    ; Let's try to figure out the metric which we want to use...

    ; First we want to determine whether we had a skimmer pass or a real
    ; pass.  To do this, we figure out whether there are two maxima
    ; in the Electron Fluxes:

    MaxLeft  = max(J4EFlux(loc(0:nTruePts/2-1)))
    LocLeft  = where(J4EFlux(loc(0:nTruePts/2-1)) eq MaxLeft)
    ; We want the right most point of this
    LocLeft  = loc(LocLeft(n_elements(LocLeft)-1))

    MaxRight  = max(J4EFlux(loc(nTruePts/2:nTruePts-1)))
    LocRight  = where(J4EFlux(loc(nTruePts/2:nTruePts-1)) eq MaxRight)
    ; We want the left most point of this
    LocRight  = loc(LocRight(0)+nTruePts/2)

    ; Left first

    TrueMax = max(J4Data(1,LocLeft-nAve/2:LocLeft+nAve/2-1))
    LocTrueMax = where(J4Data(1,LocLeft-nAve/2:LocLeft+nAve/2-1) eq TrueMax)
    LocLeft = LocTrueMax(0) + LocLeft-nAve/2

    TrueMax_Hall = max(J4Hall(LocLeft-nAve/2:LocLeft+nAve/2-1))
    LocTrueMax_Hall = where(J4Hall(LocLeft-nAve/2:LocLeft+nAve/2-1) eq TrueMax_Hall)
    LocLeft_Hall = LocTrueMax_Hall(0) + LocLeft-nAve/2

    ; Right next

    TrueMax = max(J4Data(1,LocRight-nAve/2:LocRight+nAve/2-1))
    LocTrueMax = where(J4Data(1,LocRight-nAve/2:LocRight+nAve/2-1) eq TrueMax)
    LocRight = LocTrueMax(0) + LocRight-nAve/2

    TrueMax_Hall = max(J4Hall(LocRight-nAve/2:LocRight+nAve/2-1))
    LocTrueMax_Hall = where(J4Hall(LocRight-nAve/2:LocRight+nAve/2-1) eq TrueMax_Hall)
    LocRight_Hall = LocTrueMax_Hall(0) + LocRight-nAve/2

    nPeaks = 2

    ; Test to see how close they are to each other:
    if (LocRight - LocLeft lt nAve) then nPeaks = 1

    ; if we have a single peak, then move everything to LocLeft

    if (nPeaks eq 1) and (MaxRight gt MaxLeft) then begin
      MaxLeft = MaxRight
      LocLeft = LocRight
    endif

    ; Now we have to find the AMIE peaks:

    ; Left first

    AMIEMax = max(AMIEEFlux(LocLeft-nAve:LocLeft+nAve-1))
    LocAMIEMax = where(AMIEEFlux(LocLeft-nAve:LocLeft+nAve-1) eq AMIEMax)
    AMIELocLeft = LocAMIEMax(0) + LocLeft-nAve

    AMIEMax_Hall = max(AMIEHall(LocLeft_Hall-nAve:LocLeft_Hall+nAve-1))
    LocAMIEMax_Hall = where(AMIEHall(LocLeft_Hall-nAve:LocLeft_Hall+nAve-1) eq AMIEMax_Hall)
    AMIELocLeft_Hall = LocAMIEMax_Hall(0) + LocLeft_Hall-nAve

    ; Right next

    if (nPeaks gt 1) then begin

      AMIEMax = max(AMIEEFlux(LocRight-nAve:LocRight+nAve-1))
      LocAMIEMax = where(AMIEEFlux(LocRight-nAve:LocRight+nAve-1) eq AMIEMax)
      AMIELocRight = LocAMIEMax(0) + LocRight-nAve

      AMIEMax_Hall = max(AMIEHall(LocRight_Hall-nAve:LocRight_Hall+nAve-1))
      LocAMIEMax_Hall = where(AMIEHall(LocRight_Hall-nAve:LocRight_Hall+nAve-1) eq AMIEMax_Hall)
      AMIELocRight_Hall = LocAMIEMax_Hall(0) + LocRight_Hall-nAve

    endif

    ; Now lets do some interesting things...

    metrics(0,0,iFile) = 100*(J4Data(1,LocLeft) - AMIEEFlux(AMIELocLeft))/$
	max([J4Data(1,LocLeft),AMIEEFlux(AMIELocLeft)])
    metrics(1,0,iFile) = MagPos(0,LocLeft) - MagPos(0,AMIELocLeft)

    metrics(8,0,iFile) = 100*(J4Hall(LocLeft_Hall) - AMIEHall(AMIELocLeft_Hall))/$
	max([J4Hall(LocLeft_Hall),AMIEHall(AMIELocLeft_Hall)])
    metrics(9,0,iFile) = MagPos(0,LocLeft) - MagPos(0,AMIELocLeft_Hall)

    lag = indgen(nAve+1) - nAve/2

    if (npeaks eq 1) then begin

      c = c_correlate(J4Data(1,loc),AMIEEFlux(loc), lag)
      cloc = where(c eq max(c))
      cloc = cloc(0)
      metrics(2,0,iFile) = c(cloc)
      metrics(4,0,iFile) = MagPos(0,LocLeft) - MagPos(0,LocLeft + lag(cloc))
      c = c_correlate(J4Data(1,loc),J4EFlux(loc),0)
      metrics(3,0,iFile) = c(0)

      ; Eric would like a new metric - a sum of the fluxes around the main peak

      ; We want to sort of do this right.  Lets assume that the flux is relatively
      ; constant in longitude over a 15 degree "bin".  Then lets take 10 degrees on
      ; either side of the main peak in latitude.

      distance = abs(magpos(0,*) - magpos(0,LocLeft))
      dloc = where(distance le 10.0, dcount)
      efluxsum_dmsp = 0.0
      efluxsum_amie = 0.0

      for idist = 0, dcount-1 do begin

        long_width = sin(magpos(0,dloc(idist))*!pi/180.0) * 15.0 * !pi / 180.0 * 6371000.0
	lat_diff   = abs(magpos(0,dloc(idist)) - magpos(0,dloc(idist)-1))
        lat_width  = lat_diff * !pi / 180.0 * 6371000.0
	cell_size  = long_width*lat_width

        efluxsum_dmsp = efluxsum_dmsp + j4data(1,dloc(idist)) * cell_size
        efluxsum_amie = efluxsum_amie + amieeflux(dloc(idist)) * cell_size

      endfor

      if (metrics(2,0,iFile) gt 0.0) then begin
        metrics(12,0,iFile) = 100*(efluxsum_dmsp - efluxsum_amie)/ $
                             max([efluxsum_dmsp,efluxsum_amie])
      endif else metrics(12,0,iFile) = -9999.0

      c = c_correlate(J4Hall(loc),AMIEHall(loc), lag)
      cloc = where(c eq max(c))
      cloc = cloc(0)
      metrics(10,0,iFile) = c(cloc)
      metrics(11,0,iFile) = MagPos(0,LocLeft_Hall) - MagPos(0,LocLeft_Hall + lag(cloc))

    endif else begin
      ; Left
      locl = LocLeft + indgen(nAve*4+1) - nAve*2
      c = c_correlate(J4Data(1,locl),AMIEEFlux(locl), lag)
      cloc = where(c eq max(c))
      cloc = cloc(0)
      metrics(2,0,iFile) = c(cloc)
      metrics(4,0,iFile) = MagPos(0,LocLeft) - MagPos(0,LocLeft + lag(cloc))
      c = c_correlate(J4Data(1,locl),J4EFlux(locl),0)
      metrics(3,0,iFile) = c(0)

      locl = LocLeft_Hall + indgen(nAve*4+1) - nAve*2
      c = c_correlate(J4Hall(locl),AMIEHall(locl), lag)
      cloc = where(c eq max(c))
      cloc = cloc(0)
      metrics(10,0,iFile) = c(cloc)
      metrics(11,0,iFile) = MagPos(0,LocLeft_Hall) - MagPos(0,LocLeft_Hall + lag(cloc))

      distance = abs(magpos(0,*) - magpos(0,LocLeft))
      dloc = where(distance le 10.0, dcount)
      efluxsum_dmsp = 0.0
      efluxsum_amie = 0.0

      for idist = 0, dcount-1 do begin

        long_width = sin(magpos(0,dloc(idist))*!pi/180.0) * 15.0 * !pi / 180.0 * 6371000.0
	lat_diff   = abs(magpos(0,dloc(idist)) - magpos(0,dloc(idist)-1))
        lat_width  = lat_diff * !pi / 180.0 * 6371000.0
	cell_size  = long_width*lat_width

        efluxsum_dmsp = efluxsum_dmsp + j4data(1,dloc(idist)) * cell_size
        efluxsum_amie = efluxsum_amie + amieeflux(dloc(idist)) * cell_size

      endfor

      if (metrics(2,0,iFile) gt 0.0) then begin
        metrics(12,0,iFile) = 100*(efluxsum_dmsp - efluxsum_amie)/ $
                             max([efluxsum_dmsp,efluxsum_amie])
      endif else metrics(12,0,iFile) = -9999.0

      ; Right
      locr = LocRight + indgen(nAve*4+1) - nAve*2
      c = c_correlate(J4Data(1,locr),AMIEEFlux(locr), lag)
      cloc = where(c eq max(c))
      cloc = cloc(0)
      metrics(2,1,iFile) = c(cloc)
      metrics(4,1,iFile) = MagPos(0,LocRight) - MagPos(0,LocRight + lag(cloc))
      c = c_correlate(J4Data(1,locr),J4EFlux(locr),0)
      metrics(3,1,iFile) = c(0)

      locr = LocRight_Hall + indgen(nAve*4+1) - nAve*2
      c = c_correlate(J4Hall(locr),AMIEHall(locr), lag)
      cloc = where(c eq max(c))
      cloc = cloc(0)
      metrics(10,1,iFile) = c(cloc)
      metrics(11,1,iFile) = MagPos(0,LocRight_Hall) - MagPos(0,LocRight_Hall + lag(cloc))

      distance = abs(magpos(0,*) - magpos(0,LocRight))
      dloc = where(distance le 10.0, dcount)
      efluxsum_dmsp = 0.0
      efluxsum_amie = 0.0

      for idist = 0, dcount-1 do begin

        long_width = sin(magpos(0,dloc(idist))*!pi/180.0) * 15.0 * !pi / 180.0 * 6371000.0
	lat_diff   = abs(magpos(0,dloc(idist)) - magpos(0,dloc(idist)-1))
        lat_width  = lat_diff * !pi / 180.0 * 6371000.0
	cell_size  = long_width*lat_width

        efluxsum_dmsp = efluxsum_dmsp + j4data(1,dloc(idist)) * cell_size
        efluxsum_amie = efluxsum_amie + amieeflux(dloc(idist)) * cell_size

      endfor

      if (metrics(2,1,iFile) gt 0.0) then begin
        metrics(12,1,iFile) = 100*(efluxsum_dmsp - efluxsum_amie)/ $
                             max([efluxsum_dmsp,efluxsum_amie])
      endif else metrics(12,1,iFile) = -9999.0

    endelse

;    metrics(0,0,iFile) = J4Data(1,LocLeft) - AMIEEFlux(AMIELocLeft)

    if (nPeaks gt 1) then begin
      metrics(0,1,iFile) = 100*(J4Data(1,LocRight) - AMIEEFlux(AMIELocRight))/$
	max([J4Data(1,LocRight),AMIEEFlux(AMIELocRight)])
      metrics(1,1,iFile) = MagPos(0,LocRight) - MagPos(0,AMIELocRight)

      metrics(8,1,iFile) = 100*(J4Hall(LocRight_Hall) - AMIEHall(AMIELocRight_Hall))/$
	max([J4Hall(LocRight_Hall),AMIEHall(AMIELocRight_Hall)])
      metrics(9,1,iFile) = MagPos(0,LocRight_Hall) - MagPos(0,AMIELocRight_Hall)
    endif else metrics(*,1,iFile) = -9999.0

    ; We want to find the number of mag data points in the vacinity of
    ; the DMSP data.

    pn = pnBase + 0
    get_position, ppp, space, sizes, pn, pos
    pos([1,3]) = pos([1,3]) + fac*space

    plot, x,y,/nodata, xstyle = 5, ystyle = 5,	$
	xrange = [-mr,mr],yrange=[-mr,mr], pos = pos, /noerase

    CentralTime = mean(J4Time(loc))

    ; we only want to do this if we are comparing to Ahn formulation

    if (me_ ne 13) then begin

      DataFile = amie_file + '_data'

      openr,11, DataFile

      line = ''
      itime = intarr(6)

      done = 0

      while not done do begin

        readf, 11, line

        if eof(11) then begin
          print, 'EOF in datafile ',datafile
          done = 1
          type = -1
        endif

        if (strpos(line,'#TIME') gt -1) then begin
          readf, 11, itime
          c_a_to_r, itime, rtime_amie
          if (rtime_amie ge CentralTime) then begin
            done = 1
            type = 0
          endif
        endif

      endwhile

      metrics(5,*,iFile) = 0

      while (type eq 0) do begin

        readf,11,line

        if eof(11) then type = -1

        if (strpos(line,'#TIME') gt -1) then type = -1
        if (strpos(line,'#AHN') gt -1) then begin
          type = 1
          readf, 11, npts
          for i = 0, npts-1 do begin 
            readf, 11, lat, mlt, data

  	    if (90.0-lat le mr) then begin
              oplot, [(90.0-lat)*cos(mlt*!pi/12.0-!pi/2)],$
		  [(90.0-lat)*sin(mlt*!pi/12.0-!pi/2)], psym=4
 	    endif

            if (abs(mlt-MagPos(2,LocLeft)) lt 1.0) then $
              metrics(5,0,iFile) = metrics(5,0,iFile) + 1
            if (abs(mlt-MagPos(2,LocRight)) lt 1.0) then $
              metrics(5,1,iFile) = metrics(5,1,iFile) + 1

            ; special cases:

            if (mlt lt 1.0 and MagPos(2,LocLeft) gt 23.0) then $
              metrics(5,0,iFile) = metrics(5,0,iFile) + 1
            if (mlt lt 1.0 and MagPos(2,LocRight) gt 23.0) then $
              metrics(5,1,iFile) = metrics(5,1,iFile) + 1

            if (mlt gt 23.0 and MagPos(2,LocLeft) lt 1.0) then $
              metrics(5,0,iFile) = metrics(5,0,iFile) + 1
            if (mlt gt 23.0 and MagPos(2,LocRight) lt 1.0) then $
              metrics(5,1,iFile) = metrics(5,1,iFile) + 1

          endfor
        endif
      endwhile

    endif

    metrics(6,0,iFile) = abs(ae(CompIter(LocLeft),2))
    metrics(6,1,iFile) = abs(ae(CompIter(LocRight),2))

    MetricsTimes(iFile) = AMIETime(CompIter(LocLeft))
    MetricsTimes(iFile) = AMIETime(CompIter(LocRight))

    metrics(7,0,iFile) = MagPos(2,LocLeft)
    metrics(7,1,iFile) = MagPos(2,LocRight)

    close, 11

    oplot, [J4Time(LocLeft)-stime,J4Time(LocLeft)-stime], [0,1000]
    oplot, [J4Time(LocRight)-stime,J4Time(LocRight)-stime], [0,1000]

    pn = pnBase + 4
    get_position, ppp, space, sizes, pn, pos
    pos([1,3]) = pos([1,3]) + space - 0.01
    pos([1,3]) = pos([1,3]) + fac*space

    iter = mean(CompIter(loc))

    mini = 0.0
    maxi = max(AMIEData(iter,ha_,*,*))
    range = maxi
    dc      = 10.0^fix(alog10(range/20.0))
    factor  = 0.0
    while (range gt dc*20.0*factor) do factor=factor+0.1
    dc = factor*dc
    levels = findgen(21)*dc

    ;readct, ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
    ;readct, ncolors, getenv("IDL_EXTRAS")+"white_red.ct"
    ;clevels = (ncolors-1)*findgen(21)/20.0 + 1

    contour, reform(AMIEData(iter,ha_,*,LatLoc)),x,y,/follow, 		$
	xstyle = 5, ystyle = 5,						$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase, /cell_fill, c_color = clevels

  ;  xyouts, pos(0)-0.01,pos(3)-0.01, amie_time_1(n),charsize=0.9, /norm

    plotmlt, mr, /no06, /no12
    oplot, J4x, J4y
    xyouts, J4x(0), J4y(0), 'S'

    pn = pnBase + 5
    get_position, ppp, space, sizes, pn, pos1

    pn = pnBase + 7
    get_position, ppp, space, sizes, pn, pos2

    pos(0) = pos1(0)
    pos(1) = pos1(1)
    pos(2) = pos2(2)
    pos(3) = pos1(3)

    pos([1,3]) = pos([1,3]) + space - 0.01
    pos([1,3]) = pos([1,3]) + fac*space

    plot, J4Time(loc)-stime,J4hall(loc), pos = pos,	 		$
	xstyle = 1, /noerase,		$
	xtickname = xtickname, xtickv=xtickvalue, 			$
	xticks = xtickn, xminor = xminor, xtitle = xtitle,		$
	xrange = [s_time_range, e_time_range], 				$
	ytickname = strarr(20)+' ', yrange = mm([J4hall(loc),AMIEHall(loc)])
    axis, yax=1, charsize=0.9

    xyouts, pos(0)-0.005, (pos(1)+pos(3))/2.0, 'Hall Conductance (mhos)', $
  	/norm, alignment = 0.5, orient = 90, charsize = 0.75

    oplot, J4Time(loc)-stime, AMIEHall(loc), linestyle = 2
;    oplot, J4Time(loc)-stime, J4AAveE(loc), linestyle = 1

  endif else begin

    print, "There are no data points in this hemisphere."

  endelse

endfor

closedevice

openw, 2, amie_file+'_metrics'

printf,2,nFiles
printf,2,nMetrics
printf,2,'% Diff in E-Flux'
printf,2,'Latitude Diff'
printf,2,'CC - w/AMIE'
printf,2,'CC - w/Ave'
printf,2,'CC Lat Diff'
printf,2,'nMag'
printf,2,'abs(AL)'
printf,2,'MLT'
printf,2,'% Diff Hall'
printf,2,'Latitude Diff Hall'
printf,2,'CC - w/AMIE Hall'
printf,2,'CC Lat Diff Hall'
printf,2,'DMSP Tenergy - AMIE Tenergy'
printf,2,metrics
printf,2,MetricsTimes

close,2

end
