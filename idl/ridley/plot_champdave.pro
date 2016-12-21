
ccmc = 1

if n_elements(plotdiff) eq 0 then plotdiff = 0
plotdiff = fix(ask('plot GITM contour difference: ',tostr(plotdiff)))
reread = 1
if (n_elements(GitmDensity) gt 0) then begin
    answer = ask('whether to re-read data','n')
    if (strpos(mklower(answer),'n') gt -1) then reread = 0
endif

if (reread) then begin

    if (n_elements(dir) eq 0) then dir = '.'
    dir = ask('directory',dir)
    if plotdiff then begin
        if (n_elements(dirbase) eq 0) then dirbase = '.'
        dirbase = ask('base directory', dirbase)
        filelist_base = file_search(dirbase+'/cham*.bin')
    endif

    filelist = file_search(dir+'/cham*.bin')
    
    file = filelist(0)
    length = strpos(file,'.bin')
    ls = length-13
    
    yr = strmid(file,ls,2)
    iYear = tostr(2000+fix(yr))
    iMonth = strmid(file,ls+2,2)
    iDay = strmid(file,ls+4,2)

    champdir = '/bigdisk1/Data6/Data/CHAMP/'+tostr(iyear)+'/'
;    itime = [fix(iyear), fix(imonth), fix(iday), 0,0,0]
;    date = iyear+'-'+imonth+'-'+iday
;    realdate = date_conv(date,'r')
    
    doy0 = jday(fix(iyear), fix(imonth), fix(iday))
    ndays = 1

    thermo_readsat, filelist, data, time, nTimes, Vars, nAlts, nSats, nTimes
    if plotdiff then begin
        thermo_readsat, filelist_base, data_base, time_base, nTimes_base, Vars_base, nAlts_base,$
          nSats_base, nTimes_base
        GitmAltsbase = reform(data_base(0,0,2,*))/1000.0
        GitmRhobase  = reform(data_base(0,*,3,*))
        GitmLonsbase = reform(data_base(0,*,0,*))*180.0/!pi
        GitmLatsbase = reform(data_base(0,*,1,*))*180.0/!pi
    endif

    GitmAlts     = reform(data(0,0,2,*))/1000.0
    GitmRho      = reform(data(0,*,3,*))
    GitmElectron = reform(data(0,*,33,*))
    GitmLons     = reform(data(0,*,0,*))*180.0/!pi
    GitmLats     = reform(data(0,*,1,*))*180.0/!pi

    c_r_to_a, itime, time(0)
    itime(3:5) = 0
    ndays = round((time(ntimes-1)-time(0))/3600. /24.)
    ndays = fix((time(ntimes-1)-time(0))/3600. /24.)+1
;    if (ndays eq 0) then ndays = 1

    c_a_to_r, itime, basetime
    hour = (time/3600.0 mod 24.0) + fix((time-basetime)/(24.0*3600.0))*24.0
    localtime = (reform(GitmLons(*,0))/15.0 + hour) mod 24.0

    ChampDensity    = fltarr(nTimes)
    GitmDensity     = fltarr(nTimes)
    GitmeDensity     = fltarr(nTimes)
    GitmHMF2         = fltarr(nTimes)
    GitmHMF2u        = fltarr(nTimes)
    GitmHMF2l        = fltarr(nTimes)
    GitmNMF2         = fltarr(nTimes)
    GitmDensityHigh = fltarr(nTimes)
    GitmDensityLow  = fltarr(nTimes)
    GitmHeight      = fltarr(nTimes)
    GitmHeightHigh  = fltarr(nTimes)
    GitmHeightLow   = fltarr(nTimes)
    if plotdiff then begin
        GitmDensityBase     = fltarr(nTimes)
        GitmDensityHighBase = fltarr(nTimes)
        GitmDensityLowBase  = fltarr(nTimes)
    endif
    ChampAltitude   = fltarr(nTimes)

    nChampMax = 50000L
    ChampPosition = fltarr(3,nChampMax)
    ChampTime  = dblarr(nChampMax)
    MassDensity = fltarr(nChampMax)
    ChampWind = fltarr(nChampMax)
    ChampLocalTime=fltarr(nChampMax)
    t = ' '
    line = 0L
    for cday = 0, ndays - 1 do begin
        doy = doy0 + cday
        champ_file_a = champdir+'Density_3deg_'+yr+'_'+ $
          chopr('00'+tostr(doy),3)+'.ascii'
   ;     champ_file_w = champdir+'Wind_3deg_'+yr+'_'+tostr(doy)+'.ascii'
        
        close,/all
        print, "Reading Champ File : ",champ_file_a
        openr,1,champ_file_a
    ;    openr,2,champ_file_w
        readf,1,t
        tarr = strsplit(t,/extract)
        version = float(tarr(1))
        print, version
        readf,1,t
    ;    readf,2,t

        while (not eof(1)) do begin
            readf,1,t
            tarr = strsplit(t,/extract)
            year = fix(tarr(0))
            day = fix(tarr(1))
            seconds = float(tarr(2))
            lat =float(tarr(4))
            if (version lt 2.2) then begin
               long = float(tarr(5))
               height = float(tarr(6))
               chlocaltime = float(tarr(7))
               density = float(tarr(8))
               density400 = float(tarr(9))
               density410 =float(tarr(10))
            endif else begin 
               long = float(tarr(5))
               height = float(tarr(6))
               chlocaltime = float(tarr(7))
               density = float(tarr(11))
               density400 = float(tarr(12))
               density410 = float(tarr(13))
            endelse

            itime = [Year, 1, Day, 0,0,0]
            c_a_to_r, iTime, BaseTime

            ChampTime(line) = seconds+ basetime
            ChampPosition(0,line) = long
            ChampPosition(1,line) = lat
            ChampPosition(2,line) = height
;            ChampPosition(2,line) = 400.0
            MassDensity(line) = density
            ChampLocalTime(line) = chlocaltime
;            ChampWind(line) = wind
           
            line = line + 1
        endwhile
        
        close,1,2
    endfor

    for iTime = 0, nTimes-1 do begin

        dt = abs(time(iTime)-ChampTime)
        loc = where(dt eq min(dt))

        i = loc(0)
  
        ChampDensity(iTime)  = MassDensity(i)/1.e-12
        ChampAltitude(iTime) = ChampPosition(2,i)
        
        loc = where(GitmAlts gt ChampAltitude(iTime))
        i = loc(0)
        x = (ChampAltitude(iTime) - GitmAlts(i-1)) / $
          (GitmAlts(i) - GitmAlts(i-1))
        GitmDensity(iTime) = exp((1.0 - x) * alog(GitmRho(iTime,i-1)) + $
          (      x) * alog(GitmRho(iTIme,i)))

        GitmeDensity(iTime) = exp((1.0 - x) * alog(GitmElectron(iTime,i-1)) + $
          (      x) * alog(GitmElectron(iTIme,i)))

        ed = reform(GitmElectron(iTime,*))

        l = where(ed eq max(ed))
        GitmNMF2(iTime) = ed(l(0))
        GitmHMF2(iTime) = GitmAlts(l(0))

        l = where(ed gt 0.9*max(ed),c_ed)
        GitmHMF2l(iTime) = GitmAlts(l(0))
        GitmHMF2u(iTime) = GitmAlts(l(c_ed-1))

        if plotdiff then  GitmDensityBase(iTime) = exp((1.0 - x) * alog(GitmRhoBase(iTime,i-1)) + $
          (      x) * alog(GitmRho(iTIme,i)))

        GitmHeight(iTime) = ChampAltitude(iTime)

        h = (GitmAlts(i+1) - GitmAlts(i-1))/2.0

;        print, ChampAltitude(iTime), ChampDensity(iTime), GitmDensity(iTime)/1.0e-12, $
;          (1.0 - x)*GitmAlts(i-1)+x*GitmAlts(i)

        loc = where(GitmAlts gt ChampAltitude(iTime)+h)
        i = loc(0)
        x = ((ChampAltitude(iTime)+h) - GitmAlts(i-1)) / $
          (GitmAlts(i) - GitmAlts(i-1))
        GitmDensityHigh(iTime) = (1.0 - x) * GitmRho(iTime,i-1) + $
          (      x) * GitmRho(iTIme,i)
        if plotdiff then GitmDensityHighBase(iTime) = (1.0 - x) * GitmRhoBase(iTime,i-1) + $
          (      x) * GitmRhoBase(iTIme,i)

        loc = where(GitmAlts gt ChampAltitude(iTime)-h)
        i = loc(0)
        x = ((ChampAltitude(iTime)-h) - GitmAlts(i-1)) / $
          (GitmAlts(i) - GitmAlts(i-1))
        GitmDensityLow(iTime) = (1.0 - x) * GitmRho(iTime,i-1) + $
          (      x) * GitmRho(iTIme,i)
        if plotdiff then GitmDensityLowBase(iTime) = (1.0 - x) * GitmRhoBase(iTime,i-1) + $
          (      x) * GitmRhoBase(iTIme,i)
    endfor
    
    GitmDensity = GitmDensity*1.0e12
    ChampLocalTime = ChampLocalTime(0:nTimes-1)

    if plotdiff then GitmDensityBase = GitmDensityBase*1.0e12
    
    
    GitmAve  = fltarr(nTimes)
    if plotdiff then GitmAveBase  = fltarr(nTimes)
    ChampAve = fltarr(nTimes)

    for iTime = 0, nTimes-1 do begin

        loc = where(abs(time-time(iTime)) lt 45.0*60.0, count)
        if (count gt 0) then begin
            GitmAve(iTime) = mean(GitmDensity(loc))
            if plotdiff then GitmAveBase(iTime) = mean(GitmDensityBase(loc))
            ChampAve(iTime) = mean(ChampDensity(loc))
        endif

    endfor

    nOrbits = 0
        
    test = where(localtime gt 4.0 and localtime lt 8.0, nTest)
    if (float(ntest)/nTimes lt 0.25) then begin    
       ; This means we NOT in a dawn-dusk orbit
       day   = where(localtime gt 6.0 and localtime lt 18.0,nPtsDay)
       night = where(localtime lt 6.0 or  localtime gt 18.0,nPtsNight)
    endif else begin
       ; If we are in Dawn/Dusk, define day and night differently
       day   = where(localtime lt 12.0,nPtsDay)
       night = where(localtime gt 12.0,nPtsNight)
    endelse

    for i = 1, nPtsDay-1 do begin
        
        if (day(i)-day(i-1) gt 1) then begin
            if (nOrbits eq 0) then begin
                DayOrbitStart = day(i)
                DayOrbitEnd   = day(i-1)
            endif else begin
                if (day(i)-day(i-1) gt 4) then begin
                    DayOrbitStart = [DayOrbitStart,day(i)]
                    DayOrbitEnd   = [DayOrbitEnd  ,day(i-1)]
                endif
            endelse
            if (day(i)-day(i-1) gt 4) then nOrbits = nOrbits+1
        endif

    endfor

    nY = max(DayOrbitStart - DayOrbitEnd)

    print, "day : ", nY, nOrbits

    if (nOrbits gt 0) then begin

        xDay = fltarr(nOrbits,nY)
        yDay = fltarr(nOrbits,nY)
        vDay = fltarr(nOrbits,nY)
        cDay = fltarr(nOrbits,nY)
        
        iOrbit = 0
        iY = 0
        iFound = 0
        for i = 1, nPtsDay-1 do begin
        
            if (day(i)-day(i-1) gt 1) then begin
                if (day(i)-day(i-1) gt 4) then begin
                    iOrbit = iOrbit+1
                    iY = 0
                endif
                iFound = 1
            endif else iY = iY + 1
        
            if (iFound) then begin
                xDay(iOrbit-1, iY) = hour(DayOrbitStart(iOrbit-1))
                yDay(iOrbit-1, iY) = GitmLats(Day(i))
                vDay(iOrbit-1, iY) = GitmDensity(Day(i))
                if plotdiff then $
                  vDay(iOrbit-1, iY) = $
                  (GitmDensity(Day(i))-GitmDensityBase(Day(i))) /$
                  GitmDensity(Day(i))
                cDay(iOrbit-1, iY) = ChampDensity(Day(i))
            endif

        endfor

        for iOrbit = 0, nOrbits-2 do begin
            l = where(xday(iOrbit,*) eq 0,c)
            if (c gt 0) then begin
                for j = 0,c-1 do begin
                    xDay(iOrbit,l(j)) = xDay(iOrbit,l(j)-1)
                    yDay(iOrbit,l(j)) = yDay(iOrbit,l(j)-1)
                    vDay(iOrbit,l(j)) = vDay(iOrbit,l(j)-1)
                    cDay(iOrbit,l(j)) = cDay(iOrbit,l(j)-1)
                endfor
            endif
        endfor

        DayGood = 1

    endif else DayGood = 0

    nNOrbits = 0
    
    for i = 1, nPtsNight-1 do begin

        if (night(i)-night(i-1) gt 1) then begin
            if (nNorbits eq 0) then begin
                NightorbitStart = night(i)
                NightOrbitEnd   = night(i-1)
            endif else begin
                if (night(i)-night(i-1) gt 4) then begin
                    NightOrbitStart = [NightOrbitStart,night(i)]
                    NightOrbitEnd   = [NightOrbitEnd  ,night(i-1)]
                endif
            endelse
            if (night(i)-night(i-1) gt 4) then nNorbits = nNorbits+1
        endif

    endfor

    nY = max(NightOrbitStart - NightOrbitEnd)
    
    print, "night : ", nY, nOrbits

    if (nOrbits gt 0) then begin

        xNight = fltarr(nNorbits,nY)
        yNight = fltarr(nNorbits,nY)
        vNight = fltarr(nNorbits,nY)
        cNight = fltarr(nNorbits,nY)

        iNorbit = 0
        iY = 0
        iFound = 0
        for i = 1, nPtsNight-1 do begin

           print, night(i)-night(i-1)
            if (night(i)-night(i-1) gt 1) then begin
               print, '> 1'
                if (night(i)-night(i-1) gt 4) then begin
                    iNorbit = iNorbit+1
                    iY = 0
                endif
                if (iNorbit gt 0) then iFound = 1
            endif else iY = iY + 1

            if (iFound) then begin
                xNight(iNorbit-1, iY) = hour(NightorbitStart(iNorbit-1))
                yNight(iNorbit-1, iY) = GitmLats(Night(i))
                vNight(iNorbit-1, iY) = GitmDensity(Night(i))
                if plotdiff then $
                  vNight(iOrbit-1, iY) = (GitmDensity(Night(i))-$
                                          GitmDensityBase(Night(i))) /$
                  GitmDensity(Night(i))
                cNight(iNorbit-1, iY) = ChampDensity(Night(i))
            endif

        endfor

        for iOrbit = 0, nNorbits-1 do begin
            l = where(xNight(iOrbit,*) eq 0,c)
            if (c gt 0) then begin
                for j = 0,c-1 do begin
                    xNight(iOrbit,l(j)) = xNight(iOrbit,l(j)-1)
                    yNight(iOrbit,l(j)) = yNight(iOrbit,l(j)-1)
                    vNight(iOrbit,l(j)) = vNight(iOrbit,l(j)-1)
                    cNight(iOrbit,l(j)) = cNight(iOrbit,l(j)-1)
                endfor
            endif
        endfor

        NightGood = 1

    endif else NightGood = 0

endif

print, "daynight : ", daygood, nightgood

yrange = mm([ChampDensity,GITMDensity])
;yrange = [0.0,40.0]

ppp = 2
space = 0.1
pos_space, ppp, space, sizes, ny = ppp
    
get_position, ppp, space, sizes, 0, pos, /rect
pos(0) = pos(0) + 0.05
pos(2) = pos(2) - 0.05

get_position, ppp, space, sizes, 1, pos2, /rect
pos2(0) = pos2(0) + 0.05
pos2(2) = pos2(2) - 0.05

stime = min(time)
etime = max(time)
time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn

p = strpos(dir,'/')
if p eq -1 then p = strlen(dir)
run = strmid(dir,0,p)

psfile_raw = 'champ_raw.ps'
psfile_ave = 'champ_ave.ps'
psfile_2dd = 'champ_2dd.ps'
psfile_2dn = 'champ_2dn.ps'

;--------------------------------------------------------------------
; Raw
;--------------------------------------------------------------------

if not plotdiff then begin
    setdevice, psfile_raw, 'p', 5
    
    plot, time-stime, ChampDensity, yrange = yrange, pos = pos, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = 'Mass Density (10!E-12!N kg/m!E3!N)',   $
      thick = 3
    oplot, time-stime, gitmdensity, linestyle = 2, thick = 3
    
    t1 = (etr-btr)*0.05
    t2 = (etr-btr)*0.10
    t3 = (etr-btr)*0.11
    
    oplot, [t1,t2], max(yrange) - [2.0,2.0], thick = 3, linestyle = 2
    xyouts, t3, max(yrange) - 2.0, 'GITM'
    
    oplot, [t1,t2], max(yrange) - [4.0,4.0], thick = 3
    xyouts, t3, max(yrange) - 4.0, 'CHAMP'
    
    
    plot, time-stime, ChampDensity-gitmdensity, $
      yrange = yrange-max(yrange)/2, pos = pos2, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = 'Mass Density (10!E-12!N kg/m!E3!N)',   $
      thick = 3, /noerase, ystyle = 1
    
    rmse = sqrt(mean((ChampDensity-gitmdensity)^2))
    rmsd = sqrt(mean((ChampDensity)^2))
    nrms = rmse/rmsd * 100.0
    
    pdif = mean((ChampDensity-gitmDensity)/ChampDensity) * 100.0
    
    srms = +' (nRMS: '+string(nrms,format = '(f5.1)')+'%, '
    srms = srms+string(pdif,format = '(f5.1)')+'% Difference)'
    
    oplot, [t1,t2], max(yrange-max(yrange)/2) - [2.0,2.0], thick = 3
    xyouts, t3, max(yrange-max(yrange)/2) - 2.0, 'CHAMP - GITM'+srms
    
    oplot, [btr,etr], [0.0,0.0], linestyle = 1
    
    xyouts, 0.0, -0.02, dir, /norm, charsize = 0.8
    
    closedevice
    
;--------------------------------------------------------------------
; Ave
;--------------------------------------------------------------------
    
    setdevice, psfile_ave, 'p', 5
    
    plot, time-stime, ChampAve, yrange = yrange, pos = pos, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = 'Mass Density (10!E-12!N kg/m!E3!N)',   $
      thick = 3
    oplot, time-stime, GitmAve, linestyle = 2, thick = 3
    
    oplot, [t1,t2], max(yrange) - [2.0,2.0], thick = 3, linestyle = 2
    xyouts, t3, max(yrange) - 2.0, 'GITM'
    
    oplot, [t1,t2], max(yrange) - [4.0,4.0], thick = 3
    xyouts, t3, max(yrange) - 4.0, 'CHAMP'
    
    plot, time-stime, ChampAve-gitmave, $
      yrange = yrange-max(yrange)/2, pos = pos2, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = 'Mass Density (10!E-12!N kg/m!E3!N)',   $
      thick = 3, /noerase, ystyle = 1
    
    rmse = sqrt(mean((ChampAve-gitmave)^2))
    rmsd = sqrt(mean((ChampAve)^2))
    nrms = rmse/rmsd * 100.0
    
    pdif = mean((ChampAve-gitmave)/ChampAve) * 100.0
    
    srms = +' (nRMS: '+string(nrms,format = '(f5.1)')+'%, '
    srms = srms+string(pdif,format = '(f5.1)')+'% Difference)'
    
    oplot, [t1,t2], max(yrange-max(yrange)/2) - [2.0,2.0], thick = 3
    xyouts, t3, max(yrange-max(yrange)/2) - 2.0, 'CHAMP - GITM'+srms
    
    oplot, [btr,etr], [0.0,0.0], linestyle = 1
    
    xyouts, 0.0, -0.02, dir, /norm, charsize = 0.8
    
    closedevice
endif

if (DayGood) then begin

;--------------------------------------------------------------------
; Day 2d
;--------------------------------------------------------------------

    setdevice, psfile_2dd, 'p', 5

    makect, 'wyr'

    nX = n_elements(cDay(*,0))
    nY = n_elements(cDay(0,*))

    levels = findgen(61) * 15.0/60.0
    linelevels = findgen(7) * 15.0/6.0

    ytickv = [-90,-60,-30,0,30,60,90]

    contour, cday(0:nX-2,0:nY-2), $
      xday(0:nX-2,0:nY-2)*3600.0, yday(0:nX-2,0:nY-2), $
      /fill, pos = pos, yrange = [-90,90], ystyle = 1, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytickv = ytickv, yticks = 7, yminor = 6, $
      ytitle = 'Latitude (Deg)',   $
      thick = 3, levels = levels

    ctpos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
    plotct,254,ctpos,mm(levels),$
      'CHAMP Dayside Mass Density (10!E-12!N kg/m!E3!N)',/right
    

    contour, vday(0:nX-2,0:nY-2), $
      xday(0:nX-2,0:nY-2)*3600.0, yday(0:nX-2,0:nY-2), $
      /fill, pos = pos2, yrange = [-90,90], ystyle = 1, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytickv = ytickv, yticks = 7, yminor = 6, $
      ytitle = 'Latitude (Deg)',   $
      thick = 3, /noerase, levels = levels

    contour, cday(0:nX-2,0:nY-2), $
      xday(0:nX-2,0:nY-2)*3600.0, yday(0:nX-2,0:nY-2), $
      pos = pos2, yrange = [-90,90], ystyle = 1, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytickv = ytickv, yticks = 7, yminor = 6, $
      thick = 3, /noerase, levels = linelevels, /follow, c_linestyle = 1

    ctpos = [pos2(2)+0.01,pos2(1),pos2(2)+0.03,pos2(3)]
    plotct,254,ctpos,mm(levels),$
      'GITM Dayside Mass Density (10!E-12!N kg/m!E3!N)',/right
    
    xyouts, 0.0, -0.02, dir, /norm, charsize = 0.8

    closedevice

endif

if (NightGood) then begin

;--------------------------------------------------------------------
; Night 2d
;--------------------------------------------------------------------

    setdevice, psfile_2dn, 'p', 5

    makect, 'wyr'

;    contour, cnight(0:nX-2,0:nY-2), $
;      xnight(0:nX-2,0:nY-2)*3600.0, ynight(0:nX-2,0:nY-2), $
    contour, cnight, $
      xnight*3600.0, ynight, $
      /fill, pos = pos, yrange = [-90,90], ystyle = 1, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytickv = ytickv, yticks = 7, yminor = 6, $
      ytitle = 'Latitude (Deg)',   $
      thick = 3, levels = levels, /follow

    ctpos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
    plotct,254,ctpos,mm(levels),$
      'CHAMP Nightside Mass Density (10!E-12!N kg/m!E3!N)',/right

;    contour, vnight(0:nX-2,0:nY-2), $
;      xnight(0:nX-2,0:nY-2)*3600.0, ynight(0:nX-2,0:nY-2), $
    contour, vnight, $
      xnight*3600.0, ynight, $
      /fill, pos = pos2, yrange = [-90,90], ystyle = 1, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytickv = ytickv, yticks = 7, yminor = 6, $
      ytitle = 'Latitude (Deg)',   $
      thick = 3, /noerase, levels = levels

;    contour, cnight(0:nX-2,0:nY-2), $
;      xnight(0:nX-2,0:nY-2)*3600.0, ynight(0:nX-2,0:nY-2), $
    contour, cnight, $
      xnight*3600.0, ynight, $
      pos = pos2, yrange = [-90,90], ystyle = 1, $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytickv = ytickv, yticks = 7, yminor = 6, $
      thick = 3, /noerase, levels = linelevels, /follow, c_linestyle = 1

    ctpos = [pos2(2)+0.01,pos2(1),pos2(2)+0.03,pos2(3)]
    plotct,254,ctpos,mm(levels),$
      'GITM Nightside Mass Density (10!E-12!N kg/m!E3!N)',/right

    xyouts, 0.0, -0.02, dir, /norm, charsize = 0.8

    closedevice

endif

if (ccmc) then begin

   stime = min(time)
   c_r_to_a, itime, stime
   c_a_to_ymd, itime, ymd
   outfile1 = 'champ_density_'+ymd+'.ccmc'
   outfile2 = 'champ_electron_'+ymd+'.ccmc'
   outfile3 = 'champ_nmf2_'+ymd+'.ccmc'
   outfile4 = 'champ_hmf2_'+ymd+'.ccmc'

   openw,1,outfile1
   openw,2,outfile2
;   openw,3,outfile3
;   openw,4,outfile4

;   printf,1, '#year doy hh mm ss lat lon rho'
;   printf,1, '#         hr mi s  deg deg kg/m3'

;   printf,2, '#year doy hh mm ss lat lon eden'
;   printf,2, '#         hr mi s  deg deg cm-3'
   
   for iPt = 0, nTimes-1 do begin

      c_r_to_a, itime, time(iPt)
      jd = jday(itime(0), itime(1), itime(2))
      printf,1, itime(0), jd, itime(3), itime(4), itime(5), gitmlats(iPt,0), gitmlons(iPt,0), gitmheight(iPt), gitmdensity(iPt)/1.0e12, format = '(i5,4i4,3f7.2,e12.4)'

      printf,2, itime(0), jd, itime(3), itime(4), itime(5), gitmlats(iPt,0), gitmlons(iPt,0), gitmheight(iPt), gitmedensity(iPt)/1.0e6, format = '(i5,4i4,3f7.2,e12.4)'

;      printf,3, itime(0), jd, itime(3), itime(4), itime(5), gitmlats(iPt,0), gitmlons(iPt,0), gitmnmf2(iPt)/1.0e6, format = '(i5,4i4,2f7.2,e12.4)'

;      printf,4, itime(0), jd, itime(3), itime(4), itime(5), gitmlats(iPt,0), gitmlons(iPt,0), gitmhmf2(iPt), gitmhmf2u(iPt), gitmhmf2l(iPt), format = '(i5,4i4,2f7.2,3f12.2)'

   endfor

   close,1,2,3,4

endif


;if (not plotdiff and nOrbits gt 0) then begin
;    ilats = [0,-10,-20,-30]
;    dl = 5
;    cplot = fltarr(n_elements(ilats),nOrbits)
;    gplot = cplot
;    ctime = cplot
;    cplotnorm = cplot
;    gplotnorm = gplot
;    yplot = fltarr(nOrbits)
;    setdevice,run+'line.ps','p',5,.95
;    ppp = n_elements(ilats)+1
;    space = 0.01
;    pos_space, ppp, space, sizes, ny = ppp
;    
;    for il = 0, n_elements(ilats) - 1 do begin
;        for iorbit = 0, norbits - 1 do begin
;            ind = where(yday(iorbit,*) ge ilats(il) and yday(iorbit,*) lt ilats(il)+dl)
;            cplot(il,iorbit) = max(cday(iorbit,ind))
;            gplot(il,iorbit) = max(vday(iorbit,ind))
;            it = where(champdensity eq cday(iorbit,ind(0)))
;            ctime(il,iorbit) = time(it) 
;        endfor
;    endfor
;       
;    for il = 0, n_elements(ilats) - 1 do begin
;        
;        cplotnorm(il,*) = (cplot(il,*)-min(cplot(il,*)))/$
;          (max(cplot(il,*)) - min(cplot(il,*)))*1.0 
;        gplotnorm(il,*) = (gplot(il,*)-min(gplot(il,*)))/$
;          (max(gplot(il,*)) - min(gplot(il,*)))*1.0 
;        
;        get_position, ppp, space, sizes, il, pos, /rect
;        pos(0) = pos(0) + 0.05
;        pos(2) = pos(2) - 0.05
;        
;        mini = -0.1             ;.85*min([cplot,gplot])
;        maxi = 1.1              ;1.1*max([cplot,gplot])
;        
;        if il eq n_elements(ilats) - 1 then begin
;            plot,ctime(il,*)-stime,cplotnorm(il,*),xrange = [0,24*60.*60],xstyle = 1, yrange = [mini,maxi], $
;              xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, xminor = xminor, xticks = xtickn, $
;              charsize = 1.2, thick = 3,/noerase,pos = pos,/nodata,$
;              ytitle = '                                                                                               Normalized Mass Density (x 10!U-12!N kg/m!U3!N)',ystyle = 1
;        endif else begin
;            plot,ctime(il,*)-stime,cplot(il,*),xrange = [0,24*60.*60],xstyle = 1, yrange = [mini,maxi], $
;              xtickname = strarr(10)+' ', xtickv = xtickv, xminor = xminor, xticks = xtickn, $
;              charsize = 1.2,/noerase,pos = pos,/nodata,ystyle = 1
;        endelse
;        
;        oplot, ctime(il,*)-stime,cplotnorm(il,*),psym = 4,color = 60,thick = 5
;        oplot, ctime(il,*)-stime,gplotnorm(il,*), thick = 5, psym= 1,color = 230
;        xyouts, pos(0) + .03 ,pos(3)-.03, tostr((ilats(il)+dl)/2.)+'!Uo!N Lat',/norm
;        
;    endfor
;    
;    legend,[' Champ',' GITM'],psym = [4,1],color = [60,230],thick = 5,pos = [pos(0)+.03,pos(1) - .02]$
;      ,/norm, box = 0
;    
;    closedevice
;
;endif
end

