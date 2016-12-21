
ccmc = 1

if n_elements(plotdiff) eq 0 then plotdiff = 0
plotdiff = fix(ask('plot GITM contour difference: ',tostr(plotdiff)))
reread = 1
if (n_elements(GitmDensity) gt 0) then begin
    answer = ask('whether to re-read data','n')
    if (strpos(mklower(answer),'n') gt -1) then reread = 0
endif

if (reread) then begin

;    if (n_elements(dir) eq 0) then dir = '.'
;    dir = ask('directory',dir)
;    if plotdiff then begin
;        if (n_elements(dirbase) eq 0) then dirbase = '.'
;        dirbase = ask('base directory', dirbase)
;        filelist_base = file_search(dirbase+'/cham*.bin')
;    endif

    filelist = file_search('cham*.bin')
    gitm_read_bin, filelist, data, time, nVars, Vars, version, /skip
    nTimes = n_elements(time)
    nAlts = n_elements(data(0,0,0,0,*))
    
    ; Need to reform the data
    newdata = dblarr(1,nTimes,nVars,nAlts)
    newdata(0,*,*,*) = data(*,*,0,0,*)
    data = newdata

    ; Ok, it seems that sometimes the appended files are a bit on
    ; the screwed up side.  Crap. Need to re-arrange the data if it is
    ; screwed up.  First test to see if it is messed up:
    dt = time(1:nTimes-1) - time(0:nTimes-2)
    if (min(dt) lt 0.0) then begin
       ; Data is screwed up!
       print, 'Order of data seems to be messed up.  Correcting.'
       timesave = time
       sort, time, order
       datafixed = fltarr(1,nTimes,nVars,nAlts)
       for i=0L,nTimes-1 do datafixed(0,i,*,*) = data(0,order(i),*,*)
       data = datafixed
    endif    

    nSats = 1

    c_r_to_a, itime, time(0)

    iYear = itime(0)
    iMonth = itime(1)
    iDay = itime(2)

    champdir = '/bigdisk1/Data6/Data/CHAMP/'+tostr(iyear)+'/'
    
    doy0 = jday(fix(iyear), fix(imonth), fix(iday))
    ndays = 1

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
    c_r_to_a, ietime, max(time)
    itime(3:5) = 0
    ietime(3:5) = [23,59,59]

    c_a_to_r,  itime, stime
    c_a_to_r, ietime, etime

    ndays = round((etime-stime)/3600.0/24.0)

;    ndays = round((time(ntimes-1)-time(0))/3600. /24.)
;    ndays = fix((time(ntimes-1)-time(0))/3600. /24.)+1


;    if (ndays eq 0) then ndays = 1

    c_a_to_r, itime, basetime
    hour = (time/3600.0 mod 24.0) + fix((time-basetime)/(24.0*3600.0))*24.0
    localtime = (reform(GitmLons(*,0))/15.0 + hour) mod 24.0

    ChampDensity    = fltarr(nTimes)
    ChampLats       = fltarr(nTimes)
    ChampLTs        = fltarr(nTimes)
    GitmDensity     = fltarr(nTimes)
    GitmeDensity    = fltarr(nTimes)
    GitmHMF2        = fltarr(nTimes)
    GitmHMF2u       = fltarr(nTimes)
    GitmHMF2l       = fltarr(nTimes)
    GitmNMF2        = fltarr(nTimes)
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

    nChampMax = 75000L
    ChampPosition = fltarr(3,nChampMax)
    ChampTime  = dblarr(nChampMax)
    MassDensity = fltarr(nChampMax)
    ChampWind = fltarr(nChampMax)
    ChampLocalTime=fltarr(nChampMax)
    t = ' '
    line = 0L
    for cday = 0, ndays - 1 do begin

       doy = doy0 + cday
       champ_file_a = champdir+'Density_3deg_'+tostr((iYear mod 100),2)+'_'+ $
                      chopr('00'+tostr(doy),3)+'.ascii'
   ;     champ_file_w = champdir+'Wind_3deg_'+yr+'_'+tostr(doy)+'.ascii'
        
       close,/all
       print, "Reading Champ File : ",champ_file_a

       IsGoodFile = 0
       on_ioerror, badchampfile

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
        
       IsGoodFile = 1
       badchampfile: if (not IsGoodFile) then print, "Bad CHAMP file!"

       close,1,2

    endfor

    iTimeSave = 0L

    for iTime = 0L, nTimes-1 do begin

;        dt = abs(time(iTime)-ChampTime)
;        loc = where(dt eq min(dt))
;
;        i = loc(0)
;  
;        ChampDensity(iTime)  = MassDensity(i)/1.e-12
;        ChampAltitude(iTime) = ChampPosition(2,i)
;        ChampLats(iTime)      = ChampPosition(1,i)
;        ChampLTs(iTime)       = ChampLocalTime(i)

        dt = time(iTime) - ChampTime(0:line-1)
        loc = where(dt lt 0, c)
        champdt = 0.0
        if (c eq 0) then begin
           ; haven't reached the start of the file yet
           i = 1
           x = 1.0
        endif else begin
           if (time(iTime) gt max(ChampTime)) then begin
              i = iline
              x = 0.0
           endif else begin
              i = loc(0)
              if (i eq 0) then i=1
              x = (champtime(i) - time(itime)) / (champtime(i)-champtime(i-1))
              champdt = champtime(i)-champtime(i-1)
           endelse
        endelse

        ChampDensity(iTimeSave)   = $
           ((1-x) * MassDensity(i) + x*MassDensity(i-1))/1.e-12
        ChampAltitude(iTimeSave)  = $
           (1-x) * ChampPosition(2,i) + x * ChampPosition(2,i-1)
        ChampLats(iTimeSave)      = $
           (1-x) * ChampPosition(1,i) + x * ChampPosition(1,i-1)
        ChampLTs(iTimeSave)       = $
           (1-x) * ChampLocalTime(i) + x * ChampLocalTime(i-1)
        
        loc = where(GitmAlts gt ChampAltitude(iTimeSave),c)
        if (c gt 0) then begin

           i = loc(0)
           x = (ChampAltitude(iTimeSave) - GitmAlts(i-1)) / $
               (GitmAlts(i) - GitmAlts(i-1))

           GitmDensity(iTimeSave) = $
              exp((1.0 - x) * alog(GitmRho(iTime,i-1)) + $
                  (      x) * alog(GitmRho(iTIme,i)))

           GitmeDensity(iTimeSave) = $
              exp((1.0 - x) * alog(GitmElectron(iTime,i-1))+$
                  (      x) * alog(GitmElectron(iTIme,i)))

           ed = reform(GitmElectron(iTime,*))

           l = where(ed eq max(ed))
           GitmNMF2(iTimeSave) = ed(l(0))
           GitmHMF2(iTimeSave) = GitmAlts(l(0))

           l = where(ed gt 0.9*max(ed),c_ed)
           if (c_ed gt 0) then begin
              GitmHMF2l(iTimeSave) = GitmAlts(l(0))
              GitmHMF2u(iTimeSave) = GitmAlts(l(c_ed-1))
           endif else begin 
              GitmHMF2l(iTimeSave) = GitmHMF2(iTimeSave)
              GitmHMF2u(iTimeSave) = GitmHMF2(iTimeSave)
           endelse

           if plotdiff then  $
              GitmDensityBase(iTimeSave) = $
              exp((1.0 - x) * alog(GitmRhoBase(iTime,i-1)) + $
                  (      x) * alog(GitmRho(iTIme,i)))

           GitmHeight(iTimeSave) = ChampAltitude(iTime)

           h = (GitmAlts(i+1) - GitmAlts(i-1))/2.0

           loc = where(GitmAlts gt ChampAltitude(iTimeSave)+h)
           i = loc(0)
           x = ((ChampAltitude(iTimeSave)+h) - GitmAlts(i-1)) / $
               (GitmAlts(i) - GitmAlts(i-1))
           GitmDensityHigh(iTimeSave) = (1.0 - x) * GitmRho(iTime,i-1) + $
                                        (      x) * GitmRho(iTIme,i)
           if plotdiff then $
              GitmDensityHighBase(iTimeSave) = $
              (1.0 - x) * GitmRhoBase(iTime,i-1) + $
              (      x) * GitmRhoBase(iTIme,i)

           loc = where(GitmAlts gt ChampAltitude(iTimeSave)-h)
           i = loc(0)
           x = ((ChampAltitude(iTimeSave)-h) - GitmAlts(i-1)) / $
               (GitmAlts(i) - GitmAlts(i-1))
           GitmDensityLow(iTimeSave) = (1.0 - x) * GitmRho(iTime,i-1) + $
                                       (      x) * GitmRho(iTIme,i)
           if plotdiff then $
              GitmDensityLowBase(iTimeSave) = $
              (1.0 - x) * GitmRhoBase(iTime,i-1) + $
              (      x) * GitmRhoBase(iTIme,i)

        endif else begin
           GitmDensity(iTimeSave)  = GitmDensity(iTimeSave-1)
           GitmeDensity(iTimeSave) = GitmeDensity(iTimeSave-1)
           GitmNMF2(iTimeSave) = GitmNMF2(iTimeSave-1) 
           GitmHMF2(iTimeSave) = GitmHMF2(iTimeSave-1)
           GitmHMF2l(iTimeSave) = GitmHMF2l(iTimeSave-1)
           GitmHMF2u(iTimeSave) = GitmHMF2u(iTimeSave-1)
           GitmDensityHigh(iTimeSave) = GitmDensityHigh(iTimeSave-1)
           GitmDensityLow(iTimeSave) = GitmDensityLow(iTimeSave-1)
        endelse
        
        if (champdt lt 120.0) then begin
           time(iTimeSave) = time(iTime)
           iTimeSave = iTimeSave + 1L
        endif

    endfor
    
    nTimes = iTimeSave-1L

    GitmDensity = GitmDensity(0:nTimes-1)*1.0e12
    GitmDensityHigh = GitmDensityHigh(0:nTimes-1)*1.0e12
    GitmDensityLow = GitmDensityLow(0:nTimes-1)*1.0e12
    GitmeDensity = GitmeDensity(0:nTimes-1)
    GitmNMF2 = GitmNMF2(0:nTimes-1)
    GitmHMF2 = GitmHMF2(0:nTimes-1)
    GitmHMF2u = GitmHMF2u(0:nTimes-1)
    GitmHMF2l = GitmHMF2l(0:nTimes-1)

    ChampDensity = ChampDensity(0:nTimes-1)
    ChampAltitude = ChampAltitude(0:nTimes-1)
    ChampLats = ChampLats(0:nTimes-1)
    ChampLTs = ChampLTs(0:nTimes-1)

    ChampLocalTime = ChampLocalTime(0:nTimes-1)

    if plotdiff then GitmDensityBase = GitmDensityBase*1.0e12
    
    GitmAve  = fltarr(nTimes)
    if plotdiff then GitmAveBase  = fltarr(nTimes)
    ChampAve = fltarr(nTimes)

    ; Here we want to average the data to the champ orbit:
    AltMean = mean(ChampAltitude)

    Re = 6378.0D
    r = (Re + AltMean) * 1000.0D
    G = 6.67259d-11
    Me = 5.9722d24              ; reference 1
    mu = G * Me
    v0 = sqrt(mu/r)
    period = 2*!pi*r / v0 / 60.0
    print, period

    for iTime = 0L, nTimes-1 do begin

        loc = where(abs(time-time(iTime)) lt period*60.0/2, count)
        if (count gt 0) then begin
            GitmAve(iTime) = mean(GitmDensity(loc))
            if plotdiff then GitmAveBase(iTime) = mean(GitmDensityBase(loc))
            ChampAve(iTime) = mean(ChampDensity(loc))
        endif

    endfor

endif


d = [ChampDensity,GITMDensity]
maxi = max(d)
n = n_elements(d)
l = where(d gt maxi,c)
while c lt 0.05*n do begin
   maxi = maxi*0.99
   l = where(d gt maxi,c)
endwhile

maxi = maxi*1.2

yrange = [0,maxi]

;yrange = [0.0,40.0]

ppp = 3
space = 0.04
pos_space, ppp, space, sizes, ny = ppp

get_position, ppp, space, sizes, 0, pos0, /rect
pos0(0) = pos0(0) + 0.05
pos0(2) = pos0(2) - 0.05

get_position, ppp, space, sizes, 1, pos, /rect
pos(0) = pos(0) + 0.05
pos(2) = pos(2) - 0.05

get_position, ppp, space, sizes, 2, pos2, /rect
pos2(0) = pos2(0) + 0.05
pos2(2) = pos2(2) - 0.05

;l = where(time gt 0.0)
;time = time(l)

stime = min(time)
etime = max(time)
time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn

c_r_to_a, itime, stime
c_a_to_ymd, itime, ymd
psfile_raw = 'champ_raw_'+ymd+'.ps'
psfile_ave = 'champ_ave_'+ymd+'.ps'
psfile_2dd = 'champ_2dd_'+ymd+'.ps'
psfile_2dn = 'champ_2dn_'+ymd+'.ps'

;--------------------------------------------------------------------
; Raw
;--------------------------------------------------------------------

if not plotdiff then begin

    setdevice, psfile_raw, 'p', 5
    
    if ((etr-btr)/(24.0*3600.0) lt 1.1) then begin

       plot, time-stime, ChampLats, yrange = [-90,90], pos = pos0, $
             xtickname = xtickname, xtickv = xtickv, $
             xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
             ystyle = 5, $
             thick = 3, /noerase, linestyle = 0

       axis, yaxis=0, ytitle = 'Latitude (deg)', yrange = [-90,90], $
             ystyle = 1, $
             ytickv = findgen(7)*30.0-90, yticks=7, yminor = 3

       plot, time-stime, ChampLTs, yrange = [0,24], pos = pos0, $
             xtickname = xtickname, xtickv = xtickv, $
             xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
             ystyle = 5, $
             thick = 3, /noerase, linestyle = 1

       axis, yaxis = 1, ytitle = 'Local Time (hrs)', yrange = [0,24], $
             ystyle = 1, ytickv = findgen(9)*3, yticks=9, yminor=3

    endif

    plot, time-stime, ChampDensity, yrange = yrange*1.2, pos = pos, $
          xtickname = xtickname, xtickv = xtickv, ystyle = 1, $
          xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
          ytitle = 'Mass Density (10!E-12!N kg/m!E3!N)',   $
          thick = 3, /noerase
    oplot, time-stime, gitmdensity, linestyle = 2, thick = 3
    
    t1 = (etr-btr)*0.05
    t2 = (etr-btr)*0.10
    t3 = (etr-btr)*0.11
    
    oplot, [t1,t2], max(yrange)*[1.05,1.05], thick = 3, linestyle = 2
    xyouts, t3, max(yrange)*1.05, 'GITM'
    
    t1 = (etr-btr)*0.55
    t2 = (etr-btr)*0.60
    t3 = (etr-btr)*0.61
    
    oplot, [t1,t2], max(yrange)*[1.05,1.05], thick = 3
    xyouts, t3, max(yrange)*1.05, 'CHAMP'
    
    d = ChampDensity-gitmdensity 
    yr = [-maxi,maxi]

    plot, time-stime, d, $
      yrange = yr, pos = pos2, $
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
    
    xyouts, pos2(2)+0.01, (pos2(1)+pos2(3))/2.0, srms, $
            alignment = 0.5, orient=270, /norm 

    xyouts, pos2(2)+0.035, (pos2(1)+pos2(3))/2.0, 'CHAMP - GITM', $
            alignment = 0.5, orient=270, /norm 

    oplot, [btr,etr], [0.0,0.0], linestyle = 1
    
    spawn, 'pwd', dir
    dir = dir(0)
    xyouts, 0.0, -0.05, dir, /norm, charsize = 0.8
    
    closedevice
    
;--------------------------------------------------------------------
; Ave
;--------------------------------------------------------------------
    
    setdevice, psfile_ave, 'p', 5
    
    if ((etr-btr)/(24.0*3600.0) lt 1.1) then begin
       plot, time-stime, ChampLats, yrange = [-90,90], pos = pos0, $
             xtickname = xtickname, xtickv = xtickv, $
             xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
             ystyle = 5, $
             thick = 3, /noerase, linestyle = 0

       axis, yaxis=0, ytitle = 'Latitude (deg)', yrange = [-90,90], $
             ystyle = 1, $
             ytickv = findgen(7)*30.0-90, yticks=7, yminor = 3

       plot, time-stime, ChampLTs, yrange = [0,24], pos = pos0, $
             xtickname = xtickname, xtickv = xtickv, $
             xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
             ystyle = 5, $
             thick = 3, /noerase, linestyle = 1

       axis, yaxis = 1, ytitle = 'Local Time (hrs)', yrange = [0,24], $
             ystyle = 1, $
             ytickv = findgen(9)*3, yticks=9, yminor=3

    endif

    plot, time-stime, ChampAve, yrange = yrange*1.2, pos = pos, $
      xtickname = xtickname, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
      ytitle = 'Mass Density (10!E-12!N kg/m!E3!N)',   $
      thick = 3, /noerase
    oplot, time-stime, GitmAve, linestyle = 2, thick = 3
    
    t1 = (etr-btr)*0.05
    t2 = (etr-btr)*0.10
    t3 = (etr-btr)*0.11
    
    oplot, [t1,t2], max(yrange)*[1.05,1.05], thick = 3, linestyle = 2
    xyouts, t3, max(yrange)*1.05, 'GITM'
    
    t1 = (etr-btr)*0.55
    t2 = (etr-btr)*0.60
    t3 = (etr-btr)*0.61
    
    oplot, [t1,t2], max(yrange)*[1.05,1.05], thick = 3
    xyouts, t3, max(yrange)*1.05, 'CHAMP'

    plot, time-stime, ChampAve-gitmave, $
      yrange = yr, pos = pos2, $
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
    
    xyouts, pos2(2)+0.01, (pos2(1)+pos2(3))/2.0, srms, $
            alignment = 0.5, orient=270, /norm 

    xyouts, pos2(2)+0.035, (pos2(1)+pos2(3))/2.0, 'CHAMP - GITM', $
            alignment = 0.5, orient=270, /norm 
    
    oplot, [btr,etr], [0.0,0.0], linestyle = 1
    
    xyouts, 0.0, -0.05, dir, /norm, charsize = 0.8
    
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
   
   for iPt = 0L, nTimes-1 do begin

      c_r_to_a, itime, time(iPt)
      jd = jday(itime(0), itime(1), itime(2))
      printf,1, itime(0), jd, itime(3), itime(4), itime(5), gitmlats(iPt,0), gitmlons(iPt,0), gitmheight(iPt), gitmdensity(iPt)/1.0e12, format = '(i5,4i4,3f7.2,e12.4)'

      printf,2, itime(0), jd, itime(3), itime(4), itime(5), gitmlats(iPt,0), gitmlons(iPt,0), gitmheight(iPt), gitmedensity(iPt)/1.0e6, format = '(i5,4i4,3f7.2,e12.4)'

;      printf,3, itime(0), jd, itime(3), itime(4), itime(5), gitmlats(iPt,0), gitmlons(iPt,0), gitmnmf2(iPt)/1.0e6, format = '(i5,4i4,2f7.2,e12.4)'

;      printf,4, itime(0), jd, itime(3), itime(4), itime(5), gitmlats(iPt,0), gitmlons(iPt,0), gitmhmf2(iPt), gitmhmf2u(iPt), gitmhmf2l(iPt), format = '(i5,4i4,2f7.2,3f12.2)'

   endfor

   close,1,2,3,4

endif

end

