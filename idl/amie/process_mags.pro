
common ffinfo, header

spawn, 'date -u +%Y', sYear
spawn, 'date -u +%m', sMonth
spawn, 'date -u +%d', sDay

Print, "This code will take the year, month, and day as input,"
print, "and will process all of the data 28-31 days before this time,"
print, "if it can find that much data."
print, "So, enter the LAST day to process!!"

sYear  = ask('year to process',sYear)
sMonth = ask('month to process',sMonth)
sDay   = ask('day to process',sDay)

Year  = fix(sYear)
Month = fix(sMonth)
Day   = fix(sDay)
itime = [Year, Month, Day, 12, 00, 00]

if (month eq  1) then nDays = 31
if (month eq  2) then begin
    if (Day eq 29) then nDays = 29 else nDays = 28
endif
if (month eq  3) then nDays = 31
if (month eq  4) then nDays = 30
if (month eq  5) then nDays = 31
if (month eq  6) then nDays = 30
if (month eq  7) then nDays = 31
if (month eq  8) then nDays = 31
if (month eq  9) then nDays = 30
if (month eq 10) then nDays = 31
if (month eq 11) then nDays = 30
if (month eq 12) then nDays = 31

nDays = fix(ask('number of days to process',tostr(nDays)))

if (nDays lt 5) then begin
    quietdayfiles = findfile('*quietday.save')
    nQuietDayFiles = n_elements(quietdayfiles)
    if (nQuietDayFiles gt 1) then begin
        print, "You have asked for a very small amount of days,"
        print, "and you have some quietday files.  Would you like to"
        print, "use those files for the quiet day subtraction?"
        UseQuietDayFiles = mklower(ask('(y)es or (n)o','y'))
    endif else UseQuietDayFiles = 'n'
endif else UseQuietDayFiles = 'n'

OutputRemca = ask('whether you want to output remca files or not (y or n)','y')

c_a_to_r, itime, LastTime
StartTime = LastTime - (nDays-1)*24.0*3600.0

for iDay = 0, nDays-1 do begin

    CurrentTime = StartTime + iDay * 24.0*3600.0
    c_r_to_a, itime, CurrentTime

    sYear  = tostr(itime(0))
    sMonth = chopr('0'+tostr(itime(1)),2)
    sDay   = chopr('0'+tostr(itime(2)),2)

    filelist_tmp = findfile(sYear+'/'+sMonth+'/'+sDay+'/*.hed')

    if (iDay eq 0) then filelist = filelist_tmp $
    else filelist = [filelist,filelist_tmp]

endfor

nfiles = n_elements(filelist)
filelist_master = filelist

statlist = [strmid(filelist(0),11,3)]
nstats = 1

for i=1,nfiles-1 do begin

  stat = strmid(filelist(i),11,3)

  print, stat

  ifound = 0
  istat = 0

  while (not ifound and istat lt nstats) do begin
    if (stat eq statlist(istat)) then ifound = 1
    istat = istat + 1
  endwhile

  if (not ifound) then begin
    statlist = [statlist,stat]
    nstats = nstats+1
  endif

endfor

; We need to get the latitude of the stations - so we have to read in the
; master.psi file

openr,1,'master.psi'
ns = 0
line = ''
readf,1,ns
readf,1,line
stations = strarr(ns)
lats     = fltarr(ns)
for i=0,ns-1 do begin
  readf,1,line
  stations(i) = mklower(strmid(line,4,3))
  lats(i)     = float(strmid(line,36,6))
endfor
close,1

maglats = fltarr(nstats)

for istat = 0,nstats-1 do begin

    print, "station : ",statlist(iStat), iStat, nStats-1

    loc = where(stations eq statlist(istat),count)

    if count gt 0 then latitude = lats(loc(0)) else begin
        print, 'working on station : ',statlist(istat),istat
        print, "Station not found!!!"
        latitude = 0.0
    endelse

;  print, "Station latitude : ", latitude
    maglats(iStat) = latitude

endfor

ind = sort(maglats)

itime = intarr(6)

nPtsMax = 1440

nStatsToProcess = nStats
nStatsToProcess = fix(ask('Number of stations to process',tostr(nStatsToProcess)))

if (nStatsToProcess eq 1) then begin 

    for istat = nstats-1, 0, -1 do begin

        iMag = ind(istat)
        lat = maglats(iMag)

        print, istat,'. ',statlist(ind(istat))

    endfor

    iStatStart = fix(ask('stat to process','0'))
    iStatEnd = iStatStart

endif else begin

;    fl = findfile('????01/*.remca')
;    iStatStart = nStats-n_elements(fl)
    iStatStart = nStats-1
    iStatEnd = 0

endelse

IsFirst = 1

for istat = iStatStart, iStatEnd, -1 do begin

    iMag = ind(istat)
    lat = maglats(iMag)

    n = 0
    filelist = strarr(31)

print, nFiles

    for i=0,nfiles-1 do begin
        if (strpos(filelist_master(i),statlist(ind(istat))) ge 0) then begin

print, filelist_master(i),i,n
            filelist(n) = filelist_master(i)
            n = n + 1
        endif
    endfor

    filelist = filelist(0:n-1)

    ndays = n_elements(filelist)

    print, "Station : ",statlist(ind(istat)), ndays

    for n=0,ndays-1 do begin
  
        itime(0) = fix(strmid(filelist(n),14,2))
        itime(1) = fix(strmid(filelist(n),16,2))
        itime(2) = fix(strmid(filelist(n),18,2))
        itime(3) = 0
        itime(4) = 0
        itime(5) = 0

        c_a_to_r, itime, stime

        itime(2) = itime(2) + 1
        c_a_to_r, itime, etime

        col_scal = [0,1,2]

        filename = strmid(filelist(n),0,20)
        print, "filename : ", filename

        spawn, 'wc '+filename+'.hed', wc
        if (fix(wc) gt 10) then begin
            read_flat_scalor, stime, etime, col_scal, time1f, data1f, nrows,  $
              filename = filename
        endif else begin
            nRows = [1,1,1]
            data1f = fltarr(3,1)
            time1f = dblarr(3,1)
        endelse

        if (n eq 0) then begin
            data    = fltarr(3,nRows(0),ndays) - 99999.0
            time    = dblarr(nRows(0),ndays)
            dataqdr = fltarr(3,nRows(0),ndays) - 99999.0

            if (fix(wc) gt 10) then begin
                readhed, 1, filename, timeint, variables, units, 	$
                  hedcount, hednrows, hedrowlen, vartype, heddata_struct
                components = $
                  strmid(variables(0,0),3,1)+$
                  strmid(variables(0,1),3,1)+$
                  strmid(variables(0,2),3,1)
            endif else components = 'XYZ'

        endif

        if (n_elements(data(0,*,n)) eq n_elements(data1f(0,*))) then begin
            data(*,*,n) = data1f
            time(*,n)   = reform(time1f(0,*))
        endif else begin

            npts_new = n_elements(data1f(0,*))
            npts_old = n_elements(data(0,*,n))

            if (npts_new gt npts_old) then begin

                datanew    = fltarr(3,npts_new,ndays)
                timenew    = dblarr(npts_new,ndays)
                dataqdrnew = fltarr(3,npts_new,ndays) - 99999.0

                print, "Growing Array (Including More Data)"

                for n2 = 0,n-1 do begin
                    for iComp = 0,2 do begin
                        datanew(iComp,0:npts_old-1,n2) = $
                          data(iComp,0:npts_old-1,n2)
                        datanew(iComp,npts_old:npts_new-1,n2) = -1.0e32
                    endfor
                    timenew(0:npts_old-1,n2) = time(0:npts_old-1,n2)
                    timenew(npts_old:npts_new-1,n2) = $
                      time(npts_old-1,n2) + $
                      dindgen(npts_new-npts_old)/(npts_new-npts_old)
                endfor
                data = datanew
                time = timenew

                data(*,*,n) = data1f
                time(*,n)   = reform(time1f(0,*))

            endif else begin

                print, "Filling Array (Including Less Data)"

                for iComp = 0,2 do begin
                    data(iComp,0:npts_new-1,n) = $
                      data1f(iComp,0:npts_new-1)
                    data(iComp,npts_new:npts_old-1,n) = -1.0e32
                endfor
                time(0:npts_new-1,n) = reform(time1f(0,0:npts_new-1))
                time(npts_new:npts_old-1,n) = $
                  time(npts_new-1,n) + $
                  dindgen(npts_old-npts_new)/(npts_old-npts_new)

            endelse

        endelse

    endfor

    loc = where(abs(data) lt 70000,count)

    ; Choose 100, since it is about 31 days * 3 components

    if (count gt 100) then begin

        despike_mags, data, time, lat, quality
        average_mags, data, time, lat, quality, nPtsMax, $
          averagedata, averagetime, averagequality

        data = averagedata
        time = averagetime
        quality = averagequality

        if (strpos(UseQuietDayFiles,'n') gt -1) then begin
            quiet_day, data, time, lat, quality, quietdata, QuietDay
            fix_mags, quietdata, time, LastDayOffset
            for iComp=0,2 do $
              QuietDay(iComp,*) = QuietDay(iComp,*) + LastDayOffset(iComp)
            save, QuietDay, file=statlist(ind(istat))+'.quietday.save'
        endif else begin
            QuietDayFile = findfile(statlist(ind(istat))+'.quietday.save')
            if (strlen(QuietDayFile) gt 0) then begin
                print, "Restoring File : ",QuietDayFile
                restore, QuietDayFile
                quiet_day, data, time, lat, quality, quietdata, $
                  QuietDay, /usequietday
            endif else begin
                quietdata = data*0.0 - 99999.0
            endelse
        endelse

        loc = where(abs(quietdata) lt 70000,count)

    endif

    if (count gt 100) then begin

;        if (iStat eq iStatStart) then begin

        if (IsFirst) then begin
                IsFirst = 0

            if (strpos(UseQuietDayFiles,'n') gt -1) then $
              psfile = 'stats_all.ps' $
            else psfile = 'stats_rt.ps'
            if !d.name eq 'X' then setdevice,psfile,'l',4,0.95
  
            ppp = 3
            space = 0.01
            pos_space, ppp, space, sizes, ny = ppp
  
            get_position,ppp,space,sizes,0,pos1,/rect
            get_position,ppp,space,sizes,1,pos2,/rect
            get_position,ppp,space,sizes,2,pos3,/rect
            pos1(0) = pos1(0) + 0.05
            pos2(0) = pos2(0) + 0.05
            pos3(0) = pos3(0) + 0.05
  
            pos = fltarr(3,4)
            pos(0,*) = pos1
            pos(1,*) = pos2
            pos(2,*) = pos3
  
;            stime_keep = mean(time)
            stime_keep = StartTime
            c_r_to_a, itime, stime_keep
;            itime(2) = 1
            itime(3) = 0
            itime(4) = 0
            itime(5) = 0
            c_a_to_r, itime, stime_keep

            c_r_to_a, itime, LastTime
            itime(3) = 23
            itime(4) = 59
            itime(5) = 59
            c_a_to_r, itime, etime_keep

;            if (itime(1) eq 1) then itime(2) = 31
;            if (itime(1) eq 2) then itime(2) = 29
;            if (itime(1) eq 3) then itime(2) = 31
;            if (itime(1) eq 4) then itime(2) = 30
;            if (itime(1) eq 5) then itime(2) = 31
;            if (itime(1) eq 6) then itime(2) = 31
;            if (itime(1) eq 7) then itime(2) = 30
;            if (itime(1) eq 8) then itime(2) = 31
;            if (itime(1) eq 9) then itime(2) = 30
;            if (itime(1) eq 10) then itime(2) = 31
;            if (itime(1) eq 11) then itime(2) = 30
;            if (itime(1) eq 12) then itime(2) = 31
;            itime(3) = 23
;            itime(4) = 59
;            itime(5) = 59
;            c_a_to_r, itime, etime_keep
  
            time_axis, stime_keep, etime_keep, s_time_range, e_time_range, $
              xtickname, xtitle, xtickvalue, xminor, xtickn
  
;            comp = ['X','Y','Z']
  
        endif

        i = ind(istat)

        plotdumb

        for icomp = 0,2 do begin

            if (icomp eq 2) then begin
                xtn = xtickname
                xt  = xtitle
            endif else begin
                xtn = strarr(20)+' '
                xt  = ' '
            endelse

            title = ' '
            if (icomp eq 0) then $
              title =statlist(i)+' '+tostr(i+1)+'/'+ $
              tostr(istat)+' of '+tostr(nstats)+' mlat : '+tostr(fix(lat))

            plot, time-stime_keep,quietdata(icomp,*,*), pos = pos(icomp,*),$
              xstyle = 1, /noerase,		$
              xtickname = xtn, xtickv=xtickvalue, 			$
              xticks = xtickn, xminor = xminor, xtitle = xt,		$
              xrange = [s_time_range, e_time_range],                  $
              ytitle = strmid(components,icomp,1), $
              yrange = yrange, min_val = -5000.0, $
              title = title

            for iDay = 1, 31 do begin
                oplot, 3600.0*24*[iDay, iDay],[-5000,5000], linestyle = 1
            endfor

            oplot, [s_time_range, e_time_range], [0.0,0.0], linestyle = 1

            d = reform(quietdata(icomp,*,*))
            loc = where(abs(d) lt 2500.0,count)
            
            if (count gt 0) then begin
                m = median(d(loc))
                mi = min(d(loc))
                hi = histogram(d(loc))
                ma = where(hi eq max(hi))
                m  = ma(0) + mi
            endif else m = -99999.

            ms = strcompress(string(m))

            xyouts, pos(icomp,2)+0.02, pos(icomp,1), ms, orient = 90, /norm

        endfor

        if (strpos(OutputRemca,'y') gt -1) then begin

            loc = where(quietdata lt -1.0e6,count)
            if count gt 0 then quietdata(loc) = -99999.0

            for n1 = 0, ndays-1 do begin
                c_r_to_a, itime, mean(time(*,n1))
                ymd = chopr('0'+tostr(itime(0)),2)+ $
                  chopr('0'+tostr(itime(1)),2)+ $
                  chopr('0'+tostr(itime(2)),2)

                sYear  = tostr(itime(0))
                sMonth = chopr('0'+tostr(itime(1)),2)
                
                filename = sYear+'/'+sMonth+'/'+ymd+'/'+$
                  mklower(statlist(i))+ymd+'.remca'
                print, "Writing File : ", filename
                openw,1,filename, error = iError
            
                if (iError ne 0) then begin
                    spawn, "mkdir "+sYear+'/'+sMonth+'/'+ymd
                    openw,1,filename
                endif
            
                for it=0,nptsmax-1 do begin
                    c_r_to_a, itime, time(it,n1)
                    if (itime(5) ge 30) then itime(4) = itime(4) + 1
                    printf,1,format='(5i2,a,1x,a,3f10.1)', $
                      itime(0) mod 100,itime(1),itime(2),itime(3),itime(4),$
                      statlist(i), components, $
                      quietdata(0,it,n1),quietdata(1,it,n1),quietdata(2,it,n1)
                endfor
                close,1
            endfor

        endif

    endif

endfor

closedevice

end
