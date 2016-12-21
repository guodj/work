

if (n_elements(basedir) eq 0) then basedir = '.'
basedir = ask('basedir for log file (e.g., ensemble*/data)',basedir)

logfileguess = findfile(basedir+'/log0*.dat')

display, logfileguess

nDirs = n_elements(logfileguess)
;nDirs = 4

for i=0,nDirs-1 do begin

   l =  strpos(logfileguess(i),'log0')
   print, 'reading : ',strmid(logfileguess(i),0,l)
   filelist = findfile(strmid(logfileguess(i),0,l)+'log*.dat')
   thermo_readlog, filelist, time, data, vars

   t  = time-time(0)
   if (i eq 0) then begin

      stime = time(0)
      c_r_to_a, itime, stime
      c_a_to_s, itime, strtime

      nMinutes = long(max(t)/60.0)
      cTime = 'Simulation Time'

      alldata = fltarr(3, nDirs, nMinutes)

   endif

   da = reform(data(10,*))
   db = reform(data( 9,*))
   dc = reform(data(11,*))

   ; linearly interpolate to 1-minute data

   print, 'Interpolating...'
   t1 = dindgen(nMinutes)*60.0
   d2 = da
   alldata(0,i,*) = interpolate_mine(t1, d2, t)

   d2 = db
   alldata(1,i,*) = interpolate_mine(t1, d2, t)

   d2 = dc
   alldata(2,i,*) = interpolate_mine(t1, d2, t)

endfor

setdevice, 'log_compare.ps','p',5

if (max(t) gt 7200.0) then begin
    t  = t1/3600.0
    t2 = t1/3600.0
    cTime = cTime+' (Hours)'
    nDays = max(t)/(24.0)
endif else begin
   t  = t1
   t2 = t1
   cTime = cTime+' (Seconds)'
   nDays = 0
endelse

maxi = max(alldata)*1.1

plot, t, alldata(0,0,*), linestyle = 0, ytitle = 'Temperature (K)', $
  xtitle = cTime, xstyle = 1, pos = [0.1,0.5,0.9,0.9], yrange = [0,maxi]
oplot, t, alldata(1,0,*), linestyle = 1
oplot, t, alldata(2,0,*), thick = 3, linestyle = 2

for i=1,nDirs-1 do begin

   oplot, t, alldata(0,i,*), linestyle = 0
   oplot, t, alldata(1,i,*), linestyle = 1
   oplot, t, alldata(2,i,*), linestyle = 2, thick = 3

endfor

xyouts, 0.9, 0.905, /norm, 'Start Time : '+strtime, alignment = 1

if (nDays gt 1) then begin
   for i=1, nDays do begin
      oplot, [i,i]*24.0, [0, max(data(10,*)*10.0)], linestyle = 1
   endfor
endif

plot, t, alldata(0,1,*)-alldata(0,0,*), $
      linestyle = 0, ytitle = 'Temperature (K)', $
      xtitle = cTime, xstyle = 1, pos = [0.1,0.05,0.9,0.45], /noerase
oplot, t, alldata(1,1,*)-alldata(1,0,*), linestyle = 1
oplot, t, alldata(2,1,*)-alldata(2,0,*), thick = 3, linestyle = 2

closedevice

setdevice, 'log_compare_pd.ps','p',5

maxi = max(alldata)*1.1

plot, t, alldata(0,0,*), linestyle = 0, ytitle = 'Temperature (K)', $
  xtitle = cTime, xstyle = 1, pos = [0.1,0.5,0.9,0.9], yrange = [0,maxi]
oplot, t, alldata(1,0,*), linestyle = 1
oplot, t, alldata(2,0,*), thick = 3, linestyle = 2

for i=1,nDirs-1 do begin

   oplot, t, alldata(0,i,*), linestyle = 0
   oplot, t, alldata(1,i,*), linestyle = 1
   oplot, t, alldata(2,i,*), linestyle = 2, thick = 3

endfor

xyouts, 0.9, 0.905, /norm, 'Start Time : '+strtime, alignment = 1

if (nDays gt 1) then begin
   for i=1, nDays do begin
      oplot, [i,i]*24.0, [0, max(data(10,*)*10.0)], linestyle = 1
   endfor
endif

pd = 100.0*(alldata(0,1,*)-alldata(0,0,*))/alldata(0,0,*)

plot, t, pd, $
      linestyle = 0, ytitle = 'Percentage Difference', $
      xtitle = cTime, xstyle = 1, pos = [0.1,0.05,0.9,0.45], /noerase, $
      yrange = [-15.0,15.0], ystyle = 1
pd =100.0* (alldata(1,1,*)-alldata(1,0,*))/alldata(1,0,*)
oplot, t, pd, linestyle = 1
pd = 100.0*(alldata(2,1,*)-alldata(2,0,*))/alldata(2,0,*)
oplot, t, pd, thick = 3, linestyle = 2

closedevice



end


