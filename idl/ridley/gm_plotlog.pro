ccmc =  1

f = findfile('.gm_plotlog_history')
if (strlen(f) gt 0) then begin
   openr,1, '.gm_plotlog_history'
   syy = ''
   smm = '' 
   sdd = '' 
   shh = '' 
   readf,1,syy
   readf,1,smm
   readf,1,sdd
   readf,1,shh
   close,1
endif

if n_elements(syy) eq 0 then syy = '1998'
syy = ask('year',syy)

if n_elements(smm) eq 0 then smm = '05'
smm = chopr('0'+ask('month',smm),2)

if n_elements(sdd) eq 0 then sdd = '04'
sdd = chopr('0'+ask('day',sdd),2)

if n_elements(shh) eq 0 then shh = '02'
shh = chopr('0'+ask('hour',shh),2)

if (strlen(f) eq 0) then begin
   openw,1, '.gm_plotlog_history'
   printf,1,syy
   printf,1,smm
   printf,1,sdd
   printf,1,shh
   close,1
endif

itime = [fix(syy),fix(smm),fix(sdd),fix(shh),0,0]
c_a_to_r, itime, starttime
c_a_to_ymd, itime, ymd

filelist = findfile('log_n*.log')

nfiles=n_elements(filelist)
for iFile = 0L, nFiles-1 do begin
    logfilename = filelist(iFile)
    getlogpro, logfilename, $
      nlogfile, nVars, wlognames, wlog, wlog0, wlog1, wlog2

    if (iFile eq 0) then begin
        wlogtotal = wlog
    endif else wlogtotal = [wlogtotal, wlog]
endfor

wlog = wlogtotal

time = starttime + reform(double(wlog(*,1)))

stime = min(time)
etime = max(time)

download_dst, itime, dstfile
read_dst, dstfile, dst, dsttime

; Need to see whether you need another Dst file!

c_r_to_a, itime_start, stime
c_r_to_a, itime_end, etime

if (itime_start(1) ne itime_end(1)) then begin
   download_dst, itime_end, dstfile2
   read_dst, dstfile2, dst2, dsttime2
   dst = [dst,dst2]
   dsttime = [dsttime,dsttime2]
endif

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

iVar1 = 13
iVar2 = 15
mhddst1 = wlog(*,iVar1)
mhddst2 = wlog(*,iVar2)

l = where(dsttime ge stime-3600.0 and dsttime le etime+3600.0)

mini = min([dst(l),mhddst1,mhddst2])
maxi = max([dst(l),mhddst1,mhddst2])
r = maxi-mini

yrange = [mini - r/10, maxi + r/10]

setdevice, 'dst.ps', 'p', 5
plot, time-stime, mhddst2, $
  ytitle = 'dst (nT)', charsize = 1.2, thick = 2, pos = [0.1, 0.3, 0.9, 0.7],$
  xstyle = 1, xrange = [btr,etr], $
  xtickname = xtickname, xminor=xminor, xticks=xtickn, $
  xtickv = xtickv, xtitle = xtitle, yrange = yrange, ystyle = 1
;oplot, time-stime, dst1, linestyle = 1

oplot, dsttime - stime, dst, linestyle = 2, psym = -4

mhddsti = dst*0.0-1.0e32
dt = dsttime(1) - dsttime(0)
nDst = n_elements(dst)
for i = 0L, nDst-1 do begin
    l = where(time ge dsttime(i)-dt/2 and time lt dsttime(i)+dt/2,c)
    if c gt 0 then mhddsti(i) = mean(mhddst2(l))
endfor

l = where(abs(dsttime-time(0)) eq min(abs(dsttime-time(0))))
offset = 0.0 ; dst(l(0)) - mhddsti(l(0))
oplot, dsttime - stime, mhddsti + offset, linestyle = 1, psym = -5

l = where(mhddsti gt -10000.0,c)
;l = where(mhddsti gt -10000.0 and dsttime-stime le 25200.0,c)
;print, "ONLY USING FIRST 7 HOURS!!!"
if (c gt 0) then rms = sqrt(mean((mhddsti(l)-dst(l))^2))/sqrt(mean((dst(l))^2))

xyouts, 0.89, 0.67, 'nRMS : '+string(rms, format = '(f6.3)'), /norm, $
  alignment = 1.0

closedevice

setdevice, 'log.ps', 'p', 5

iVarDt = 2
mhddt = wlog(*,iVarDt)

nPts = n_elements(wlog(*,0))
nVars = 12
data = fltarr(nVars,nPts)
for i=0,nVars-1 do data(i,*) = reform(wlog(*,i+2))
; pressure
data(nVars-3) = alog10(data(nVars-3))
data(nVars-2) = alog10(data(nVars-2))
overplots = [0,0,0,1,1,0,0,1,1,0,1,0]
titles = wlognames(2)
for i=1,nVars-1 do if overplots(i) eq 0 then titles = [titles,wlognames(i+2)]

plot_data, data, time, nvars, titles, overplots

closedevice

if ccmc then begin

   close,1
   openw,1,'dst'+ymd+'.ccmc'

   nPts = n_elements(mhddst2)
   for i=0L,nPts-1 do begin

      if (time(i) gt time(0)) then begin
         c_r_to_a, itime, time(i)
         printf,1, itime, mhddst2(i), format = '(i5,5i3,f8.1)'
      endif

   endfor

   close,1

endif

end


