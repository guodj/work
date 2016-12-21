
filelist = findfile('log*.dat')

nFiles = n_elements(filelist)

nLines = 0L

for i=0,nFiles-1 do begin

    spawn, 'wc '+filelist(i), nL
    nLines = nLines+long(nL)-2

endfor

line = ''

iLine = 0L

for i=0,nFiles-1 do begin

    openr, 1, filelist(i)

    while (strpos(line,"START") lt 0) do readf,1,line
    readf,1,line

    if (i eq 0) then begin

        l = strsplit(line)
        nVars = n_elements(l)
        l = [l,strlen(line)]
        Vars = strarr(nVars)
      
        for iVar = 0, nVars-1 do $
           Vars(iVar) = strmid(line,l(iVar),l(iVar+1)-l(iVar))

        data = fltarr(nVars,nLines)
        tmp = fltarr(nVars)

    endif

    while not eof(1) do begin

        readf,1,tmp
        data(*,iLine) = tmp
        iLine = iLine + 1L

    endwhile

    close,1

endfor

setdevice, 'log.ps','p',5

nLines = iLine
time = dblarr(nLines)

for i=0L,nLines-1 do begin

    itime = fix(data(1:6,i))
    c_a_to_r, itime, rtime
    time(i) = rtime

endfor

stime = time(0)
c_r_to_a, itime, stime
c_a_to_s, itime, strtime

t = time-time(0)
cTime = 'Simulation Time'

if (max(t) gt 7200.0) then begin
    t = t/3600.0
    cTime = cTime+' (Hours)'
    nDays = max(t)/(24.0)
endif else begin
    cTime = cTime+' (Seconds)'
    nDays = 0
endelse

maxi = min([4000.0,max(data(10,*))])

dt = reform(data(8,*))

plot, t, data(10,*), linestyle = 2, ytitle = 'Temperature (K)', $
  xtitle = cTime, xstyle = 1, pos = [0.1,0.55,0.9,1.0], yrange = [0,maxi]
oplot, t, data(9,*), linestyle = 2
oplot, t, data(11,*), thick = 3

plot, t, dt, thick = 3, pos = [0.1,0.05,0.9,0.5], $
      xtitle = cTime, xstyle = 1, /noerase, ytitle = 'dT (s)'

xyouts, 0.9, 0.705, /norm, 'Start Time : '+strtime, alignment = 1

if (nDays gt 1) then begin
   for i=1, nDays do begin
      oplot, [i,i]*24.0, [0, max(data(10,*)*10.0)], linestyle = 1
   endfor
endif

closedevice

setdevice, 'gitm_hpi.ps','p',5

stime = min(time)
etime = max(time)
t = time - stime

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

plot, t, data(21,*), linestyle = 0, ytitle = 'Hemispheric Power (GW)', $
      xtickname = xtickname,			$
      xtitle = xtitle,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn,   $
      pos = [0.1,0.3,0.95,0.8], /noerase, $
      yrange = [0,1000], $
      thick = 4

oplot, t, data(20,*), linestyle = 2, color = 180, thick = 3
oplot, t, data(21,*), thick = 4

dx = 0.1*(etr-btr)
oplot, [dx, 2*dx], [900,900], thick = 3
xyouts, dx*2.1, 900, 'GITM Hemispheric Power'
oplot, [dx, 2*dx], [825,825], thick = 3, linestyle = 2, color = 180
xyouts, dx*2.1, 825, 'NOAA Hemispheric Power'


closedevice

end


