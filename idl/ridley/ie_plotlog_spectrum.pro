
pro reform_data, data, time, dt

  dttotal = max(time)-min(time)
  nTotal = long(dttotal/dt)+1L

  datanew = fltarr(nTotal)
  timenew = min(time)+findgen(nTotal)*dt

  datanew(0) = data(0)
  for i=1,nTotal-2 do begin
      l = where(time gt timenew(i))
      l = l(0)
      if (l le 0) then l = 1
      r = (time(l)-time(l-1))/dt
      datanew(i) = (1-r)*data(l) + r*data(l-1)
  endfor

  data = datanew
  time = timenew

end

pro fill_spectrum, data, time, spectrum, x, y

    DtStep   = 60.0
    DtWindow = 60.0 * 60.0
    Dt       = time(1) - time(0)

    nTimesTotal = max(time) / DtStep + 1
    nPts        = DtWindow/Dt + 1

    spectrum = fltarr(nTimesTotal, nPts/2) + 1.0e32
    x        = fltarr(nTimesTotal, nPts/2)
    y        = fltarr(nTimesTotal, nPts/2)

    df = 0.5/(dt*float(npts))
    freq = findgen(nPts/2) * df
    c_window, 2*nPts, wind, 1

    for i = 0, nTimesTotal-1 do begin

        deltat = abs(time - DtStep*i)
        l = where(deltat le DtWindow, c)
        x(i,*) = DtStep*i
        y(i,*) = freq*1000.0

        if (c gt nPts-2) then begin
            signal = reform(data(l)) - mean(data(l))
            dumamp   = fft(signal*wind,-1)
            power = (2.0*abs(dumamp(0:nPts/2)))^2 / df
            spectrum(i,0:nPts/2-1) = power(0:nPts/2-1)
        endif

    endfor

end

filelist = findfile('IE_??.log')

itime = intarr(6)
nSolve = 0
t = 0.0
imilli = 0
tilt1 = 0.0
tilt2 = 0.0
cpcptemp1 = 0.0
cpcptemp2 = 0.0

nMaxPoints = 30000

line = ''

iLine = 0
time = dblarr(nMaxPoints)
cpcp = fltarr(nMaxPoints,2)
iter = intarr(nMaxPoints)

close,1

for i=0, n_elements(filelist)-1 do begin

    openr,1,filelist(i)
    readf,1,line
    readf,1,line

    while not eof(1) do begin
        readf,1,nSolve,t,itime,imilli,tilt1,tilt2,cpcptemp1,cpcptemp2
        iter(iLine) = nSolve
        c_a_to_r, itime, rtime
        time(iLine) = rtime
        cpcp(iLine,0) = cpcptemp1
        cpcp(iLine,1) = cpcptemp2
        iLine = iLine + 1
    endwhile
    close,1

endfor

l1 = where(time eq time(0), nSteady)
if (nSteady gt 0) then begin
    l = where(time gt time(0), nLines)
    l = [l(0)-1,l]
    nLines = nLines + 1
endif else begin
    l = where(time ge time(0), nLines)
endelse

iter = iter(l)
time = time(l)
cpcp = cpcp(l,*)

imffile = findfile('../imf*.dat')
read_imf_input, imffile(0), imftime, mag, vel, den, temp, nPts

stime = min(time)
etime = max(time)

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

setdevice, 'cpcp.ps','p', 5

plotdumb

yr = mm(cpcp(*,1))
r = yr(1)-yr(0)
yr(0) = yr(0)-r*0.1
yr(1) = yr(1)+r*0.1
plot, time-stime, cpcp(*,1), xstyle = 1, xrange = [btr,etr], $
  xtickname = xtickname, xminor=xminor, xticks=xtickn, $
  xtickv = xtickv, ystyle = 1, pos = [0.1,0.03, 0.9,0.3], $
  xtitle = xtitle, ytitle = 'CPCP (kV)', yrange = yr, /noerase

MeanCpcp = mean(cpcp(*,1))

iMag = 2
std = stddev(mag(2,*))
if (std gt 0) then begin
    ytitle = 'IMF Bz (nT)'
endif else begin
    ytitle = 'IMF By (nT)'
    iMag = 1
endelse

if (stddev(mag(iMag,*)) gt 0) then begin
    MeanImf  = mean(mag(iMag,*))
    Var = reform(mag(iMag,*))
    Delay = 15.0*60.0
endif else begin

    if (stddev(den) gt 0) then begin
        MeanSW  = mean(den)
        ytitle = 'Density (/cm3)'
        Var = den
    endif else begin
        MeanSW  = mean(vel(0,*))
        ytitle = 'Velocity (km/s)'
        Var = reform(vel(0,*))
    endelse
    Delay = 10.0*60.0

endelse

yr = mm(Var)
r = yr(1)-yr(0)
yr(0) = yr(0)-r*0.1
yr(1) = yr(1)+r*0.1
plot, imftime-imftime(0)+Delay, Var, xstyle = 1, xrange = [btr,etr], $
  xtickname = xtickname, xminor=xminor, xticks=xtickn, $
  xtickv = xtickv, ystyle = 1, pos = [0.1,0.33, 0.9,0.5], $
  ytitle = ytitle, yrange = yr, /noerase

makect,'mid'

t = imftime-imftime(0)
v = Var

dt = imftime(1)-imftime(0)
reform_data, v, t, dt

;yfreq = [0.001,2.0,4.0,6.0,8.0,10.0,12.0,14.0]*0.001
;yfreq = [0.00001,1.0,2.0,3.0,4.0,5.0]*0.001
yfreq = [0.00001,0.5,1.0,1.5,2.0,2.5]*0.001
yperiod = 1.0/yfreq/60.0
ytickname = string(yperiod, format='(f4.1)')
ytickname(0) = ' '

fill_spectrum, v, t, imf_spectrum, imf_x, imf_y

l = where(imf_spectrum gt 1.0e31, c)
if (c gt 0) then imf_spectrum(l) = min(imf_spectrum)

nPtsX = n_elements(imf_spectrum(*,0))
nPtsY = n_elements(imf_spectrum(0,*))
maxvals = fltarr(nPtsX)
for i=0,nPtsX-1 do begin
    l = where(imf_spectrum(i,1:nPtsY-1) eq max(imf_spectrum(i,1:nPtsY-1)))
    maxvals(i) = imf_y(i,l(0)+1)
endfor

contour, alog10(imf_spectrum+1.0e-7), imf_x+Delay, imf_y, $
  /fill, nlevels = 31, $
  xstyle = 1, xrange = [btr,etr], yrange =mm(yfreq)/0.001, $
  xtickname = xtickname, xminor=xminor, xticks=xtickn, $
  xtickv = xtickv, ystyle = 1, pos = [0.1,0.78, 0.9,1.00], $
  /noerase, ytitle = 'Frequency (mHz) IMF'

axis, yaxis=1, yrange = mm(yfreq)/0.001, ystyle = 1, $
  ytitle = 'Period (min)', ytickname = ytickname

t = time-stime
v = cpcp(*,1)

fill_spectrum, v, t, ion_spectrum, x, y

l = where(ion_spectrum gt 1.0e31, c)
if (c gt 0) then ion_spectrum(l) = min(imf_spectrum)

contour, alog10(ion_spectrum), x, y, /fill, nlevels = 31, $
  xstyle = 1, xrange = [btr,etr], yrange =mm(yfreq)/0.001, $
  xtickname = xtickname, xminor=xminor, xticks=xtickn, $
  xtickv = xtickv, ystyle = 1, pos = [0.1,0.53, 0.9,0.75], $
  /noerase, ytitle = 'Frequency (mHz) Iono'

oplot, imf_x(*,0)+Delay, maxvals, thick = 3, linestyle = 1

axis, yaxis=1, yrange = mm(yfreq)/0.001, ystyle = 1, $
  ytitle = 'Period (min)', ytickname = ytickname

closedevice


end
