
psfile = ask("ps filename", "idl.ps")

setdevice, psfile, "p", 4, 0.95

if (n_elements(logfilename) eq 0) then begin
    filelist = findfile("*.sat")
    if (strlen(filelist(0)) gt 0) then logfilename = filelist(0) $
    else logfilename = "idontknow.sat"
endif 
logfilename = ask("MHD satellite output file",logfilename)

getlogpro, logfilename, nlogfile, nVars, wlognames, wlog, wlog0, wlog1, wlog2


;------------------------------------------------------
; Get delay
;------------------------------------------------------

if (n_elements(delay) eq 0) then delay = " 0.0"
delay = float(ask("delay in seconds",string(delay)))

npts = n_elements(wlog(*,0))

time = dblarr(npts)

if (strpos(wlognames(0),'it') eq 0) then iStart = 1 else iStart = 0

for i=0,npts-1 do begin
    itime = reform(fix(wlog(i,iStart:iStart+5)))
    c_a_to_r, itime, rtime
    time(i) = rtime + delay
endfor

IsSat = 1

if (IsSat) then begin
    x = reform(wlog(*,8))
    y = reform(wlog(*,9))
    z = reform(wlog(*,10))
endif else begin
    nTimes = n_elements(wlog(*,0))
    x = fltarr(nTimes) + 32.0
    y = fltarr(nTimes) + 0.0
    z = fltarr(nTimes) + 0.0
endelse

xmini = floor(min([min(x),-10.0])/10.0)*10.0
xmaxi = round(max([max(x), 10.0])/10.0)*10.0

ymaxi = floor(max([max(abs(y)),10.0])/10.0)*10.0
zmaxi = floor(max([max(abs(z)),10.0])/10.0)*10.0

xrange = [ xmini, xmaxi]
yrange = [-ymaxi, ymaxi]
zrange = [-zmaxi, zmaxi]

stime = min(time)
etime = max(time)

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

tmptime = time - stime

dt = (etime-stime)/float(npts-1)

; Bx By Bz for an IMF file is 7,8,9, rho and ux = 13,10

if (IsSat) then vars = [15,16,17] else vars = [7,8,9] ; vars = [13,10]

nv = n_elements(vars)

  timefft = 0.0
  print, 'Enter time interval for each FFT (in minutes) :'
  read, timefft
  timefft = 60.0*timefft
  ntime = fix(timefft/dt)

  step = 0.0
  print, 'Enter step in between FFTs (in minutes) :'
  read, step
  nstep = fix(step*60.0/dt)
  nfft = (npts-ntime)/nstep + 1

  window_shape = 0
  print,'  '
  print,'   1 - no window (rectangular)'
  print,'   2 - trianglular'
  print,'   3 - blackman'
  print,'   4 - hamming'
  print,'   5 - hanning'
  print,'   6 - kaiser (not available)'
  print,'   7 - lanczos (not available)'
  print,'   8 - tukey (not available)'

  while (window_shape lt 1) or (window_shape gt 5) do begin

    print,'Enter window type :'
    read, window_shape

    if window_shape lt 1 or window_shape gt 5 then 			$
      print, 'Sorry, not available.... try again'

  endwhile  

  c_window, ntime, wind, window_shape

  signal  = fltarr(nv,npts)
  for iVar = 0,nv-1 do signal(iVar,*) = wlog(*,vars(iVar))
  t       = fltarr(npts)
  power   = fltarr(nv,nfft,ntime)
  sigmean = fltarr(nv,nfft)
  freq    = fltarr(ntime)
  period  = fltarr(ntime)
  df      = 1.0/(dt*float(ntime))

    for ii = 0, ntime-1 do begin

      if ii le ntime/2 then freq(ii) = float(ii)*df 			$
      else freq(ii) = float(ntime/2 - ii)*df

    endfor

    period(1:ntime/2) = 1.0/freq(1:ntime/2)
    period(0) = timefft
    period = period(0:ntime/2)

    dumsig = fltarr(ntime)

    for i=0,nv-1 do for k=0,npts-1-ntime,nstep do begin

	kb = k
	ke = k + ntime - 1

	dumsig   = signal(i,kb:ke)
	remove_mean, dumsig, mean
	sigmean(i,k/nstep) = mean
	dumsig = dumsig + mean
	dumamp   = fft(dumsig*wind,-1)
	power(i,k/nstep,*) = (2.0*abs(dumamp))^2.0/df

    endfor

  signal = signal(*,ntime/2:npts-1-ntime/2)

  maxf = 0.0
  print, 'The maximum frequency is :', 1.0/(dt*2.0)
  print, 'Enter the freqency that you would like to cut off at :'
  read, maxf 
  floc = where(freq le maxf and freq ge 0.0, fcount)


  minspec = 0.0
  maxspec = 0.0
  print, ' '
  print, 'The max and min values of the power are : ',mm(power(*,*,floc))
  for i=-1,9 do begin
   ftemploc = where(power(*,*,floc) gt 10^i and 			     $
		    power(*,*,floc) lt 10^(i+1), ftempcount)
   print, fcount*ntime, ftempcount
   print, tostr(fix(100.0*float(ftempcount)/(float(fcount)*float(ntime)))),$
     ' % between 10e',tostr(i),' and 10e', tostr(i+1)
  endfor
  print, 'Enter min,max values for color table (ex 1.0,1.0e8) :'
  read, minspec, maxspec
  minspec = alog10(minspec)
  maxspec = alog10(maxspec)

  print, ' '
  print, 'Npts : ',npts,'  dt : ',dt

nvars = n_elements(vars)
ppp = nvars+2
space = 0.01
pos_space, ppp, space, sizes, ny = ppp

if strpos(!d.name,'X') gt -1 then thick = 1 else thick = 3

plotdumb

get_position, ppp, space, sizes, 0, posm

yspace = (posm(3) - posm(1))/2.0

posl = posm
posl(2) = posl(2)-posl(0)+space/2 + 0.075
posl(0) = space/2 + 0.075

prx = posl(2) - posl(0)
pmx = (posl(2)+posl(0))/2.0
pry = posl(3) - posl(1)
pmy = (posl(3)+posl(1))/2.0

if (xrange(1)-xrange(0) gt yrange(1)-yrange(0)) then begin
  pry = pry * (yrange(1)-yrange(0))/(xrange(1)-xrange(0))
  posl(3) = pmy + pry/2.0
  posl(1) = posl(3) - pry
endif else begin
  prx = prx * (xrange(1)-xrange(0))/(yrange(1)-yrange(0))
  posl(2) = pmx + prx/2.0
  posl(0) = posl(2) - prx
endelse

plot, xrange, yrange, /nodata, pos = posl, /noerase, $
  xtitle = 'X-GSM (Re)', ytitle = 'Y-GSM (Re)', $
  xstyle =1, ystyle = 1, yrange = -yrange, xrange = [xrange(1), xrange(0)], $
  xtickname = ['',' ','',' ','',' ','',' ','',' ']

oplot, x, y, thick = thick
oplot, [x(1*npts/4)], [y(1*npts/4)], psym = 4
oplot, [x(2*npts/4)], [y(2*npts/4)], psym = 2
oplot, [x(3*npts/4)], [y(3*npts/4)], psym = 5
oplot, [-1000.0,1000.0], [0.0,0.0], linestyle = 1
oplot, [0.0,0.0], [-1000.0,1000.0], linestyle = 1
theta = findgen(19)/18.0*2.0*!pi
oplot, cos(theta), sin(theta)

posr = posm
posr(0) = 1.0 - space/2.0 - (posr(2)-posr(0))
posr(2) = 1.0 - space/2

prx = posr(2) - posr(0)
pmx = (posr(2)+posr(0))/2.0
pry = posr(3) - posr(1)
pmy = (posr(3)+posr(1))/2.0

if (xrange(1)-xrange(0) gt yrange(1)-yrange(0)) then begin
  pry = pry * (yrange(1)-yrange(0))/(xrange(1)-xrange(0))
  posr(3) = pmy + pry/2.0
  posr(1) = posr(3) - pry
endif else begin
  prx = prx * (xrange(1)-xrange(0))/(yrange(1)-yrange(0))
  posr(2) = pmx + prx/2.0
  posr(0) = posr(2) - prx
endelse

plot, xrange, zrange, /nodata, pos = posr, /noerase, $
  xtitle = 'X-GSM (Re)', ytitle = 'Z-GSM (Re)', $
  xstyle =1, ystyle = 1,  xrange = [xrange(1), xrange(0)], $
  xtickname = ['',' ','',' ','',' ','',' ','',' ']

oplot, x, z, thick = thick
oplot, [x(1*npts/4)], [z(1*npts/4)], psym = 4
oplot, [x(2*npts/4)], [z(2*npts/4)], psym = 2
oplot, [x(3*npts/4)], [z(3*npts/4)], psym = 5
oplot, [-1000.0,1000.0], [0.0,0.0], linestyle = 1
oplot, [0.0,0.0], [-1000.0,1000.0], linestyle = 1
oplot, cos(theta), sin(theta)

makect, 'mid'

for nv = 0,nvars-1 do begin

    if nv eq 0 then title = logfilename else title = ""

    ivar = vars(nv)

    get_position, ppp, space, sizes, nv+1, pos, /rect

    pos(0) = pos(0) + 0.075
    pos(1) = pos(1) - yspace
    pos(3) = pos(3) - yspace
    pos(2) = pos(2) - 0.05

    if nv lt nvars-1 then begin
      xtn = strarr(60)+' '
      xt  = ' '
    endif else begin
      xtn = xtickname
      xt  = xtitle
    endelse

    yrange = [min(wlog(*,ivar)),max(wlog(*,ivar))]

;    plot, tmptime, wlog(*,ivar), xstyle=1,		$
;      ytitle = wlognames(ivar),		$
;      xtickname = xtn,			$
;      xtitle = xt,			$
;      xtickv = xtickv,			$
;      xminor = xminor,			$
;      xticks = xtickn,   $
;      xrange = [btr, etr], $
;      pos = pos, /noerase, $
;      thick = thick, yrange = yrange
;
;   oplot, [tmptime(1*npts/4)], [wlog(1*npts/4,ivar)], psym = 4
;   oplot, [tmptime(2*npts/4)], [wlog(2*npts/4,ivar)], psym = 2
;   oplot, [tmptime(3*npts/4)], [wlog(3*npts/4,ivar)], psym = 5

    locnoz = where(power(nv,0:nfft-1,1) gt 0, nnoz)

;stop

    spec = fltarr(nfft,fcount)
    spec(0:nnoz-1,0:fcount-1) = 			$
      (alog10(power(nv,0:nnoz-1,floc))-minspec)*	$
      255.0/(maxspec-minspec)
    loc = where(spec lt 1.0, count)
    if count gt 0 then spec(loc) = 1.0
    loc = where(spec gt 255.0, count)
    if count gt 0 then spec(loc) = 254.0
    xsize = pos(2)-pos(0)-0.001
    ysize = pos(3)-pos(1)-0.001
    tv,spec, pos(0)+0.0005,pos(1)+0.0005,	$
      xsize=xsize,ysize=ysize,/norm

    plot, [btr,etr], mm(freq(floc))*1000.0, pos = pos, /noerase, $
      /nodata, ytitle = 'Frequency (mhz)', $
      xtickname = xtn,			$
      xtitle = xt,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn,   $
      xrange = [btr, etr], xstyle = 1, title = title

    posct = [pos(2)+0.01, pos(1), pos(2)+0.05, pos(3)]
    plotct, 256, posct, [minspec, maxspec], $
      "log(Power("+wlognames(ivar)+"))", /right


endfor

closedevice


end

