pro remove_mean, data, mean

  npts = n_elements(data)

  mean = 0.0

  for i=0,npts-1 do mean = mean + data(i)

  if npts gt 0 then begin
    mean = mean/npts
    data = data - mean
  endif

  return

end

pro p_color_t, min, max

  xp1 = 0.97
  xp2 = 1.0
  yp1 = 0.5
  yp2 = 0.9

  pos = [xp1,yp1,xp2,yp2]

  minsp = 10.0^min
  maxsp = 10.0^max

  plot_io, [0,1], [minsp,maxsp],	$
	xstyle=5,			$
	ystyle=1,			$
	pos = pos,			$
	charsize = 0.75,		$
	/nodata,			$
	/noerase

  oplot, [0,1], [minsp,minsp]
  oplot, [0,1], [maxsp,maxsp]

  xp1 = 0.9705
  xp2 = 0.9995
  yp1 = 0.5005
  yp2 = 0.8995

  pos = [xp1,yp1,xp2,yp2]

  plot, [0,1], [0,201],			$
	xstyle=5,			$
	ystyle=5,			$
	pos = pos,			$
	/nodata,			$
	/noerase

   x = [0,0,1,1,0]
   y = [0,1,1,0,0]

   for i=0,200 do polyfill, x, y+i, color = i

   xyouts, 1.5, 100.0, 'Power (nT!E2!N/hz)', 	$
	alignment = 0.5, orientation = 270.0, 	$
	charsize = 0.75

  return

end

;pro ffttest

  data = fltarr(100)
  utdata = double(0.0)

  infile = ''
  print, 'Enter file to input (return for widget findfile) : '
  read, infile
  if strlen(infile) eq 0 then begin
    df = strarr(1)
    nf = 1
    mx = 1
    ft = 'Choose file to FFT'
    ip = ''
    ex = 'HED'
    get_files, df,nf,mx,ft,ip,ex
    infile = df(0)
    print, infile
  endif
  setup_var_menu, 1, infile, timeint, vars, units, nvars, nrows
  infile = pad(infile,100)
  fdnfile = 0.0
  fdids = fltarr(6,6)
  inunit = 0.0
  z = call_external('fopen_exe','fopen',infile,fdnfile,fdids,inunit)

  print, 'Enter ps file name :'
  psfile = ''
  read,psfile

  print, ' '
  print, 'Start time : ', timeint(0,0)
  stime = ''
  read, 'Enter starting time (return for same) : ',stime
  if strlen(stime) eq 0 then stime = timeint(0,0)
  c_s_to_a, sitime, stime
  c_a_to_r, sitime, suttime
  srow = 0.0
  z = call_external('fsearch_exe','fsearch',fdnfile,fdids,		$
	inunit, suttime, srow)
  z = call_external('fget_exe','fget',fdnfile,fdids,			$
	inunit,data,utdata,srow)
  suttime = utdata

  print, ' '
  print, 'Ending time : ', timeint(0,1)
  etime = ''
  read, 'Enter ending time (return for same) : ',etime
  if strlen(etime) eq 0 then etime = timeint(0,1)
  c_s_to_a, eitime, etime
  c_a_to_r, eitime, euttime
  erow = 0.0
  z = call_external('fsearch_exe','fsearch',fdnfile,fdids,		$
	inunit, euttime, erow)
  z = call_external('fget_exe','fget',fdnfile,fdids,			$
	inunit,data,utdata,erow)
  euttime = utdata
  npts = fix(erow - srow) + 1
  dt = (euttime-suttime)/float(npts-1)

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

  for i=0,nvars(0) do print, i+2,' - ',vars(i)

  nv = 0
  done = 0

  while not done do begin

    que = ''
    print, 'Enter variable ',tostr(nv+1), ' (return to end) : '
    read, que
    if strlen(que) gt 0 then begin
      if nv eq 0 then col = fix(que)					$
      else col = [col, fix(que)]
      nv = nv + 1
    endif else done = 1

  endwhile

  if nv gt 0 then begin

    signal  = fltarr(nv,npts)
    t       = fltarr(npts)
    power   = fltarr(nv,nfft,ntime)
    sigmean = fltarr(nv,nfft)
    freq    = fltarr(ntime)
    period  = fltarr(ntime)
    df      = 1.0/(dt*float(ntime))

    for i=srow, erow do begin

      ii = i - srow
      row = float(i)

      z = call_external('fget_exe','fget',fdnfile,fdids,		$
		inunit,data,utdata,row)
      t(ii) = utdata - suttime
      for j=0,nv-1 do							$
	signal(j,ii) = data(col(j))

    endfor

    print, 'Done reading data'

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

  endif

  signal = signal(*,ntime/2:npts-1-ntime/2)

  if strlen(psfile) gt 0 then begin
    !p.font = 0
    orient = 'p'
    psfont = 28
    percent = 0.92
    setdevice, psfile, orient, psfont, percent
  endif

  set_up_color

  plotdumb

  times = intarr(1,2,6)
  placement = 0
  curvar = 0
  st = double(t(ntime/2)) + suttime
  c_r_to_a, itime, st
  times(0,0,*) = itime
  basetime = st
  et = double(t(npts-1-ntime/2)) +suttime
  c_r_to_a, itime, et
  times(0,1,*) = itime
  moncheck = 'JanFebMarAprMayJunJulAugSepOctNovDec'
  compute_axis, times, placement, btr, etr, curvar, xlab, xtl,		$
		basetime, moncheck, nl, ticktime, nminor

  pos = fltarr(nv,4)
  space = 0.0075
  xdis = (1.0 - space)/float(nv)

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
   print, tostr(100*ftempcount/fcount),' % between 10e',tostr(i),' and 10e', $
	tostr(i+1)
  endfor
  print, 'Enter min,max values for color table (ex 1.0,1.0e8) :'
  read, minspec, maxspec
  minspec = log(minspec)
  maxspec = log(maxspec)

  print, ' '
  print, 'Npts : ',npts,'  dt : ',dt

  title = ''
  print, 'Enter title (return for no title) :'
  read, title

  for i=0,nv-1 do begin

    if i ne 0 then begin
      relxlab = strarr(n_elements(xlab))+' '
      relxtl = ' '
      relylab = strarr(10)
      relylab(0) = ' '
    endif else begin
      relxlab = xlab
      relxtl = xtl
      relylab = strarr(10)
    endelse

    if i eq nv-1 then reltl = title else reltl = ' '

    xp1 = 0.05
    xp2 = 0.9
    yp1 = float(i)*xdis+space
    yp2 = float(i+1)*xdis

    tl = strcompress(vars(col(nv-1-i)-2),/remove_all)
    ytl = 'Frequency (mhz)'

    pos(i,*) = [xp1,yp1,xp2,yp2]

    if fcount gt 0 then begin

      plot, [btr,etr], mm(freq(floc))*1000.0,	$
	xticks = nl,				$
	xtickname = relxlab,			$
	xtitle = relxtl,			$
	xstyle = 1,				$
	xminor = nminor,			$
	xticklen = -0.02,			$
	xtickv = ticktime,			$
	ytickname = relylab,			$
	ytitle = ytl,				$
	ystyle = 1,				$
	yticklen = -0.005,			$
	pos = pos(i,*),				$
	title = reltl,				$
	/noerase,				$
	/nodata

      xp = 0.91
      yp = (yp2+yp1)/2.0
      xyouts, xp, yp, tl, alignment=0.5, orientation=270, /norm

      locnoz = where(power(nv-1-i,0:nfft-1,1) gt 0, nnoz)

      spec = fltarr(nfft,fcount)
      spec(0:nnoz-1,0:fcount-1) = 			$
	(log(power(nv-1-i,0:nnoz-1,floc))-minspec)*	$
	200.0/(maxspec-minspec)
      loc = where(spec lt 1.0, count)
      if count gt 0 then spec(loc) = 1.0
      loc = where(spec gt 200.0, count)
      if count gt 0 then spec(loc) = 200.0
      xsize = pos(i,2)-pos(i,0)-0.001
      ysize = pos(i,3)-pos(i,1)-0.001
      tv,spec, pos(i,0)+0.0005,pos(i,1)+0.0005,	$
	  xsize=xsize,ysize=ysize,/norm

    endif

  endfor

  p_color_t, minspec, maxspec

  z = call_external('fclose_exe','fclose',fdnfile,fdids,inunit)
  if !d.name eq 'PS' then begin

    device, /close
    set_plot, 'X'

    b = byte(psfile)
    loc = where(b eq 46, count)
    if count eq 0 then loc(0) = strlen(psfile)
    psfile = strmid(psfile,0,loc(0))+'_mean.ps'

    !p.font = 0
    orient = 'p'
    psfont = 28
    percent = 0.92
    setdevice, psfile, orient, psfont, percent

  endif

  plotdumb

  x1 = findgen(n_elements(signal(0,*)))*etr/float(n_elements(signal(0,*)))
  x2 = findgen(nfft)*etr/float(nfft)
  
  for i=0,nv-1 do begin

    if i ne 0 then begin
      relxlab = strarr(n_elements(xlab))+' '
      relxtl = ' '
      relylab = strarr(10)
      relylab(0) = ' '
    endif else begin
      relxlab = xlab
      relxtl = xtl
      relylab = strarr(10)
    endelse

    if i eq nv-1 then reltl = title else reltl = ' '

    xp1 = 0.05
    xp2 = 0.9
    yp1 = float(i)*xdis+space
    yp2 = float(i+1)*xdis

    ytl = strcompress(vars(col(nv-1-i)-2),/remove_all)

    pos(i,*) = [xp1,yp1,xp2,yp2]

    y = signal(nv-1-i,*)
    if !d.name eq 'PS' then color = 0 else color = 255

      plot, x1, y,				$
	xticks = nl,				$
	xtickname = relxlab,			$
	xtitle = relxtl,			$
	xstyle = 1,				$
	xminor = nminor,			$
	xticklen = -0.02,			$
	xtickv = ticktime,			$
	ytickname = relylab,			$
	ytitle = ytl,				$
	yticklen = -0.005,			$
	title = reltl,				$
	pos = pos(i,*),				$
	color = color,				$
	/noerase

    y = sigmean(nv-1-i,*)
    oplot, x2,y, color = color, linestyle = 2

      xp = 0.91
      yp = (yp2+yp1)/2.0
      xyouts, xp, yp, 'nT', alignment=0.5, orientation=270, /norm

  endfor

  if !d.name eq 'PS' then begin

    device, /close
    set_plot, 'X'

  endif

end    
