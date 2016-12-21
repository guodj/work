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

  print, ' '
  print, 'Npts : ',npts,'  dt : ',dt

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
    amp     = fltarr(nv,npts)
    power   = fltarr(nv,npts)
    t       = fltarr(npts)
    freq    = fltarr(npts)
    period  = fltarr(npts)
    df      = 1.0/(dt*float(npts))

    for i=srow, erow do begin

      ii = i - srow
      row = float(i)

      z = call_external('fget_exe','fget',fdnfile,fdids,		$
		inunit,data,utdata,row)
      t(ii) = utdata - suttime
      for j=0,nv-1 do							$
	signal(j,ii) = data(col(j))

      if ii le npts/2 then freq(ii) = float(ii)*df 			$
      else freq(ii) = float(npts/2 - ii)*df
    
    endfor

    period(1:npts/2) = 1.0/freq(1:npts/2)

    fileout = ''
    print, 'Enter output flat-file name :'
    read, fileout

    tname = strarr(nv+2)
    tname(0) = 'Frequency'
    tname(1) = 'Period'
    for i=0,nv-1 do tname(i+2) = vars(col(i)-2)

    tunit = strarr(nv+2)
    tunit(0) = 'Hz'
    tunit(1) = 'sec'
    for i=0,nv-1 do tunit(i+2) = units(col(i)-2)+'^2/Hz'

    source = strarr(nv+2)+'fftspec.pro'

    outunit = ffcreate(fileout,tname,tunit,source,fdnfile,fdids)

    M = n_elements(signal(0,*))-1
    c_window, M, wind,3

    for i=0,nv-1 do begin

      dumsig = signal(i,*)

      dumamp   = fft(dumsig*wind,-1)
      dumpower = (2.0*abs(dumamp))^2.0/df
      amp(i,*)   = 2.0*abs(dumamp(*))
      power(i,*) = dumpower(*)

    endfor

    for i=1,npts/2 do begin

      data(2) = freq(i)
      data(3) = period(i)

      for j=0,nv-1 do data(j+4) = power(j,i)

      row = float(i)
      utdata =  suttime + double(dt*2.0*(i-1))

      z = call_external('fput_exe','fput',				$
	fdnfile,fdids,outunit,data,utdata,row)

    endfor

    z = call_external('fclose_exe','fclose',fdnfile,fdids,outunit)

  endif

  z = call_external('fclose_exe','fclose',fdnfile,fdids,inunit)

  !p.multi=0

end    
