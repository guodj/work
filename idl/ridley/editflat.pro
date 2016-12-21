pro get_command, pymax, xmid, pos, npts, nsubpts, command, range

  out = 0

  while not out do begin

    cursor,x,y,2, /norm
    mouse = !err

    if mouse eq 1 then begin

      if y gt pymax and x gt xmid then begin
	command = 0
	out = 1
      endif

      if y gt pymax and x lt xmid then begin
	command = 1
	out = 1
      endif

      if y lt pymax then begin

        l = where( (x gt pos(*,0)) and (x lt pos(*,2)) and 		$
		   (y gt pos(*,1)) and (y lt pos(*,3)), count)

	if count gt 0 then begin

	  command = l(0)+2
	  midr = (x-pos(0,0))/(pos(0,2)-pos(0,0))*float(npts)
	  if midr lt nsubpts then midr = nsubpts
	  range = intarr(2)
	  range(0) = midr-nsubpts
	  range(1) = midr+nsubpts
	  if range(1) ge npts then range(1) = npts-1

	  out = 1

	endif
  
      endif

      clear_mouse

    endif

  endwhile

  return

end

;  This procedure is designed to take output files from imf_dbdt.for and
;  find the times in a magnetometer file, then output 1 hour before the
;  events and 2 hours after the events.

;pro edit_flat

;  magfile = ''
;  print,'Enter File to edit : '
;  read,magfile

  get_files, datafile,n,1,'Enter flatfile to edit',ip,'HED'
  magfile = datafile(0)

;  fileout = ''
;  print,'Enter file to output : '
;  read,fileout

  fileout = getstr('Output Filename')

  setup_var_menu, 1, magfile, timeint, name, unit, nvar, nrows
  magfile = pad_string(magfile, 100)
  fileout = pad_string(fileout, 100)

  source = strarr(nvar(0)+1)+'edited data'

  fdnfile = 0.0
  fdids = fltarr(6,6)
  data = fltarr(100)
  inunit = 0.0

; open all of the data files :
;  first one, is the new one, which we have to create, brand new:

  outunit = ffcreate(fileout,name,unit,source,fdnfile,fdids)

;  The next one is the old magnetometer data file:

  z = call_external('fopen_exe','fopen',		$
	magfile, fdnfile, fdids, inunit)

;  the final one is the outyymm.dat file, which contains the times of the
;    events which we want to take out of the magnetometer file :

  nvart = nvar(0)
  setup_file, fbase, 'Choose columns to edit', name, 0
  get_var, fbase, name, 0, nvart
  print, 'Variables selected :'
  nsel = 0
  vsel = intarr(nvar(0)+1)
  for i=0,nvar(0) do if strpos(name(i),'(selected)') gt 0 then begin
    vsel(nsel) = i
    nsel = nsel + 1
    print, '     ',strcompress(strmid(name(i),0,11))
  endif
  nsel = nsel - 1
  if nsel gt -1 then vsel = vsel(0:nsel)

  space = 0.005
  ybot = 0.05
  ytop = 0.91
  yint = ((ytop-ybot) - space*float(nsel))/float(nsel+1)

  xleft = 0.05
  xright = 0.995

  ut1 = double(0.0)
  ut2 = double(0.0)
  row1 = 1.0
  row2 = 2.0

  pos = fltarr(nsel+1,4)
  for i=0,nsel do					$
    pos(i,*) = [xleft, ytop-(yint+space)*float(i)-yint,	$
		xright,ytop-(yint+space)*float(i)]

  z = call_external('fget_exe','fget',			$
	  fdnfile, fdids, inunit, data, ut1, row1)
  z = call_external('fget_exe','fget',			$
	  fdnfile, fdids, inunit, data, ut2, row2)

  dt = ut2-ut1
  print, 'data is at a resolution of : ',tostr(dt),' seconds
  ppd = 24.0*3600.0/dt
  print, 'npts per day : ',tostr(ppd)

  snptsmax = getstr('Number of minutes at a time to edit')
  nptsmax = float(snptsmax)*60.0/dt

  day = 1
  range = [1,fix(nrows(0)/ppd)]
  if range(0) ne range(1) then 				$
    day = get_integer(range, day,title='Choose start day')

  done = 0

  daydata = fltarr(nvar(0)+1,ppd+1)
  daytime = dblarr(ppd+1)
  missing = -1.0e32

  row = float(day-1)*ppd+1.0
  outrow = 1.0

  window,1,xsize=1000,ysize=800

  while not done do begin

    rowe = row + ppd
    if rowe gt nrows(0) then rowe = float(nrows(0))

    for crow = row, rowe do begin

      z = call_external('fget_exe','fget',			$
	fdnfile, fdids, inunit, data, ut1, crow)
      daydata(0:nvar(0),crow-row) = data(2:nvar(0)+2)
      daytime(crow-row) = ut1

    endfor

    next = 0

    while not next do begin

     plotdumb

      pymax = (max(pos(*,3))+1.0)/2.0
    
      xmid = (xleft+xright)/2.0

      x = [xleft,xleft,xmid-0.002,xmid-0.002,xleft]
      y = [pymax,0.999,0.999,pymax,pymax]
      oplot, x, y
      xyouts, mean([xleft,xmid-0.002]),				$
	    pymax+0.004,   					$
	    'Press Here to Go to NEXT time period',		$
	    alignment = 0.5, charsize = 1.0

      x = [xright,xright,xmid+0.002,xmid+0.002,xright]
      y = [pymax,0.999,0.999,pymax,pymax]
      oplot, x, y
      xyouts, mean([xright,xmid+0.002]),				$
	    pymax+0.004, 					$
	    'Press Here to END',				$
	    alignment = 0.5, charsize = 1.0

      xyouts, xmid, max(pos(*,3))+0.004,	 			$
	'Select time period to edit by pressing mouse near center of time', $
	    alignment = 0.5, charsize = 1.5

      loc = where(daytime gt 0, count)
      if count gt 0 then begin
        c_r_to_a, itime, daytime(0)
        c_a_to_s, itime, sdate
      endif

      for i=0,nsel do begin

        dumbdata = daydata(vsel(i),0:ppd-1)
        loc = where(dumbdata eq missing,count)
        if count gt 0 then dumbdata(loc) = -1.0*missing
        xlab = strarr(20)
        if i ne nsel then begin
	  xlab = xlab + ' '
	  xt = ' '
        endif else xt = 'Minutes after '+sdate+' UT'
        loc = where(daytime gt 0, count)

        if count gt 0 then begin

	  loc2 = where(dumbdata(loc) lt -0.2*missing,count)
	  if count eq 0 then loc2 = [0]
	  mmdata = mm(dumbdata(loc(loc2)))
	  drange = mmdata(1) - mmdata(0)
	  mmdata(0) = mmdata(0) - 0.1 * drange
	  mmdata(1) = mmdata(1) + 0.1 * drange
          plot, (daytime(loc)-daytime(loc(0)))/60.0, dumbdata(loc), 	$
	    /noerase, pos = pos(i,*), xstyle = 1, 			$
	    max_value = -0.2*missing,	$
	    xtickn = xlab, xtitle = xt, yrange = mmdata, ystyle=1,	$
	    ytitle = strcompress(strmid(name(vsel(i)),0,11))

        endif

      endfor

      get_command, pymax, xmid, pos, count, nptsmax, command, range
      if command le 1 then next = 1
      if command eq 0 then done = 1

      if command ge 2 then begin

	ddata = daydata(vsel(command-2), range(0):range(1))
	dtime = daytime(range(0):range(1))
	edit_data, ddata, dtime-dtime(0)
	daydata(vsel(command-2), range(0):range(1)) = ddata

      endif

    endwhile

    for crow = row, rowe do begin

      data(2:nvar(0)+2) = daydata(0:nvar(0),crow-row)
      data(0:1) = 0.0
      loc = where(abs(data) gt 100000.0, count)
      if count gt 0 then data(loc) = missing

      ut1 = daytime(crow-row)

      z = call_external('fput_exe','fput',			$
	  fdnfile,fdids,outunit,data,ut1,outrow)

      outrow=outrow+1.0

    endfor

    if not done then row = row + ppd

  endwhile

  z = call_external('fclose_exe','fclose',fdnfile,fdids,inunit)
  z = call_external('fclose_exe','fclose',fdnfile,fdids,outunit)

;  return

end

