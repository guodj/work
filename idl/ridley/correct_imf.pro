
common ffinfo, header

dt = 50.0*60.0
min_slope  = 1.0
min_slopez = 1.0
missing = -1.0e32
re = 6371.2

infile = ''
read, 'Enter input file name : ',infile

outfile = ''
read, 'Enter output file name : ',outfile

openr,2, infile
openw,1, outfile

missing = -1.0e32

ppp = 3
space = 0.03

window,0,xsize=1000,ysize=800

pos_space,ppp,space,sizes, nx = 1

pos = fltarr(ppp+3,4)

for i=0,ppp-1 do begin
  get_position, ppp, space, sizes, i, tpos, 			$
                 xmargin = 0.05, ymargin=0.09, /rect
  pos(i,*) = tpos(*)
endfor

pos(ppp,*)   = [0.00,0.0,0.30,0.05]
pos(ppp+1,*) = [0.33,0.0,0.63,0.05]
pos(ppp+2,*) = [0.70,0.0,1.00,0.05]

col_scal = [0,1,2]
col_imp  = [13,14,15]

itime = intarr(6)

done_with_large_loop = 0

while not done_with_large_loop do begin

  readf,2,itime

  marked = itime(5)
  itime(5) = 0

  c_a_to_r, itime, midtime

  year = itime(0)
  if year gt 1900 then year = year - 1900
  if year gt 100 then year = year - 100

  done = 0

  while (not done) do begin

    c_r_to_a, itime, midtime
    stime = midtime - dt
    etime = midtime + dt

    filename = 'imf_'+chopr('0'+tostr(year),2)

    read_flat_scalor, stime, etime, col_scal, time, data_scal, nrows,	$
	filename = filename

    nrows = nrows(0)

    data = fltarr(nrows)
    time = time - stime

    time_axis, stime, etime, srtime, ertime,        			$
        xtickname1, xtitle1, xtickvalue, xminor, xtickn
    xtickname2 = strarr(10)+' '
    xtitle2 = ' '

    mtime = midtime-stime
    mtime_save = mtime

    plotdumb

    for col = 0, 2 do begin

      if col lt 2 then begin
        xtitle = xtitle2
        xtickname = xtickname2
      endif else begin
        xtitle = xtitle1
        xtickname = xtickname1
      endelse

      data(*) = data_scal(col,*)
      loc = where(data lt -20.0 or data gt 20.0,count)
      if count gt 0 then data(loc) = 100.0

      plot, time(col,*), data,						$
	ystyle = 1, xstyle = 1,						$
	/noerase, pos = pos(col,*), 					$
        yrange = [-20,20], max_value = 20.0,				$
        xtickname = xtickname,						$
        xtickv = xtickvalue, xticks = xtickn,				$
        xminor = xminor, xtitle=xtitle,					$
	xrange = [srtime, ertime]

      oplot, [srtime,ertime], [0.0,0.0], linestyle = 2
      oplot, [mtime,mtime], [-20.0,20.0]

    endfor

    for i=ppp,ppp+2 do begin

      x = [0,1]
      y = [0,1]
      plot, x, y, /nodata, /noerase, xstyle = 5, ystyle = 5, pos = [0,0,1,1]

      x = [pos(i,0), pos(i,0), pos(i,2), pos(i,2), pos(i,0)]
      y = [pos(i,1), pos(i,3), pos(i,3), pos(i,1), pos(i,1)]
      oplot, x, y

      if i eq ppp then title = 'Reduce 2 min'
      if i eq ppp+1 then title = 'Add 2 min'
      if i eq ppp+2 then title = 'Done'

      xyouts, mean(x(0:3)), mean(y(0:3)), title, alignment = 0.5

    endfor

    command = 0

    basetime = min(time)
    picked   = 0
    utsec    = fltarr(100)

    get_mouse_command, pos, command, ret_pos

    if command eq ppp then midtime = midtime - 2.0*60.0
    if command eq ppp+1 then midtime = midtime + 2.0*60.0
    if command eq ppp+2 then done = 1

  endwhile

  filename_2 = 'imf_'+chopr('0'+tostr(year),2)+'_cross'

  col_vect_2 = [46]+marked

  mtime = midtime

  done = 0

  direction = -1.0
  ntimes = 0

  while not done do begin

    stime = mtime - 60.0
    etime = mtime + 60.0

    read_flat_vector, stime, etime, col_vect_2, time_v, data_vect, 	$
	nrows, filename = filename_2

    loc = where(data_vect ne missing, count)

    if count gt 0 then done = 1						$
    else begin
      if ntimes lt 30 then begin
        print, 'Could not find time, stepping Backwards.'
      endif else begin
        if ntimes eq 30 then mtime = midtime
        print, 'Could not find time, stepping Forward.'
        direction = 1.0  
      endelse
      ntimes = ntimes + 1
      mtime = mtime + direction*120.0
    endelse

  endwhile

  col_scal_2 = [45,4,5,6,17,18,19,7]

  read_flat_scalor, stime, etime, col_scal_2, time_2, data_scal_2,	$
	nrows_2, filename = filename_2

  n_e = n_elements(data_vect(0,*,0))
  maxc = max(data_vect(0,n_e-1,*))
  loc = where(data_vect(0,n_e-1,*) eq maxc) 
  t = 2.0*(loc(0)-30.0)*60.0 + data_scal_2(0,n_e-1)*60.0

  mx_wind = data_scal_2(1,n_e-1)
  my_wind = data_scal_2(2,n_e-1)
  mz_wind = data_scal_2(3,n_e-1)

  mx_imp8 = data_scal_2(4,n_e-1)
  my_imp8 = data_scal_2(5,n_e-1)
  mz_imp8 = data_scal_2(6,n_e-1)

  dx = (mx_wind - mx_imp8)*re
  dy = (my_wind - my_imp8)*re
  dz = (mz_wind - mz_imp8)*re

  vx = data_scal_2(7,n_e-1)

  loc = where(time(0,*) eq mtime_save,count)

  if count gt 0 then begin

    bx = mean(data_scal(0,loc(0)-5:loc(0)),missing)
    by = mean(data_scal(1,loc(0)-5:loc(0)),missing)
    bz = mean(data_scal(2,loc(0)-5:loc(0)),missing)

    bx_after = mean(data_scal(0,loc(0)+3:loc(0)+8),missing)
    by_after = mean(data_scal(1,loc(0)+3:loc(0)+8),missing)
    bz_after = mean(data_scal(2,loc(0)+3:loc(0)+8),missing)

  endif else begin

    bx = 0.0
    by = 0.0
    bz = 0.0

    bx_after = 0.0
    by_after = 0.0
    bz_after = 0.0

  endelse

  b_before = [bx,by,bz]
  b_after  = [bx_after,by_after,bz_after]
  b_cross  = crossp(b_before, b_after)

; compute propagation times:

; first x distance

  dt_x = (dx/(-1.0*vx))/60.0

; then parker-spiral

  slope = -1.0

  magx = (-1.0*dy + slope*mx_wind*re)/slope
  dt_p = ((magx-mx_imp8*re)/(-1.0*vx))/60.0

; then magnetic field

  dt_b = fltarr(6)

  for i=0,n_elements(dt_b)-1 do begin

; Use the xy magnetic field (10 minutes before the event):

    if i eq 0 then begin
      bx = b_before(0)
      by = b_before(1)
      bz = 100000.0
    endif

; Use the total magnetic field (10 minutes before the event):

    if i eq 1 then begin
      bx = b_before(0)
      by = b_before(1)
      bz = b_before(2)
    endif

; Use the xy magnetic field (after the event):

    if i eq 2 then begin
      bx = b_after(0)
      by = b_after(1)
      bz = 100000.0
    endif

; Use the total magnetic field (after the event):

    if i eq 3 then begin
      bx = b_after(0)
      by = b_after(1)
      bz = b_after(2)
    endif

; Use the xy crossed magnetic field :

    if i eq 4 then begin
      bx = b_cross(0)
      by = b_cross(1)
      bz = 100000.0
    endif

; Use the total crossed magnetic field :

    if i eq 5 then begin
      bx = b_cross(0)
      by = b_cross(1)
      bz = b_cross(2)
    endif

    if (bx ne 0.0) then begin
      yslope = by/bx
      zslope = bz/bx
    endif else begin
      yslope = by/0.01
      zslope = bz/0.01
    endelse

    if (abs(yslope) gt min_slope) and (bx ne missing) then		$
      ymagx = (-1.0*dy + yslope*mx_wind*re)/yslope			$
    else begin
      if yslope eq 0.0 then yslope = 10000.0				$
      else yslope = min_slope*yslope/abs(yslope)
      ymagx = (-1.0*dy + yslope*mx_wind*re)/yslope
    endelse

    if (abs(zslope) gt min_slope) and (bx ne missing) then		$
      magx = (-1.0*dz + zslope*ymagx)/zslope				$
    else begin
      if zslope eq 0.0 then zslope = 10000.0				$
      else zslope = min_slopez*zslope/abs(zslope)
      magx = (-1.0*dz + zslope*ymagx)/zslope
    endelse

    if (bx ne missing) then begin
      dt_b(i) = ((magx-mx_imp8*re)/(-1.0*vx))/60.0
    endif

  endfor

  done = 0

  while (not done) do begin

    stime = midtime - dt + t
    etime = midtime + dt + t

    filename = 'imf_'+chopr('0'+tostr(year),2)

    read_flat_scalor, stime, etime, col_imp, time_imp, data_imp, nrows,	$
	filename = filename

    nrows = nrows(0)

    imp = fltarr(nrows)
    time_imp = time_imp - stime

    mtime = midtime - (stime - t)

    plotdumb

    for col = 0, 2 do begin

      if col lt 2 then begin
        xtitle = xtitle2
        xtickname = xtickname2
      endif else begin
        xtitle = xtitle1
        xtickname = xtickname1
      endelse

      data(*) = data_scal(col,*)
      imp(*) = data_imp(col,*)

      loc = where(data lt -20.0 or data gt 20.0,count)
      if count gt 0 then data(loc) = 100.0

      loc = where(imp lt -20.0 or imp gt 20.0,count)
      if count gt 0 then imp(loc) = 100.0

      plot, time(col,*), data,						$
	ystyle = 1, xstyle = 1,						$
	/noerase, pos = pos(col,*), 					$
        yrange = [-20,20], max_value = 20.0,				$
        xtickname = xtickname,						$
        xtickv = xtickvalue, xticks = xtickn,				$
        xminor = xminor, xtitle=xtitle,					$
	xrange = [srtime, ertime]

      oplot, time_imp(col,*), imp, linestyle = 1, max_value = 20.0

      oplot, [srtime,ertime], [0.0,0.0], linestyle = 2
      oplot, [mtime,mtime], [-20.0,20.0]

    endfor

    for i=ppp,ppp+2 do begin

      x = [0,1]
      y = [0,1]
      plot, x, y, /nodata, /noerase, xstyle = 5, ystyle = 5, pos = [0,0,1,1]

      x = [pos(i,0), pos(i,0), pos(i,2), pos(i,2), pos(i,0)]
      y = [pos(i,1), pos(i,3), pos(i,3), pos(i,1), pos(i,1)]
      oplot, x, y

      if i eq ppp then title = 'Reduce 2 min'
      if i eq ppp+1 then title = 'Add 2 min'
      if i eq ppp+2 then title = 'Done'

      xyouts, mean(x(0:3)), mean(y(0:3)), title, alignment = 0.5

    endfor

    command = 0

    get_mouse_command, pos, command, ret_pos

    if command eq ppp then t = t - 2.0*60.0
    if command eq ppp+1 then t = t + 2.0*60.0
    if command eq ppp+2 then done = 1

    if command ge ppp+3 then begin
      done = 1
      done_with_large_loop = 1
    endif

  endwhile

  t = t/60.0

  c_r_to_a, itime, midtime
  itime(5) = marked
  printf,1, itime
  printf,1, dx/re, dy/re, dz/re, b_before(0), b_before(1), b_before(2)
  printf,1, t, dt_x, dt_p
  printf,1, dt_b

  print, t, dt_x, dt_p
  print, dt_b

  if eof(2) then done_with_large_loop = 1

endwhile

close,1,2
wdelete

end
