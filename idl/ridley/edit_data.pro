pro edit_data, data, time

  window, 1, xsize=1000, ysize = 800

  missing = -1.0e32
  mv = -0.2*missing

  loc = where(data lt mv,count)
  if count eq 0 then loc = [0]
  mmdata = mm(data(loc))
  range = mmdata(1) - mmdata(0)
  mmdata(0) = mmdata(0) - 0.1*range
  mmdata(1) = mmdata(1) + 0.1*range

  plot, time, data, psym = -4, yrange = mmdata, ystyle=1, xstyle=1,	$
	pos = [0,0,1,1], max_value=mv
  x = [0,-1,0,1,0]
  y = [-1,0,1,0,-1]
  usersym,x,y,/fill

  out = 0
  first = 0

  pts = intarr(2)

  xyouts, 0.05, 0.95, 'Right mouse to Select first point', /norm
  xyouts, 0.5, 0.95, 'Middle mouse to skip to delete points', 		$
	alignment = 0.5, /norm
  xyouts, 0.95, 0.95, 'Left mouse to end', alignment = 1.0, /norm

  while out eq 0 do begin

    cursor,x,y,2
    mouse = !err

    if mouse eq 1 then begin

      if first gt 0 then out = 1
      xd = time - x
      yd = data - y
      dis = sqrt(xd^2.0 + yd^2.0)
      loc = where(dis eq min(dis), count)
      if count ne 0 then begin
	oplot, [time(loc(0))],[data(loc(0))], psym = 8
	pts(first) = loc(0)
	first = first + 1
        xyouts, 0.05, 0.95, 'Right mouse to Select first point', /norm,color=0
        xyouts, 0.05, 0.95, 'Right mouse to Select second point', /norm
	while mouse eq 1 do begin
	  cursor,x,y,0, /device
	  mouse = !err
	endwhile

      endif else print, 'Error in routine...'

    endif

    if mouse eq 2 then begin
      pts(0) = n_elements(data)/3
      pts(1) = 2*n_elements(data)/3
      out = 1
      xyouts, 0.05, 0.95, 'Right mouse to Select first point', /norm, color=0
    endif

    if mouse gt 2 then out = 2

  endwhile

  out = out-1

  if out eq 0 then begin

    nt = n_elements(time)-1
    mn1 = mean(data(0:pts(0)))
    mn = mean(data(pts(0):pts(1)))
    mn_old = mn
    mn2 = mean(data(pts(1):nt))

    xs = float(!d.x_size)
    x1 = xs*(time(pts(0)) - time(0))/(time(nt)-time(0)) - 2
    x2 = xs*(time(pts(1)) - time(0))/(time(nt)-time(0)) - x1 + 4

    if x2 gt xs then x2 = xs

    ys = float(!d.y_size)
    y1 = ys*(mn-mmdata(0))/(mmdata(1)-mmdata(0)) - 2
    y2 = 4

    pict = tvrd(x1,y1,x2,y2)

    oplot, [time(0),time(pts(0))],[mn1,mn1]
    oplot, [time(pts(0)),time(pts(1))],[mn,mn]
    oplot, [time(pts(1)),time(nt)],[mn2,mn2]

    xyouts, 0.05, 0.95, 'Right mouse to Select second point', /norm, color = 0
    xyouts, 0.95, 0.95, 'Left mouse to end', alignment = 1.0, /norm, color = 0
    xyouts, 0.5, 0.95, 'Middle mouse to skip to delete points', 	$
	alignment = 0.5, /norm, color = 0
    xyouts, 0.05, 0.95, 'Left mouse to edit mean', /norm
    xyouts, 0.5, 0.95, 'Second mouse to delete data point', 		$
	alignment = 0.5, /norm
    xyouts, 0.95, 0.95, 'Right mouse to move to next time', 		$
	alignment = 1.0, /norm

    y1_old = y1

  endif

  while not out do begin

    cursor,x,y,2, /device
    mouse = !err

    if mouse eq 1 then begin

      while mouse eq 1 do begin

	if y-2 ne y1 then begin
	  tv, pict, x1, y1 
	  y1 = y - 2
	  pict = tvrd(x1,y1,x2,y2)
	  mn = (y1 + 2)*(mmdata(1)-mmdata(0))/ys + mmdata(0)
	  oplot, [time(pts(0)),time(pts(1))],[mn,mn]
	endif
	cursor,x,y,0, /device
	mouse = !err

      endwhile

      data(pts(0):pts(1)) = data(pts(0):pts(1)) - mn_old + mn
      mn_old = mn
      plot, time, data, psym = -4, yrange = mmdata, ystyle=1, xstyle=1,	$
	pos = [0,0,1,1], max_value=mv
      xyouts, 0.05, 0.95, 'Left mouse to edit mean', /norm
      xyouts, 0.5, 0.95, 'Second mouse to delete data point', 		$
	alignment = 0.5, /norm
      xyouts, 0.95, 0.95, 'Right mouse to move to next time', 		$
	alignment = 1.0, /norm

      pict = tvrd(x1,y1,x2,y2)

      oplot, [time(0),time(pts(0))],[mn1,mn1]
      oplot, [time(pts(0)),time(pts(1))],[mn,mn]
      oplot, [time(pts(1)),time(nt)],[mn2,mn2]

    endif

    if mouse eq 2 then begin

      cursor,x,y,0

      if first gt 0 then out = 1
      xd = time - x
      yd = data - y
      dis = sqrt(xd^2.0 + yd^2.0)
      loc = where(dis eq min(dis), count)
      if count ne 0 then data(loc(0)) = -1.0*missing

      plot, time, data, psym = -4, yrange = mmdata, ystyle=1, xstyle=1,	$
	pos = [0,0,1,1], max_value=mv
      xyouts, 0.05, 0.95, 'Left mouse to edit mean', /norm
      xyouts, 0.5, 0.95, 'Second mouse to delete data point', 		$
	alignment = 0.5, /norm
      xyouts, 0.95, 0.95, 'Right mouse to move to next time', 		$
	alignment = 1.0, /norm

      pict = tvrd(x1,y1,x2,y2)

      oplot, [time(0),time(pts(0))],[mn1,mn1]
      oplot, [time(pts(0)),time(pts(1))],[mn,mn]
      oplot, [time(pts(1)),time(nt)],[mn2,mn2]

      clear_mouse

      out = 0

    endif

    if mouse eq 4 then out = 1

  endwhile

  return

end

