
; range is distance away from the center point (i.e. colatitude)
; theta must be in radians!!

pro plot_vectors_polar, east, north, range, theta, normalization, 	$
	arrow = arrow, no0 = no0, length = length, pos = pos,           $
        title = title, skip = skip

  if n_elements(arrow) gt 0 then arrow = 1 else arrow = 0
  if n_elements(no0) gt 0 then no0 = 1 else no0 = 0
  if n_elements(length) eq 0 then unit_length = 0 else unit_length = length
  if n_elements(title) eq 0 then title = ''
  if n_elements(pos) eq 0 then pos = [-1e32,-1e32]
  if n_elements(skip) eq 0 then skip = 1

  nx = n_elements(range(*,0))
  ny = n_elements(range(0,*))

  for i=0,nx-1, skip do for j=0,ny-1, skip+1 do begin

    e = east(i,j)*normalization
    n = north(i,j)*normalization
    t = theta(i,j)
    r = range(i,j)

    if (r gt 0) or (not no0) then begin

      x0 = r*cos(t)
      y0 = r*sin(t)

      x = -n*cos(t) - e*sin(t)
      y = -n*sin(t) + e*cos(t)

      plots, [x0,x0+x], [y0,y0+y]

      if (arrow) then begin

        length = (e^2+n^2)^0.5
        t2     = asin(y/length)
        if (x lt 0.0) then t2 = !pi - t2

        x1 = 0.5*length*cos(t2 + 15.0*!pi/180.0)
        y1 = 0.5*length*sin(t2 + 15.0*!pi/180.0)

        x2 = 0.5*length*cos(t2 - 15.0*!pi/180.0)
        y2 = 0.5*length*sin(t2 - 15.0*!pi/180.0)

        plots, [x0+x, x0+x-x1], [y0+y,y0+y-y1]
        plots, [x0+x, x0+x-x2], [y0+y,y0+y-y2]

      endif

    endif

  endfor

  if (unit_length gt 0) and (pos(0) gt -1e32) then begin

    x0 = pos(0)
    y0 = pos(1)

    if (unit_length eq -1.0e32) then begin
      l = (e^2.0 + n^2.0)^0.5
      unit_lenth = max(l)
    endif

    y = -unit_length*normalization
    x = 0.0

    plots, [x0,x0+x], [y0,y0+y]

    if (arrow) then begin

      length = (x^2+y^2)^0.5
      t2     = asin(y/length)
      if (x lt 0.0) then t2 = !pi - t2

      x1 = 0.5*length*cos(t2 + 15.0*!pi/180.0)
      y1 = 0.5*length*sin(t2 + 15.0*!pi/180.0)

      x2 = 0.5*length*cos(t2 - 15.0*!pi/180.0)
      y2 = 0.5*length*sin(t2 - 15.0*!pi/180.0)

      plots, [x0+x, x0+x-x1], [y0+y,y0+y-y1]
      plots, [x0+x, x0+x-x2], [y0+y,y0+y-y2]

    endif

    if strlen(title) gt 0 then 						$
      xyouts, x0+x-unit_length*normalization/5.0,y0+y/2.0, title, 	$
	alignment = 1.0, charsize = 0.8

  endif

  return

end
