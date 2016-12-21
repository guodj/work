
PRO draw_curve, final_x, final_y, opos

  if n_elements(opos) eq 0 then opos = [0,0,1,1]
  loc = where(opos ne 0.0, count)
  if count lt 2 then opos = [0,0,1,1]

  plot, [0.0,1.0], pos = opos,xstyle=5,ystyle=5, /nodata, /noerase

  yposr = opos(3) - opos(1)
  xposr = opos(2) - opos(0)
  ybot = opos(1)
  xbot = opos(0)

  drange = [!d.x_size,!d.y_size]

  out = 0
  firstout=0

  x_pos = -1
  y_pos = -1

  first = 0
  xrangeo=fltarr(1)
  yrangeo=fltarr(1)

  while not out do begin

    cursor,x,y,2
    mouse = !err

    letoff = 0

    while mouse eq 1 do begin

      if letoff eq 0 then begin

        if x_pos(0) eq -1 then begin
          x_pos = x
          y_pos = y
        endif else begin
          x_pos = [x_pos,x]
          y_pos = [y_pos,y]
	endelse

        x_old = x_pos
        y_old = y_pos

      endif else begin

	cursor,x,y,2
        mouse = !err

	x_pos(nele-1) = x
	y_pos(nele-1) = y

      endelse

      nele = n_elements(x_pos)

      if ((x_old(nele-1) ne x_pos(nele-1)) or				      $ 
         (y_old(nele-1) ne y_pos(nele-1)) or				      $
	 (letoff eq 0)) and (mouse ne 0) then begin

        x_old = x_pos
        y_old = y_pos

        if nele ge 3 then begin
          sort, x_pos, order
	  dummy = y_pos
          for i=0,nele-1 do dummy(order(i)) = y_pos(i)
	  y_pos = dummy
        endif

        if nele ge 3 then begin
          natural_cubic, x_pos, cubic_x
          natural_cubic, y_pos, cubic_y
        endif

        if (nele ge 3) and (letoff eq 1) then begin
	  xa = [x_pos, cubic_x]
	  ya = [y_pos, cubic_y]
        endif else begin
	  xa = x_pos
	  ya = y_pos
        endelse

	xmin = xposr*min([xa])+xbot
	xmax = xposr*max([xa])+xbot
	ymin = yposr*min([ya])+ybot
	ymax = yposr*max([ya])+ybot
	
	xrange = [xmin*drange(0),(xmax-xmin)*drange(0)+1]		
	yrange = [ymin*drange(1),(ymax-ymin)*drange(1)+1]

	if xrange(0) lt 0 then begin
	  xrange(1) = xrange(1) - xrange(0)
	  xrange(0) = 0
	endif
	if yrange(0) lt 0 then begin
	  yrange(1) = yrange(1) - yrange(0)
	  yrange(0) = 0
	endif
	if xrange(0)+xrange(1) ge !d.x_size then			      $ 
	  xrange(1) = !d.x_size - xrange(0)-1
	if yrange(0)+yrange(1) gt !d.y_size then			      $ 
	  yrange(1) = !d.y_size - yrange(0)-1

        if first eq 1 then tv, pict, xrangeo(0), yrangeo(0) else first = 1
        pict = tvrd(xrange(0),yrange(0),xrange(1),yrange(1))

        if nele eq 1 then begin
	  oplot, [x_pos], [y_pos], psym = 3
        endif

        if nele eq 2 then begin
	  oplot, x_pos, y_pos
        endif

        if nele ge 3 then begin
          oplot, cubic_x, cubic_y
        endif

        x_pos = x_old
        y_pos = y_old
        xrangeo=xrange
        yrangeo=yrange

      endif

      letoff = 1

      firstout = 0

    endwhile

    if (mouse eq 2) and (x_pos(0) ne -1) then begin

      firstout = 0

      while mouse ne 0 do begin
	cursor,x,y,2
	mouse = !err
      endwhile

      if n_elements(x_pos) eq 1 then begin
	x_pos = -1
	y_pos = -1
        tv, pict, xrangeo(0), yrangeo(0)
      endif else begin
	x_pos = x_pos(0:n_elements(x_pos)-2)
	y_pos = y_pos(0:n_elements(y_pos)-2)

        nele = n_elements(x_pos)

        x_old = x_pos
        y_old = y_pos

        if nele ge 3 then begin
          sort, x_pos, order
	  dummy = y_pos
          for i=0,nele-1 do dummy(order(i)) = y_pos(i)
	  y_pos = dummy
        endif

        if nele ge 3 then begin
          natural_cubic, x_pos, cubic_x
          natural_cubic, y_pos, cubic_y
        endif

        if (nele ge 3) and (letoff eq 1) then begin
	  xa = [x_pos, cubic_x]
	  ya = [y_pos, cubic_y]
        endif else begin
	  xa = x_pos
	  ya = y_pos
        endelse

	xmin = xposr*min([xa])+xbot
	xmax = xposr*max([xa])+xbot
	ymin = yposr*min([ya])+ybot
	ymax = yposr*max([ya])+ybot
	
	xrange = [xmin*drange(0),(xmax-xmin)*drange(0)+1]		
	yrange = [ymin*drange(1),(ymax-ymin)*drange(1)+1]

	if xrange(0) lt 0 then begin
	  xrange(1) = xrange(1) - xrange(0)
	  xrange(0) = 0
	endif
	if yrange(0) lt 0 then begin
	  yrange(1) = yrange(1) - yrange(0)
	  yrange(0) = 0
	endif
	if xrange(0)+xrange(1) ge !d.x_size then			      $ 
	  xrange(1) = !d.x_size - xrange(0)-1
	if yrange(0)+yrange(1) gt !d.y_size then			      $ 
	  yrange(1) = !d.y_size - yrange(0)-1

        tv, pict, xrangeo(0), yrangeo(0)
        pict = tvrd(xrange(0),yrange(0),xrange(1),yrange(1))

        if nele eq 1 then begin
	  oplot, [x_pos], [y_pos], psym = 3
        endif

        if nele eq 2 then begin
	  oplot, x_pos, y_pos
        endif

        if nele ge 3 then begin
          oplot, cubic_x, cubic_y
        endif

        x_pos = x_old
        y_pos = y_old
        xrangeo=xrange
        yrangeo=yrange

      endelse

    endif

    if mouse eq 4 then begin

      while mouse ne 0 do begin
	cursor,x,y,2
	mouse = !err
      endwhile

      if firstout eq 1 then begin

	out = 1 
	print, 'exiting line drawer'

      endif else begin

        firstout = 1
	if n_elements(final_x) eq 0 then begin
	  final_x = [x_pos,-1]
	  final_y = [y_pos,-1]
        endif else begin
	  final_x = [final_x,x_pos,-1]
	  final_y = [final_y,y_pos,-1]
        endelse

        x_pos = -1
        y_pos = -1
	print, 'line saved'

      endelse

    endif

  endwhile

END

