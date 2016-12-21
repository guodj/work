PRO draw_line, final_x, final_y, opos

  if n_elements(opos) eq 0 then opos = [0,0,1,1]
  loc = where(opos ne 0, count)
  if count lt 2 then opos = [0,0,1,1]
  plot, [0.0,1.0], pos = opos,xstyle=5,ystyle=5, /nodata, /noerase

  drange = [!d.x_size,!d.y_size]
  oposd = float(!p.clip)
  oposd = oposd(0:3)

  plot, [0.0,drange(0)],[0.0,drange(1)], pos = [0,0,1,1] ,		$
	xstyle=5,ystyle=5, /nodata, /noerase

  out = 0
  firstout=0

  x_pos = -1
  y_pos = -1
  xlast = -1
  ylast = -1

  first = 0
  xrangeo=fltarr(1)
  yrangeo=fltarr(1)

  orient = 0
  xfir = 5
  yfir = 5
  xsec = 105
  ysec = 5

  a = findgen(17)*!pi/8
  usersym, 1.0*sin(a), 1.0*cos(a), /fill
  x1r = xfir-5
  y1r = yfir-5
  x2r = xsec+5
  y2r = ysec+5
  if x1r lt 0.0 then x1r = 0.0
  if y1r lt 0.0 then y1r = 0.0
  if x2r gt drange(0) then x2r = drange(0)
  if y2r gt drange(1) then y2r = drange(1)
  pict = tvrd(x1r,y1r,x2r-x1r,y2r-y1r)

  oplot, [xfir, xsec], [yfir, ysec]
  oplot, [xsec], [ysec], psym=8
  oplot, [xfir], [yfir], psym=8

  dx = xsec-xfir
  dy = ysec-yfir
  xfiro = x1r
  yfiro = y1r


  while not out do begin

    cursor,x,y,2, /device
    mouse = !err

    if mouse eq 1 then begin

      first = 1
      d1 = sqrt((xfir-x)^2.0+(yfir-y)^2.0)
      d2 = sqrt((xsec-x)^2.0+(ysec-y)^2.0)
      len = sqrt(dx^2.0+dy^2.0)
      if (d1 le len/5.) and (d2 gt len/5.0) then begin
	xfir = x
	yfir = y
      endif else begin
        if (d1 gt len/5.0) and (d2 le len/5.0) then begin
	  xsec = x
	  ysec = y
        endif else begin
	  xfir = x - dx/2.0
	  yfir = y - dy/2.0
	  xsec = x + dx/2.0
	  ysec = y + dy/2.0
        endelse
      endelse

      dx = xsec-xfir
      dy = ysec-yfir

    endif

    if mouse eq 2 then begin

      xfir = 5
      yfir = 5
      xsec = 105
      ysec = 5

    endif
    if (xlast ne x) and (ylast ne y) and 				      $
       ((mouse ne 0) or (first eq 0)) then begin

      if first eq 0 then begin

	xfir = x - dx/2.0
	yfir = y - dy/2.0
	xsec = x + dx/2.0
	ysec = y - dy/2.0

      endif

      tv, pict, xfiro, yfiro

      x1r = xfir-5
      y1r = yfir-5
      x2r = xsec+5
      y2r = ysec+5
      if x2r le x1r+7 then begin
	dum = x1r
	x1r = x2r-12
	x2r = dum+12
      endif
      if y2r le y1r+7 then begin
	dum = y1r
	y1r = y2r-12
	y2r = dum+12
      endif
      if x1r lt 0.0 then x1r = 0.0
      if y1r lt 0.0 then y1r = 0.0
      if x2r gt drange(0) then x2r = drange(0)
      if y2r gt drange(1) then y2r = drange(1)
      pict = tvrd(x1r,y1r,x2r-x1r,y2r-y1r)

      oplot, [xfir, xsec], [yfir, ysec]
      oplot, [xsec], [ysec], psym=8
      oplot, [xfir], [yfir], psym=8

      xfiro = x1r
      yfiro = y1r

    endif

    if mouse eq 4 then begin

      while mouse ne 0 do begin
	cursor,x,y,2, /device
	mouse = !err
      endwhile

      out = 1

      fx1 = (float(xfir)-oposd(0))/(oposd(2)-oposd(0))
      fy1 = (float(yfir)-oposd(1))/(oposd(3)-oposd(1))
      fx2 = (float(xsec)-oposd(0))/(oposd(2)-oposd(0))
      fy2 = (float(ysec)-oposd(1))/(oposd(3)-oposd(1))

      if final_x(0) eq -1 then begin
	final_x  = [fx1,fx2,-1]
	final_y  = [fy1,fy2,-1]
      endif else begin
	final_x  = [final_x,fx1,fx2,-1]
	final_y  = [final_y,fy1,fy2,-1]
      endelse

    endif

    if mouse eq 3 then begin

      out = 1
      tv, pict, xfiro, yfiro

    endif

  endwhile

  !p.font = 0

END

