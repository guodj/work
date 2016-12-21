PRO draw_arrow, final_x, final_y, opos

  if n_elements(opos) eq 0 then opos = [0,0,1,1]
  loc = where(opos ne 0.0, count)
  if count lt 2 then opos = [0,0,1,1]
  fac = 0.5
  chi = !pi/6.0

  plot, [0.0,1.0], pos = opos, xstyle=5,ystyle=5, /nodata, /noerase
  oposd =float(!p.clip)
  oposd = oposd(0:3)

  plot, [0.0,1.0], pos = [0,0,1,1],xstyle=5,ystyle=5, /nodata, /noerase

  drange = float([!d.x_size,!d.y_size])

  out = 0
  firstout=0

  x_pos = -1
  y_pos = -1
  xlast = -1
  ylast = -1

  first = 0
  xrangeo=fltarr(1)
  yrangeo=fltarr(1)

  while not out do begin

    cursor,x,y,2
    mouse = !err

    if mouse eq 1 then begin

      if x_pos(0) eq -1 then begin
        x_pos = x
        y_pos = y
        first = 1
        xold = x*drange(0)
        yold = y*drange(1)
        pict = tvrd(xold,yold,1,1)
      endif else begin
        x_pos = [x_pos,x]
        y_pos = [y_pos,y]
        first = 2
      endelse

      while mouse ne 0 do begin

	cursor,x,y,2
        mouse = !err

      endwhile

    endif

    if (first eq 1) and (xlast ne x) and (ylast ne y) then begin

      r = sqrt((x_pos(0)-x)^2.0+(y_pos(0)-y)^2.0)
      if r gt 0.0 then begin
        theta = asin((y_pos(0)-y)/r)
        if (x_pos(0)-x) lt 0.0 then theta = !pi - theta
        x_cor = [x,x_pos(0),x+fac*r*cos(theta+chi),x+fac*r*cos(theta-chi)]
        y_cor = [y,y_pos(0),y+fac*r*sin(theta+chi),y+fac*r*sin(theta-chi)]
        tv, pict, xold, yold
	xnew = [min(x_cor)*drange(0)-1,(max(x_cor)-min(x_cor))*drange(0)+3]
	ynew = [min(y_cor)*drange(1)-1,(max(y_cor)-min(y_cor))*drange(1)+3]
        pict = tvrd(xnew(0),ynew(0),xnew(1),ynew(1))
	xold = xnew(0)
	yold = ynew(0)
        oplot, [x_pos(0),x],[y_pos(0),y]
        oplot, [x,x_cor(2)],[y,y_cor(2)]
        oplot, [x,x_cor(3)],[y,y_cor(3)]
      endif

      xlast = x
      ylast = y

      firstout = 0

    endif

    if (mouse eq 2) and (x_pos(0) ne -1) then begin

      while mouse ne 0 do begin
	cursor,x,y,2
	mouse = !err
      endwhile

      x_pos = -1
      y_pos = -1
      tv, pict, xold, yold
      firstout = 0

    endif

    if mouse eq 4 then begin

      while mouse ne 0 do begin
	cursor,x,y,2
	mouse = !err
      endwhile

      if firstout eq 1 then out = 1 else begin

        firstout = 1
	x_pos = x_pos*drange(0)
	y_pos = y_pos*drange(1)
	fx = (x_pos-oposd(0))/(oposd(2)-oposd(0))
	fy = (y_pos-oposd(1))/(oposd(3)-oposd(1))
	if n_elements(final_x) le 1 then begin
	  final_x = [fx,-1]
	  final_y = [fy,-1]
        endif else begin
	  final_x = [final_x,fx,-1]
	  final_y = [final_y,fx,-1]
        endelse

        x_pos = -1
        y_pos = -1
	print, 'line saved'

      endelse

    endif

  endwhile

  return

END

