PRO str_get, final_x, final_y, char_s, or_list, str_list, string, opos

  if n_elements(opos) eq 0 then opos = [0,0,1,1]
  loc = where(opos ne 0.0, count)
  if count lt 2 then opos = [0,0,1,1]  

  plot, [0.0,1.0], pos = opos,xstyle=5,ystyle=5, /nodata, /noerase

  !p.font = -1

  drange = [!d.x_size,!d.y_size]
  oposd = float(!p.clip)
  oposd = oposd(0:3)

  out = 0
  firstout=0

  x_pos = -1
  y_pos = -1
  xlast = -1
  ylast = -1

  first = 0
  xrangeo=fltarr(1)
  yrangeo=fltarr(1)
  xchar = !d.x_ch_size*strlen(string)
  ychar = !d.y_ch_size
  cs = 1.0

  orient = 0
  xold = 0
  yold = 0
  pict = tvrd(xold,yold,xchar+2,ychar+2)
  xyouts, 0, 0, string, orient = orient,/device

  while not out do begin

    cursor,x,y,2, /device
    mouse = !err

    if mouse eq 1 then begin

      if x_pos eq -1 then begin
        x_pos = x
        y_pos = y
        first = 1
      endif else begin

	pix = abs(float(x_pos-x))/strlen(string)
        xchar = pix*strlen(string)*2.0
        ychar = pix*3.5
	cs = pix/3.0

      endelse

    endif

    if (mouse eq 2) and (x_pos(0) ne -1) then begin

      orient = asin((y-y_pos)/sqrt((x-x_pos)^2.0+(y-y_pos)^2.0))*180.0/!pi
      if x-x_pos lt 0.0 then orient = 180.0-orient

    endif

    if (xlast ne x) and (ylast ne y) and 				      $
       ((mouse ne 0) or (first eq 0)) then begin

        if first eq 1 then begin
	  x = x_pos 
	  y = y_pos
        endif

        tv, pict, xold, yold
        if xchar ne 0.0 then 						      $
	  ang = orient*!pi/180.0 + atan(float(ychar)/float(xchar))	      $
	else								      $
	  ang = orient*!pi/180.0 + !pi/2.0

        ang2 = orient*!pi/180.0
        x_cor = [x,x+sqrt(xchar^2.0+ychar^2.0)*cos(ang),	      $
		 x-ychar*sin(ang2), x+xchar*cos(ang2), x+ychar*sin(ang2)]
        y_cor = [y,y+sqrt(xchar^2.0+ychar^2.0)*sin(ang),	      $
		 y+xchar*sin(ang2), y+ychar*cos(ang2), y-ychar*cos(ang2)]
	xnew = [min(x_cor)-1,(max(x_cor)-min(x_cor))+3]
	ynew = [min(y_cor)-1,(max(y_cor)-min(y_cor))+3]
        if xnew(0)+xnew(1) gt drange(0) then xnew(1) = drange(0)-xnew(0)
        if ynew(0)+ynew(1) gt drange(1) then ynew(1) = drange(1)-ynew(0)
        if xnew(0) lt 0 then begin
	  xnew(1) = xnew(1)-xnew(0)
	  xnew(0) = 0
	endif
        if ynew(0) lt 0 then begin
	  ynew(1) = ynew(1)-ynew(0)
	  ynew(0) = 0
	endif
        pict = tvrd(xnew(0),ynew(0),xnew(1),ynew(1))
        xyouts, x, y, string, orient = orient, charsize = cs, /device
	xold = xnew(0)
	yold = ynew(0)

        xlast = x
        ylast = y

    endif

    if mouse eq 3 then begin

      tv, pict, xold, yold

      while mouse ne 0 do begin
	cursor,x,y,2, /device
	mouse = !err
      endwhile

      out = 1

    endif

    if mouse eq 4 then begin

      while mouse ne 0 do begin
	cursor,x,y,2, /device
	mouse = !err
      endwhile

      out = 1

      fx = (float(x_pos)-oposd(0))/(oposd(2)-oposd(0))
      fy = (float(y_pos)-oposd(1))/(oposd(3)-oposd(1))

      if final_x(0) eq -1 then begin
	final_x  = [fx]
	final_y  = [fy]
	char_s   = [cs/2.0]
	or_list  = [orient]
	str_list = [string]
      endif else begin
	final_x  = [final_x,fx]
	final_y  = [final_y,fy]
	char_s   = [char_s,cs/2.0]
	or_list  = [or_list,orient]
	str_list = [str_list,string]
      endelse

    endif

  endwhile

  !p.font = 0

END

