PRO draw_arr_final, final_x, final_y, opos, line

  if n_elements(line) eq 0 then line = 0 else line = 1
  if n_elements(opos) eq 0 then opos = [0,0,1,1]
  loc = where(opos ne 0.0, count)
  if count lt 2 then opos = [0,0,1,1]

  if (line eq 1) and (n_elements(final_x) gt 1) then begin
    loc = where(final_x eq -1, count)
    if count gt 0 then begin
      final_x = final_x(0:loc(count-1))
      final_y = final_y(0:loc(count-1))
    endif
  endif

  fac = 0.5
  chi = !pi/6.0

  plot, [0.0,1.0], pos = opos,xstyle=5,ystyle=5, /nodata, /noerase

  if n_elements(final_x) gt 1 then begin

    list = where(final_x eq -1, count)
    list = [-1,list]

    for ii=0,count-1 do begin

      x_pos = final_x(list(ii)+1:list(ii+1)-1)
      y_pos = final_y(list(ii)+1:list(ii+1)-1)

      r = sqrt((x_pos(0)-x_pos(1))^2.0+(y_pos(0)-y_pos(1))^2.0)

      if r gt 0.0 then begin

	if line eq 0 then begin
          theta = asin((y_pos(0)-y_pos(1))/r)
          if (x_pos(0)-x_pos(1)) lt 0.0 then theta = !pi - theta
          x_cor = [x_pos(1),x_pos(0),x_pos(1)+fac*r*cos(theta+chi),	      $
		   x_pos(1)+fac*r*cos(theta-chi)]
          y_cor = [y_pos(1),y_pos(0),y_pos(1)+fac*r*sin(theta+chi),	      $
		   y_pos(1)+fac*r*sin(theta-chi)]
          oplot, [x_cor(0),x_cor(2)],[y_cor(0),y_cor(2)]
          oplot, [x_cor(0),x_cor(3)],[y_cor(0),y_cor(3)]
	endif

        oplot, [x_pos(0),x_pos(1)],[y_pos(0),y_pos(1)]

      endif

    endfor

  endif

END

