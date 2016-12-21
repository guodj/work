PRO draw_final, final_x, final_y, opos

  if n_elements(opos) eq 0 then opos = [0,0,1,1]
  loc = where(opos ne 0.0, count)
  if count lt 2 then opos = [0,0,1,1]

  plot, [0.0,1.0], pos = opos,xstyle=5,ystyle=5, /nodata, /noerase

  if n_elements(final_x) gt 1 then begin

    list = where(final_x eq -1, count)
    list = [-1,list]

    for ii=0,count-1 do begin

      x_pos = final_x(list(ii)+1:list(ii+1)-1)
      y_pos = final_y(list(ii)+1:list(ii+1)-1)

      nele = n_elements(x_pos)

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

      if nele eq 1 then begin
	oplot, [x_pos], [y_pos], psym = 3
      endif

      if nele eq 2 then begin
	oplot, x_pos, y_pos, thick = 6
      endif

      if nele ge 3 then begin
        oplot, cubic_x, cubic_y, thick = 6
      endif

    endfor

  endif

END

