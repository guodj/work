PRO str_put, final_x, final_y, char_s, or_list, str_list, opos

  if n_elements(opos) eq 0 then opos = [0,0,1,1]
  loc = where(opos ne 0.0, count)
  if count lt 2 then opos = [0,0,1,1]

  plot, [0.0,1.0], pos = opos, xstyle=5 ,ystyle=5, /nodata, /noerase

  if final_x(0) ne -1 then begin

    for i=0,n_elements(final_x)-1 do begin
     
      xyouts,final_x(i),final_y(i),str_list(i),				      $
	     charsize=char_s(i), orient=or_list(i)

    endfor

  endif

END

