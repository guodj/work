function loc_cdf_row, id, time, srow=srow, erow=erow, goodrow = goodrow, zvar

  if n_elements(srow) eq 0 then srow = 0
  if n_elements(erow) eq 0 then begin
    r = cdf_inquire(id)
    erow = r.maxrec
    if erow lt 0 then begin
      cdf_control, id, variable=0, /zvariable,get_var_info=v
      erow = v.maxrec
    endif
  endif

  if n_elements(goodrow) eq 0 then goodrow = -1
  if srow ge erow then goodrow = erow

  if (goodrow eq -1) then begin

    if (zvar eq 0) then begin
      cdf_varget1, id, 0, stime, rec_start=srow
      cdf_varget1, id, 0, etime, rec_start=erow
    endif else begin
      cdf_varget1, id, 0, stime, rec_start=srow, /zvariable
      cdf_varget1, id, 0, etime, rec_start=erow, /zvariable
    endelse

    if stime eq time then goodrow = srow				$
    else if etime eq time then goodrow = erow				$
    else begin

      mrow = (srow+erow)/2

      if (mrow eq srow) or (mrow eq erow) then goodrow = srow		$
      else begin

        if (zvar eq 0) then begin
	  cdf_varget1, id, 0, mtime, rec_start=mrow
        endif else begin
	  cdf_varget1, id, 0, mtime, rec_start=mrow, /zvariable
        endelse

	if time eq mtime then goodrow = mrow				$
	else begin

	  if time gt mtime then srow = mrow				$
	  else erow = mrow

	  goodrow = loc_cdf_row(id,time,srow=srow,erow=erow,		$
		goodrow=goodrow, zvar)

	endelse

      endelse

    endelse

  endif

  return, goodrow

end  

  

  
