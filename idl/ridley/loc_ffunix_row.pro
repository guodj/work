function loc_ffunix_row, id, time, srow=srow, erow=erow, goodrow = goodrow

  common ffinfo, header

  if n_elements(srow) eq 0 then srow = long(0)
  if n_elements(erow) eq 0 then begin
    erow = long(header.nr(id)-1)
  endif

  if n_elements(goodrow) eq 0 then goodrow = long(-1)
  if srow ge erow then goodrow = erow
  if (erow eq -1) then goodrow = 0

  if (goodrow eq -1) then begin

    stime = get_ffunix_time(id, srow)
    etime = get_ffunix_time(id, erow)
    while (etime lt -1.0e31) do begin
        erow = erow - 1
        etime = get_ffunix_time(id, erow)
    endwhile

    if stime eq time then goodrow = srow				$
    else if etime eq time then goodrow = erow				$
    else begin

      mrow = long((srow+erow)/2)

      if (mrow eq srow) or (mrow eq erow) then begin
          stime = get_ffunix_time(id, srow)
          etime = get_ffunix_time(id, erow)
          if (abs(stime-time) lt abs(etime-time)) then goodrow = srow else goodrow = erow
      endif else begin

	mtime = get_ffunix_time(id, mrow)

	if time eq mtime then goodrow = mrow				$
	else begin

	  if time gt mtime then srow = mrow				$
	  else erow = mrow

	  goodrow = loc_ffunix_row(id,time,srow=srow,erow=erow,		$
		goodrow=goodrow)

	endelse

      endelse

    endelse

  endif

  return, goodrow

end  

  

  
