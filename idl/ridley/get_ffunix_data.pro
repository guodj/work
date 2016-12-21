pro get_ffunix_data, single, t, dt, basetime, numofvar, numofpts, 	$
	plottype, times, datafile, grafvar, same

  common ffinfo, header

; if we keep this as general as possible, we can have up to 6 files with
; basically as many variables with as many times as we want, therefore,
; we need to figure out which variables have the same times associated
; with them, and which ones are in the same file.

; this is the getting time part :

  stime = dblarr(numofvar+1)
  etime = dblarr(numofvar+1)

  recnum = lonarr(numofvar+1,2)

  for i=0,numofvar do begin

    c_a_to_r, times(i,0,0:5), rtime
    stime(i) = rtime
    c_a_to_r, times(i,1,0:5), rtime
    etime(i) = rtime

  endfor

; number of file we have to open :

  nfiles =  max(grafvar(0:numofvar,0))

; cf - current file

  for cf = 0, nfiles do begin

; floc - which variables are associated which the current file
; fcount - how many variables are associated with that file

    floc = where(grafvar(0:numofvar,0) eq cf, fcount)

; cv  - current variable
; cvi - actual variable number

    if fcount gt 0 then for cv = 0, fcount-1 do begin

      cvi = floc(cv)

; if it's the first time through the loop, open the file and get the start
; and stop rows

      if cv eq 0 then begin

	id = cf
	header.unit(id) = 51
	openr, header.unit(id), header.df(id)+'.dat'
	recnum(cvi,0) = loc_ffunix_row(id,stime(cvi))
	recnum(cvi,1) = loc_ffunix_row(id,etime(cvi))

      endif else begin

; if it is not the first time through, check to see if we have search for
; the start and stop times on previous occasions. If we have, then we can
; just assign the already found values to the current inquires. If not,
; search.

	locv = where(stime(floc(0:cv-1)) eq stime(cvi), vcount)
	if vcount gt 0 then 						$
	  recnum(cvi,0) = recnum(floc(locv(0)),0)			$
	else recnum(cvi,0) = loc_ffunix_row(id,stime(cvi))

	locv = where(etime(floc(0:cv-1)) eq etime(cvi), vcount)
	if vcount gt 0 then 						$
	  recnum(cvi,1) = recnum(floc(locv(0)),1)			$
	else recnum(cvi,1) = loc_ffunix_row(id,etime(cvi))

      endelse

; close the current cdf file if it is the last time through the loop

      if cv eq fcount-1 then close, header.unit(id)

    endfor

  endfor

  maxne = max(recnum(*,1) - recnum(*,0))

  single = fltarr(numofvar+1, maxne+1)
  t = fltarr(numofvar+1, maxne+1)
  dt = fltarr(numofvar+1, maxne+1)
  basetime = fltarr(numofvar+1)
  numofpts = intarr(numofvar+1)

  timemin = min(stime)
  basetime(0:numofvar) = timemin

  missing = -1.0e32

  for cf = 0, nfiles do begin

    floc = where(grafvar(0:numofvar,0) eq cf, fcount)

    if fcount gt 0 then begin

      id = cf
      header.unit(id) = 51
      openr, header.unit(id), header.df(id)+'.dat' 
      ncol = (header.rl(id)-8)/4
      basedata = fltarr(ncol)

      for cr = recnum(floc(0),0), recnum(floc(0),1) do begin

        dumtime = get_ffunix_time(id,cr)
        readu, header.unit(id), basedata

        for cv = 0, fcount-1 do begin

          cvi = floc(cv)

	  dumdata = basedata(grafvar(cvi,1)-1)

	  if (dumdata ne missing) or				$
	     ((plottype ge 8) and (plottype ne 20)) or		$
	     (plottype eq 5) then begin

	    if dumdata ne missing then				$
	      single(cvi,numofpts(cvi)) = dumdata			$
	    else single(cvi,numofpts(cvi)) = -1.0e10
	    t(cvi,numofpts(cvi)) = dumtime - timemin
	    if numofpts(cvi) gt 0 then 				$
	      dt(cvi,numofpts(cvi)-1) = 				$
	        t(cvi,numofpts(cvi)) - t(cvi,numofpts(cvi)-1)
	    numofpts(cvi) = numofpts(cvi) + 1

	  endif

;          numofpts(cvi) = numofpts(cvi)-1

        endfor

      endfor

      close, header.unit(id)

    endif

  endfor

  return

end

