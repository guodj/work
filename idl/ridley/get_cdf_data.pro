pro get_cdf_data, single, t, dt, basetime, numofvar, numofpts, 		$
	plottype, times, datafile, grafvar, same, descript

; if we keep this as general as possible, we can have up to 6 files with
; basically as many variables with as many times as we want, therefore,
; we need to figure out which variables have the same times associated
; with them, and which ones are in the same file.

; this is the getting time part :

  stime = dblarr(numofvar+1)
  etime = dblarr(numofvar+1)

  recnum = intarr(numofvar+1,2)

  for i=0,numofvar do begin

    sy  = times(i,0,0)+1900
    sm  = times(i,0,1)
    sd  = times(i,0,2)
    sh  = times(i,0,3)
    smi = times(i,0,4)
    ss  = times(i,0,5)
    ey  = times(i,1,0)+1900
    em  = times(i,1,1)
    ed  = times(i,1,2)
    eh  = times(i,1,3)
    emi = times(i,1,4)
    es  = times(i,1,5)

    cdf_epoch, dtime, sy, sm, sd, sh, smi, ss, /compute_epoch
    stime(i) = dtime
    cdf_epoch, dtime, ey, em, ed, eh, emi, es, /compute_epoch
    etime(i) = dtime

  endfor

; number of file we have to open :

  nfiles =  max(grafvar(0:numofvar,0))

; cf - current file

  for cf = 0, nfiles do begin

    if descript(cf,n_elements(descript(0,*))-1) eq 'Z' then zvar = 1	$
    else zvar = 0

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

        id = cdf_open(datafile(cf))
	recnum(cvi,0) = loc_cdf_row(id,stime(cvi), zvar)
	recnum(cvi,1) = loc_cdf_row(id,etime(cvi), zvar)
		
      endif else begin

; if it is not the first time through, check to see if we have search for
; the start and stop times on previous occasions. If we have, then we can
; just assign the already found values to the current inquires. If not,
; search.

	locv = where(stime(floc(0:cv-1)) eq stime(cvi), vcount)
	if vcount gt 0 then 						$
	  recnum(cvi,0) = recnum(floc(locv(0)),0)			$
	else recnum(cvi,0) = loc_cdf_row(id,stime(cvi), zvar)

	locv = where(etime(floc(0:cv-1)) eq etime(cvi), vcount)
	if vcount gt 0 then 						$
	  recnum(cvi,1) = recnum(floc(locv(0)),1)			$
	else recnum(cvi,1) = loc_cdf_row(id,etime(cvi), zvar)

      endelse

; close the current cdf file if it is the last time through the loop

      if cv eq fcount-1 then cdf_close, id

    endfor

  endfor

  maxne = max(recnum(*,1) - recnum(*,0))

  single = fltarr(numofvar+1, maxne+1)
  t = fltarr(numofvar+1, maxne+1)
  dt = fltarr(numofvar+1, maxne+1)
  basetime = fltarr(numofvar+1)

  timemin = min(stime)
  cdf_epoch, timemin, sy, sm, sd, sh, smi, ss, /break
  itime = [sy, sm, sd, sh, smi, ss]
  c_a_to_r, itime, tmin
  basetime(0:numofvar) = tmin

  for cf = 0, nfiles do begin

    if descript(cf,n_elements(descript(0,*))-1) eq 'Z' then zvar = 1	$
    else zvar = 0

    floc = where(grafvar(0:numofvar,0) eq cf, fcount)

    if fcount gt 0 then for cv = 0, fcount-1 do begin

      cvi = floc(cv)

      if cv eq 0 then id = cdf_open(datafile(cf))

      missing = get_cdf_fill(id, grafvar(cvi,1), zvar)

      for cr = recnum(cvi,0), recnum(cvi,1) do begin

        if (zvar eq 0) then begin
	  cdf_varget1, id, 0, dumtime, rec_start = cr
	  cdf_varget, id, grafvar(cvi,1), dumdata, rec_start = cr
        endif else begin
	  cdf_varget1, id, 0, dumtime, rec_start = cr, /zvariable
	  cdf_varget, id, grafvar(cvi,1), dumdata, rec_start = cr, /zvariable
        endelse

	dumdata = dumdata(grafvar(cvi,2), grafvar(cvi,3))

	if (dumdata ne missing) or				$
	   ((plottype ge 8) and (plottype ne 20)) or		$
	   (plottype eq 5) then begin

	  if dumdata ne missing then				$
	    single(cvi,numofpts(cvi)) = dumdata			$
	  else single(cvi,numofpts(cvi)) = -1.0e10
	  t(cvi,numofpts(cvi)) = (dumtime - timemin)/1000.0
	  if numofpts(cvi) gt 0 then 				$
	    dt(cvi,numofpts(cvi)) = 				$
	      t(cvi,numofpts(cvi)) - t(cvi,numofpts(cvi)-1)
	  numofpts(cvi) = numofpts(cvi) + 1

	endif

      endfor

    endfor

  endfor

  return

end

