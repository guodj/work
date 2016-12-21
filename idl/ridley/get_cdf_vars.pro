pro get_vars, id, att, vars, me, n, titles, location, cf, nv

  n = n + 1
  if n le me then begin
    value = '                      '
    cdf_attget, id, att, vars(n), value
    if strlen(value(0)) gt 1 then begin
      bv = byte(value(0))
      ts = string(bv)
      titles(cf,nv) = ts+'              '
      location(cf,nv,0) = n
      location(cf,nv,1) = 0
      nv = nv + 1
    endif
  endif

  return

end

pro get_units, id, att, vars, me, n, unitname, location, cf

  n = n + 1
  if n le me then begin
    value = ''
    cdf_attget, id, att, vars(n), value
    if strlen(value(0)) gt 1 then begin
      loc = where(location(cf,*,0) eq n, count)
      if count gt 0 then unitname(cf,loc) = value(0)
    endif
  endif

  return

end

pro get_ptr, id, att, vars, me, n, titles, location, cf, nv

  n = n + 1
  if n le me then begin
    value = ''
    cdf_attget, id, att, vars(n), value
    if strlen(value(0)) gt 1 then begin
      cdf_varget, id, value, vardum
      vardum = string(vardum)
      nx = n_elements(vardum(*,0))
      ny = n_elements(vardum(0,*))
      if ny eq 0 then begin
        for i = 0, nx-1 do begin
	  titles(cf,nv) = vardum(i,0)+'               '
          location(cf,nv,0) = n
          location(cf,nv,1) = i
          location(cf,nv,2) = 0
	  nv = nv + 1
        endfor
      endif else begin
        for i = 0, nx-1 do for k = 0, ny-1 do begin
	  titles(cf,nv) = vardum(i,k)+'               '
          location(cf,nv,0) = n
          location(cf,nv,1) = i
          location(cf,nv,2) = k
	  nv = nv + 1
        endfor
      endelse
    endif
  endif

  return

end

pro get_unit_ptr, id, att, vars, me, n, unitname, location, cf

  n = n + 1
  if n le me then begin
    value = ''
    cdf_attget, id, att, vars(n), value
    if strlen(value(0)) gt 1 then begin
      cdf_varget, id, value, vardum
      for i = 0, n_elements(vardum(0,*,0))-1 do  			$
	for k = 0, n_elements(vardum(0,0,*))-1 do begin
	  fakename = ''
	  for j=0, n_elements(vardum(*,0,0))-1 do			$
	    fakename = fakename+string(vardum(j,i))
          loc = where(location(cf,*,0) eq n, count)
          if count eq n_elements(vardum(0,*)) then		$
	    unitname(cf,loc(i)) = fakename
      endfor
    endif
  endif

  return

end

pro get_cdf_vars, numoffil, datafile, timeint, titles, units, 	$
	ncol, nrows, location, descript

  catch, /cancel

  nrows = intarr(numoffil)
  ncol = intarr(numoffil)
  titles = strarr(numoffil,150)
  units = strarr(numoffil,150)
  descript = strarr(numoffil,151)
  timeint = strarr(numoffil,2)
  location = intarr(numoffil, 150, 3)

  for curfile = 0,numoffil-1 do begin

    id = cdf_open(datafile(curfile))
    info = cdf_inquire(id)

    natt = info.natts
    nvar = info.nvars
    if nvar eq 0 then begin
      nvar = info.nzvars
      zvar = 1
    endif else zvar = 0

    atts = strarr(natt)
    vars = strarr(nvar)
    maxe = intarr(natt)
    nvaratt = 0

    nrows(curfile) = info.maxrec
    if (nrows(curfile) lt 1) then begin
      cdf_control, id, variable=0, /zvariable,get_var_info=v
      nrows(curfile) = v.maxrec
    endif

    for i=0,natt-1 do begin

      cdf_attinq, id, i, name, scope, me, mz

      if (strmid(scope,0,1) eq 'v') or (strmid(scope,0,1) eq 'V') then begin
        atts(nvaratt) = name
        if (zvar eq 0) then maxe(nvaratt) = me 		$
        else maxe(nvaratt) = mz
        nvaratt = nvaratt + 1
      endif

    endfor

    for i=0,nvar-1 do begin

      if (zvar eq 0) then r = cdf_varinq(id,i)		$
      else r = cdf_varinq(id,i,/zvariable)
      vars(i) = r.name

    endfor

    nv = 0

    loc = where(atts eq 'LABLAXIS',count)

    if count gt 0 then begin

      i = loc(0)
      me = maxe(i)
      n = -1
      catch, error_status
  
      while n le me do begin
	get_vars, id, atts(i), vars, me, n, titles, location, curfile, nv
      endwhile

    endif

    catch, /cancel

    loc = where(atts eq 'LABL_PTR_1',count)

    if count gt 0 then begin

      i = loc(0)
      me = maxe(i)
      n = -1

      catch, error_status
  
      while n le me do begin
	get_ptr, id, atts(i), vars, me, n, titles, location, curfile, nv
      endwhile

    endif

    catch, /cancel

    loc = where(atts eq 'UNITS',count)
    if count gt 0 then begin

      i = loc(0)
      me = maxe(i)
      n = -1
      catch, error_status
  
      while n le me do begin
	get_units, id, atts(i), vars, me, n, units, location, curfile
      endwhile

    endif

    catch, /cancel

    loc = where(atts eq 'UNIT_PTR',count)
    if count gt 0 then begin

      i = loc(0)
      me = maxe(i)
      n = -1

      catch, error_status
  
      while n le me do begin
	get_unit_ptr, id, atts(i), vars, me, n, units, location,curfile
      endwhile

    endif

    catch, /cancel

    loc = where(atts eq 'FIELDNAM',count)
    if count gt 0 then begin

      i = loc(0)
      me = maxe(i)
      n = -1
      catch, error_status
  
      while n le me do begin
	get_units, id, atts(i), vars, me, n, descript, location,curfile
      endwhile

    endif

    catch, /cancel

    ncol(curfile) = nv-1

    if (zvar eq 0) then begin
      cdf_varget1, id, 0, btime, rec_start = 0
      cdf_varget1, id, 0, etime, rec_start = nrows(curfile)
    endif else begin
      cdf_varget1, id, 0, btime, rec_start = 0, /zvariable
      cdf_varget1, id, 0, etime, rec_start = nrows(curfile), /zvariable
    endelse

    cdf_epoch, btime, by, bm, bd, bh, bmi, bs, /break
    cdf_epoch, etime, ey, em, ed, eh, emi, es, /break

    c_a_to_s, [by, bm, bd, bh, bmi, bs], bstring
    c_a_to_s, [ey, em, ed, eh, emi, es], estring
    timeint(curfile,0) = bstring
    timeint(curfile,1) = estring

    cdf_close, id

    if (zvar eq 0) then descript(curfile,ncol(curfile)+1) = 'R'		$
    else descript(curfile,ncol(curfile)+1) = 'Z'

  endfor

  max = max(ncol)

  titles = titles(*,0:max)
  units = units(*,0:max)
  descript = descript(*,0:max+1)

  for i=0,numoffil-1 do begin
    descript(i,ncol(i)+1:n_elements(descript(i,*))-1) = descript(i,ncol(i)+1)
  endfor

  return

end

