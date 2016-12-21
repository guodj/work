pro read_flat_scalor, stime, etime, col, time, data, nrows,	$
                      filenum = filenum,			$
                      filename = filename

  common ffinfo, header
 
  if (n_elements(filename) gt 0) then begin

    numoffil = n_elements(filename)
    datafile = strarr(6)
    datafile(0:numoffil-1) = filename
    Setup_Var_Menu, numoffil, datafile, timeint, titles, units,  $
                    count, nrows, rowlen, colloc, coltype, colnele

    ntotrows = nrows

    header = {nf : numoffil, df : datafile, ti : timeint,  $  
   	      na : titles, un : units, nr : ntotrows, rl : rowlen,   $
              loc : colloc, type : coltype, nele : colnele,   $
   	      unit : intarr(numoffil+1), ncol : count+1}

  endif

  if (n_elements(filenum) eq 0) then begin

    filenum = intarr(n_elements(col)) + 0

  endif

; make sure times are in double form and not array form:

  if (n_elements(stime) gt 1) then 			$
    c_a_to_r, stime, rstime				$
  else rstime = stime

  if (n_elements(etime) gt 1) then 			$
    c_a_to_r, etime, retime				$
  else retime = etime

; how many variables do we have :

  nvar = n_elements(col)

; figure out how many files we have:

  nfiles = max(filenum)

  id    = intarr(nvar)
  srow  = lonarr(nvar)
  erow  = lonarr(nvar)

  for file=0,nfiles do begin

    loc = where(filenum eq file, count)

    if count gt 0 then begin

      id(loc) = file
      header.unit(file) = 51+file
      openr, header.unit(file), header.df(file)+'.dat'

      srow(loc) = loc_ffunix_row(file, rstime)
      erow(loc) = loc_ffunix_row(file, retime)

    endif

  endfor

; set up data array

  nrow = erow - srow + 1

  max_nrow = max(nrow)

  data = fltarr(nvar, max_nrow)
  time = dblarr(nvar, max_nrow)

  for i=0,nvar-1 do begin

    start_row = srow(i)
    end_row   = erow(i)

    for j=start_row,end_row do begin

      fdget_scalor, id(i), j, col(i), temptime, tempdata

      time(i,j-start_row) = temptime
      data(i,j-start_row) = tempdata

    endfor

  endfor

  nrows = erow - srow + 1

  for file=0,nfiles do begin

    loc = where(filenum eq file, count)

    if count gt 0 then close, header.unit(file)

  endfor

  return

end
