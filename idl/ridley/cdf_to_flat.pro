pro cdf_to_flat, listin, fileout

  flatmissing = -1.0e32

  if n_elements(listin) eq 0 then begin

    genlist = ''
    print, 'Enter file names of cdf files (ex: *.cdf) :'
    read, genlist

    listin = findfile(genlist, count=count)

  endif else count = n_elements(listin)

  if n_elements(fileout) eq 0 then begin

    fileout=''
    print, 'Enter flatfile to output :'
    read, fileout

  endif

  if count gt 0 then begin

    get_cdf_vars, 1, listin(0), timeint, variables, units,		$
	nvars, nrows, location, description

    fdnfile = 0.0
    fdids = fltarr(6,6)

    unit = ffcreate(fileout,variables,units,description,fdnfile,fdids)

    data = fltarr(100)
    uttime = double(0.0)
    itime = intarr(6)
    outrow = 1.0

    for i=0,count-1 do begin

      print, 'Working on CDF file : ',listin(i)

      get_cdf_vars, 1, listin(i), timeint, variables, units,		$
	nvars, nrows, location, description

      id = cdf_open(listin(i))

      cdfmissing = fltarr(nvars(0)+1)

      for j=1,nvars(0)+1 do cdfmissing(j-1) = get_cdf_fill(id,location(0,j,0))

      for j=0,nrows(0) do begin

	cdf_varget1, id, 0, dumtime, rec_start = j

        for k=1,nvars(0)+1 do begin

	  cdf_varget, id, location(0,k,0), dumdata, rec_start = j
          if dumdata(location(0,k,1),location(0,k,2)) ne cdfmissing(k-1) then $
		data(k+2)=dumdata(location(0,k,1),location(0,k,2))	      $
	  else data(k+2)=flatmissing

        endfor

        cdf_epoch, dumtime, y, m, d, h, mi, s, /break
        itime(0) = y
        itime(1) = m
        itime(2) = d
        itime(3) = h
        itime(4) = mi
        itime(5) = s

        c_a_to_r, itime, uttime

        z = call_external('fput_exe','fput',				$
		fdnfile,fdids,unit,data,uttime,outrow)
        outrow=outrow+1.0

      endfor

      cdf_close, id

    endfor

    z = call_external('fclose_exe','fclose',fdnfile,fdids,unit)

  endif

  return

end
