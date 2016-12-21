
pro get_var_and_pos, grafvar, titlenam, unitnam, times, screen

  common timeblock, sid, time, bftime, eftime
  common varblock, nfile, lid, titles, units, lfn
  COMMON type, plottype, types
  common extras, unconpts, columns, symtype, maxminque, rangepl, numofrang,   $
	         titleofplot, plotwhere, psfile, placement, numofgra,	      $
		 psfont, offdiv, offmul, maxran, checkinvar, filename, nb,    $
		 publish, checkstat, statfile, xyzpos, stat, lat, lon, alpha, $
		 mln, remmean, velocity, beam, step, crb, clr, draw,	      $
		 final_x, final_y, arrow_x, arrow_y, arrow, str_x, str_y,     $
		 str_s, str_o,str_l, str, checkselect, percent, inarow, ext
  common rep_data, change
  common data, single, t, dt, basetime, numofvar, numofpts

  grafvar = intarr(100,4)

  sel = where(strmid(titles,15,8) eq 'selected', count)
  j = sel/nfile
  i = sel mod nfile

  if count gt 0 then begin

    titlenam = strmid(titles(sel),0,14)
    unitnam = units(sel)

;  type_10
;  type_20
;  type_21

    if (plottype ne 10) and 						      $
       (plottype ne 20) and						      $
       (plottype ne 21) then begin

      numofvar = count - 1
      grafvar = intarr(count,4)
      grafvar(0:count-1,0) = i
      for k = 0,count-1 do begin
        if change.flat(i(k)) eq 1 then grafvar(0:count-1,1) = j+1
        if change.cdf(i(k)) eq 1 then begin
	  grafvar(k,1) = change.location(i(k),j(k),0)
	  grafvar(k,2) = change.location(i(k),j(k),1)
	  grafvar(k,3) = change.location(i(k),j(k),2)
        endif
      endfor

    endif else begin

      numofvar = 3*count - 1
      x = indgen(count)
      grafvar = intarr(count*3,4)
      grafvar(3*x,0) = i
      grafvar(3*x,1) = xyzpos(j,0)
      grafvar(3*x+1,0) = i 
      grafvar(3*x+1,1) = xyzpos(j,1)
      grafvar(3*x+2,0) = i 
      grafvar(3*x+2,1) = xyzpos(j,2)

    endelse

  endif else begin

    numofvar = -1
    titlenam = strarr(1)
    unitnam = strarr(1)

  endelse

  if numofvar ge 0 then begin

    grafvar = grafvar(0:numofvar,*)

    if change.place eq 1 then begin

      placement = intarr(numofvar+3, numofvar+3)

      Get_Positions, plottype, placement, numofgra, numofvar, titlenam, screen

      change.place = 0

    endif

    times = intarr(numofvar+1,2,6)
    times(0,0,0:5) = time(0:5)
    times(0,1,0:5) = time(6:11)
    times(0,0,0) = times(0,0,0)-1900
    times(0,1,0) = times(0,1,0)-1900
    for j=1, numofvar do times(j,*,*) = times(0,*,*)

    if checkinvar eq 2 then begin

      if (change.data eq 1) and (change.invar eq 1) then begin

	if n_elements(screen) eq 0 then begin
          datafile = strarr(1)
          numoffil = 0
          maxnumfil = 1
          filtitle = 'Get File For Invarient Latitudes'
	  initpath = ''
          exten = 'HED'
          get_files, datafile, numoffil, maxnumfil, filtitle, initpath, exten
	endif else begin
	  datafile = ''
	  print, 'Enter File For Invarient Latitudes'
	  read, datafile
	endelse

        filename(nfile) = datafile(0)
	change.invar = 0

      endif

      numofvar = numofvar + 3
      grafvar = [grafvar,intarr(3,4)]
      times = [times,intarr(3,2,6)]

      for curvar=numofvar-2,numofvar do begin

        grafvar(curvar,0) = nfile
        grafvar(curvar,1) = curvar - (numofvar-2) + 1

        for i = 0,1 do begin

	  times(curvar,i,0) = times(0,i,0)
	  times(curvar,i,1) = 1
	  times(curvar,i,2) = 1
	  times(curvar,i,3) = 1
	  if i eq 0 then times(curvar,i,4) = 0 else times(curvar,i,4) = 10
	  times(curvar,i,5) = 0

	endfor

      endfor

    endif

  endif

  RETURN

END
