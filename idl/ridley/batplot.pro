
;MAIN PROCEDURE


pro batplot

  common fileblock, starttime, endtime, varcount, ftime 
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
  common flag_structure, flags
  common rep_data, change
  common data, single, t, dt, basetime, numofvar, numofpts
  common text_id, intext, inps, instring, inmaxx, inmaxy, inminx, inminy
  common months, mon
  common control, getplottype, getfiles, stoppro
  common comdis, displays
  common descriptions, descript
  common ffinfo, header

  mon='JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'

  getplottype = 1
  getfiles = 1
  stoppro = 0

  while stoppro eq 0 do begin
  
    initialize

    if getplottype eq 1 then get_plot_type, /screen

    if getfiles eq 1 then begin

      maxnumfil = 6
      info = get_file_batch(maxnumfil)
      filename = info.files 
      numoffil = info.nf

    endif

    for i=0,numoffil-1 do begin

      flist = findfile(filename(i)+'.hed', count = nf)
      if nf ne 0 then change.flat(i) = 1			$
      else begin
	clist = findfile(filename(i)+'.cdf',count=nc)
	if nc ne 0 then change.cdf(i) = 1
      endelse

    endfor

    getplottype = 0
    getfiles = 0
    stopplot = 1

    if numoffil gt 0 then begin

      timeint=strarr(numoffil,2)
      titles=strarr(numoffil,150)
      units=strarr(numoffil,150)
      varcount=intarr(numoffil)

;  We want to get all of the variable names out of the flat file header
;  and the units and the start and end time.
;  titles = variable names
;  units = units
;  timeint = start time and end time

      if change.flat(0) eq 1 then begin
        setup_var_menu, numoffil, filename, timeint, titles, 		$
	  units, varcount, ntotrows, rowlen
	header = {nf : numoffil, df : filename, ti : timeint,		$
	  na : titles, un : units, nr : ntotrows, rl : rowlen, 		$
	  unit : intarr(numoffil+1)}
      endif else begin
        get_cdf_vars, numoffil, filename, timeint, titles,              $
          units, varcount, nrows, location, descript
        for i=0,numoffil-1 do for j=0,varcount(i) do                    $
          change.location(i,j,*) = location(i,j,*)
      endelse
      nfile = numoffil

;  type_10
;  type_20
;  type_21

      if (plottype eq 10) or						      $
         (plottype eq 20) or						      $
         (plottype eq 21) then begin

	print, 'In order to continue, I need to get a station file name.'
        print, 'This file should contain station names and locations,'
        print, 'along with the position of each of the columns in the'
        print, 'flat file.'
        maxnumstat = 1
	info = get_file_batch(maxnumstat)
	statfile = info.files(0)
        get_station, statfile, xyzpos, stat, lat, lon, alpha, mln

      endif

;  type_10
;  type_20
;  type_21

      if (plottype eq 10) or						      $
         (plottype eq 20) or						      $
         (plottype eq 21) then begin

        padding = '               '

        titles = strarr(1,n_elements(stat))
        titles(0,*) = stat + padding
        titles(0,*) = strmid(titles(0,*),0,14)+'(selected)'
        units  = strarr(1,n_elements(titles))+ units(0,0)
        varcount = n_elements(stat)-1
        numofvar = varcount

      endif

      convert_times, timeint, starttime, endtime, ftime, 		$
		     bftime, eftime, time, dt

      d_vars_bat, filename
      d_menu_bat

      bat_manager

    endif else stoppro = 1

  endwhile

  while !D.window gt -1 do wdelete

  get_temp_ps, dummy, broadps

  x = findfile(broadps)
  xl = strlen(x)
  if xl(0) gt 0 then spawn, 'delete '+broadps, list

end
