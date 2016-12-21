
;MAIN PROCEDURE


pro widplot

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
  common ffinfo, header
  common descriptions, descript

  mon='JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'

  getplottype = 1
  getfiles = 1
  stoppro = 0

  while stoppro eq 0 do begin
  
    initialize

    if getplottype eq 1 then get_plot_type

    if getfiles eq 1 then begin

      maxnumfil = 6
      filtitle = 'Choose Flatfiles for Plotting'
      if !version.os eq 'vms' then 					$
	if n_elements(datafile) gt 0 then 				$
	  if strlen(datafile(0)) gt 0 then 				$
	    initpath = strmid(datafile(0),0,strpos(datafile(0),']')+1)	$
	  else initpath = ''						$
	else initpath = ''						$
      else initpath = ''
      if !version.os eq 'vms' then exten = 'HED CDF'			$
      else exten = 'hed cdf'
      get_files, datafile, numoffil, maxnumfil, filtitle, initpath, exten

    endif

    for i=0,numoffil-1 do begin

      flist = findfile(datafile(i)+'.hed', count=nf)
      if nf ne 0 then change.flat(i) = 1  			$
      else begin
	clist = findfile(datafile(i)+'.cdf', count=nc)
	if nc ne 0 then change.cdf(i) = 1
      endelse

    endfor

    getplottype = 0
    getfiles = 0
    stopplot = 1

    if numoffil gt 0 then begin

      print, 'filename(s) : ', datafile(0:numoffil-1)

      timeint=strarr(numoffil,2)
      titles=strarr(numoffil,500)
      units=strarr(numoffil,500)
      varcount=intarr(numoffil)

;  type_10
;  type_20
;  type_21

      if (plottype eq 10) or						      $
         (plottype eq 20) or						      $
         (plottype eq 21) then begin

        statfile = strarr(1)
        numofstat = 0
        maxnumstat = 1
        filtitle = 'Get Station file'
	if (!version.os eq 'VMS') then begin
          initpath = 'zasu$dkb500:[software.plotting.idl.data_files]'
          exten = 'DAT'
	endif else begin
          initpath = '/home/ridley/d.idl/d.datafiles/'
          exten = 'dat'
	endelse

        get_files, statfile, numofstat, maxnumstat, filtitle, initpath, exten
        get_station, statfile, xyzpos, stat, lat, lon, alpha, mln

      endif

      if change.flat(0) eq 1 then begin
        setup_var_menu, numoffil, datafile, timeint, titles, 		$
	  units, varcount, ntotrows, rowlen, colloc, coltype, colnele
	header = {nf : numoffil, df : datafile, ti : timeint,		$
	  na : titles, un : units, nr : ntotrows, rl : rowlen, 		$
          loc : colloc, type : coltype, nele : colnele,			$
	  unit : intarr(numoffil+1)}
      endif else begin
	get_cdf_vars, numoffil, datafile, timeint, titles, 		$
	  units, varcount, nrows, location, descript
	for i=0,numoffil-1 do for j=0,varcount(i) do			$
	  change.location(i,j,*) = location(i,j,*)
      endelse
      nfile = numoffil

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

      filename = datafile

      set_up_main, base, datafile

      widget_control, base, /realize 

      XMANAGER, 'widplot', base, event_handler='widplote'

    endif else stoppro = 1

  endwhile

  while !D.window gt -1 do wdelete

  get_temp_ps, dummy, broadps

  x = findfile(broadps)
  xl = strlen(x)
  if xl(0) gt 0 then begin
    if !version.os eq 'vms' then spawn, 'delete '+broadps, list		$
    else spawn, '/bin/rm '+broadps, list
  endif

end
