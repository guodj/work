pro plotvar, grafvar, titlenam, unitnam, times, screen, ani 

  common timeblock, sid, time, bftime, eftime
  common varblock, nfile, lid, titles, units, lfn
  COMMON type, plottype, types
  common fileblock, starttime, endtime, varcount, ftime 
  common extras, unconpts, columns, symtype, maxminque, rangepl, numofrang,   $
	         titleofplot, plotwhere, psfile, placement, numofgra,	      $
		 psfont, offdiv, offmul, maxran, checkinvar, filename, nb,    $
		 publish, checkstat, statfile, xyzpos, stat, lat, lon, alpha, $
		 mln, remmean, velocity, beam, step, crb, clr, draw,	      $
		 final_x, final_y, arrow_x, arrow_y, arrow, str_x, str_y,     $
		 str_s, str_o,str_l, str, checkselect, percent, inarow, ext
  common rep_data, change
  common data, single, t, dt, basetime, numofvar, numofpts
  common months, mon
  common flag_structure, flags
  common comdis, displays
  common descriptions, descript

  if !d.name eq 'PS' then plotwhere = 1 else plotwhere = 2

  if checkinvar eq 2 then same = 0 else same = 1

;  if ani eq 1 then begin
;    base = widget_base(title='Animation')
;    animator = cw_animate(base, 600,600, inarow)
;    widget_control, base, /realize
;    winnum = !d.window
;  endif

  if numofvar ge 0 then begin

    if (maxminque ne 2) and (maxminque ne 4) then begin
      rangepl = intarr(numofvar+1)
      numofrang = numofvar
    endif

    if change.range eq 1 then begin
      change.range = 0
      rs = change.rangescale
      get_range_int, maxminque, rangepl, numofvar, numofrang, 		$
	titlenam, rs, screen
      change.rangescale = rs
    endif

    seconds = [' ',strcompress(string(indgen(60)),/remove_all)]
    seconds(1:10) = '0'+seconds(1:10)

    rows = fix(float(numofgra+1)/float(columns) + 1.0 - 1.0/float(columns))
    !P.Multi = [0,columns,rows]

    if inarow eq 1 then begin
      times(0,0,0:5) = time(0:5)
      times(0,1,0:5) = time(6:11)
      times(0,*,0) = times(0,*,0) - 1900
      if checkinvar eq 2 then last = numofvar-3 else last = numofvar
      if last ge 1 then for i=1,last do times(i,*,*) = times(0,*,*)
    endif

    for np = 0, inarow-1 do begin

      if (change.data eq 1) or (inarow gt 1) then begin

	if inarow gt 1 then begin
	  if np eq 0 then times(0,1,*) = times(0,0,*)
          timinc = intarr(6)
	  timinc(0) = time(12)
	  timinc(1) = 0
	  timinc(2:5) = time(13:16)
	  if checkinvar eq 2 then					$
	    dumtime = times(0:numofvar-3,*,*)				$
	  else dumtime = times
	  increment, dumtime, timinc
	  if checkinvar eq 2 then					$
	    times(0:numofvar-3,*,*) = dumtime				$
	  else times = dumtime
	endif

        numofpts = intarr(numofvar+1)
        basetime = dblarr(numofvar+1)
	if change.flat(0) eq 1 then begin
	  if mklower(!version.os) eq 'vms' then 			$
            get_data, single, t, dt, basetime, numofvar, numofpts, 	$
		plottype, times, filename, grafvar, same		$
	  else				 				$
            get_ffunix_data, single, t, dt, basetime, numofvar, 	$
		numofpts, plottype, times, filename, grafvar, same
	endif else							$
          get_cdf_data, single, t, dt, basetime, numofvar, numofpts, 	$
		  plottype, times, filename, grafvar, same, descript

        change.data = 0

      endif

      case (plottype) of

        -1 : begin

	  make_ascii, single, t, basetime, numofvar, titlenam, 		$
		unitnam, placement, numofpts, filename, psfile

        end

        1 : begin

	  line_plt, single, t, dt, times, mon, seconds, same, 	$
		basetime, rows, columns, numofvar, titlenam, 		$
		unitnam, placement, maxminque, numofpts, unconpts, 	$
		rangepl, numofrang, filename, titleofplot, 		$
		plotwhere, symtype, publish, percent, pos,		$
		change.fontsize, change.titlesize, change.rangescale        
          ext.opos = pos

        end

        2 : begin

	  multi_v, single, t, dt, times, mon, seconds, same, 	$
		basetime, rows, columns, numofgra, titlenam, 		$
		unitnam, placement, maxminque, numofpts, unconpts, 	$
		filename, titleofplot, plotwhere, numofrang, 		$
		rangepl, symtype, percent, pos,				$
		change.fontsize, change.titlesize        
          ext.opos = pos

        end

        3 : begin

	  offset = 0.0

	  stack_pl, single, t, dt, times, mon, seconds, same, 	$
		basetime, rows, columns, numofvar, titlenam, 		$
		unitnam, placement, maxminque, numofpts, unconpts, 	$
		offset, filename, titleofplot, plotwhere, offdiv, 	$
		offmul, symtype, publish, percent, pos,			$
		change.fontsize, change.titlesize, remmean        
          ext.opos = pos

        end

        4 : begin

	  if symtype ne 0 then dot = 2 else dot = 1
	  xy_plt, single, t, dt, times, numofgra, titlenam, 	$
		unitnam, placement, maxminque, numofpts, unconpts, 	$
		dot, seconds, filename, titleofplot, plotwhere, 	$
		symtype, flags.maxx, flags.maxy, flags.minx, 		$
		flags.miny, flags.logx, flags.logy, publish, 		$
		percent, pos, change.fontsize, change.titlesize        
          ext.opos = pos

        end

        8 : begin

	  mullen = offmul*offdiv
	  survey = 0.0

	  az_scan, single, t, times, mon, seconds, basetime, 	$
		titlenam, unitnam, numofpts, placement, mullen, 	$
		maxran, plotwhere, survey, filename, titleofplot, 	$
		checkinvar, nb, publish, checkstat, crb, clr, 		$
		displays.dshape, displays.dswitch,			$
		displays.dstring, displays.dchars, percent, pos,	$
		change.fontsize, change.titlesize, displays.ddate,	$
		change.minvel
          ext.opos = pos

        end

        9 : begin

	  mullen = offmul*offdiv
          publish = publish + 1
	  dwell, single, t, times, mon, seconds, basetime, 	$
		titlenam, unitnam, numofpts, placement, mullen, 	$
		maxran, plotwhere, filename, titleofplot, publish, 	$
		crb, clr, percent, pos,	change.fontsize, change.titlesize
          publish = publish - 1
          ext.opos = pos

        end

        10 : begin

	  mullen = offmul*offdiv

	  eq_con, single, t, times, mon, seconds, basetime, 	$
		titlenam, unitnam, numofpts, placement, mullen, 	$
		plotwhere, lat, lon, velocity, numofvar, stat, 		$
		filename, titleofplot, alpha, beam, percent, pos,	$
		change.fontsize, change.titlesize        
          ext.opos = pos

        end

        13 : begin

	  mullen = offmul*offdiv
	  survey = 0.0
	  givenaz = 361.0
	  rotation = 0.0

	  radar_ef, single, t, times, mon, seconds, basetime, 	$
		titlenam, unitnam, numofpts, placement, mullen, 	$
		maxran, plotwhere, survey, filename, titleofplot, 	$
		givenaz, rotation, clr, nb, publish, displays.dshape, 	$
		displays.dswitch, displays.dstring, displays.dchars, 	$
		percent, pos, change.fontsize, change.titlesize,	$
		displays.ddate
          ext.opos = pos
	  ;checkinvar, nb, publish,checkstat, crb, clr

        end

        20 : begin

	  offset = 0.0

	  if numofvar ge 11*3 then firstnum = (numofvar+1)/6  		$
	  else firstnum=(numofvar+1)/3

	  mag_qick, single, t, dt, times, mon, seconds, same, 	$
		basetime, numofvar, titlenam, unitnam, placement,	$
		maxminque, numofpts, unconpts, offset, filename,	$
		titleofplot, stat, firstnum, plotwhere, offdiv, 	$
		offmul, symtype, remmean, percent, pos,			$
		change.fontsize, change.titlesize
          ext.opos = pos

        end

        21 : begin

	  mullen = offmul*offdiv

	  clock, single, t, times, mon, basetime,titlenam, 	$
		numofpts, placement, mullen, plotwhere, lat, mln, 	$
		numofvar, stat, filename, titleofplot, step, 		$
		publish, nb, percent, pos, ani,				$
		change.fontsize, change.titlesize
          ext.opos = pos

        end

        22 : begin

	  if symtype ne 0 then dot = 2 else dot = 1
	  xy_plt, single, t, dt, times, numofgra, titlenam, 	$
		unitnam, placement, maxminque, numofpts, unconpts, 	$
		dot, seconds, filename, titleofplot, plotwhere, 	$
		symtype, flags.maxx, flags.maxy, flags.minx, 		$
		flags.miny, flags.logx, flags.logy, publish, 		$
		percent, pos, change.fontsize, change.titlesize, /t22	
          ext.opos = pos

        end

        else : print, 'Plot type not available yet.'

      endcase

      if (checkinvar eq 2) and (np eq inarow-1) then begin

        numofvar = numofvar - 3
        grafvar = grafvar(0:numofvar,*)
        times = times(0:numofvar,*,*)

      endif

      if (np ne inarow-1) and (!d.name eq 'X') and (ani eq 0) then 	$
	prompt_for_next 

;      if ani eq 1 then cw_animate_load, animator, window=winnum, frame = np

    endfor

    if !d.name eq 'Z' then begin
      bytemap = tvrd()
      tvlct, rr, gg, bb, /get
      write_gif, psfile, bytemap, rr, gg, bb
      device, /close
      set_plot,'X'
    endif

;    if ani eq 1 then begin
;      cw_animate_run, animator
;      xmanager, 'cw_animate thing', base, event_handler='animate_event'
;    endif

  endif

  if inarow gt 1 then change.data = 1

  return

end

