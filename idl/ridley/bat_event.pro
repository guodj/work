pro bat_event, ev

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
  common rep_data, change
  common data, single, t, dt, basetime, numofvar, numofpts
  common text_id, intext, inps, instring, inmaxx, inmaxy, inminx, inminy 
  common control, getplottype, getfiles, stoppro

  bs = byte(ev)
  loc = where((bs ge 97) and (bs lt 122), count)
  if count gt 0 then bs(loc) = bs(loc)-byte(32)
  ev = string(bs)

  case (ev) of

    'PLOT' : begin
      ev = {value : ev}
      widplote, ev, /screen
    end

    'SFLD' : begin
      ev = {value : ev}
      widplote, ev, /screen
    end

    'PSNAME' : begin
      print, 'Old ps file name : ',psfile
      print, 'Enter new ps file name :'
      read, psfile
      print, 'PS file now set to ',psfile
    end

    'PSL' : begin
      ev = {value : 'Make PS File (landscape)'}
      widplote, ev, /screen
    end

    'EPSL' : begin
      ev = {value : 'Make EPS File (landscape)'}
      widplote, ev, /screen
    end

    'PSP' : begin
      ev = {value : 'Make PS File (portrait)'}
      widplote, ev, /screen
    end

    'EPSP' : begin
      ev = {value : 'Make EPS File (portrait)'}
      widplote, ev, /screen
    end

    'OL' : begin
      ev = {value : 'Print plot to oosik (landscape)'}
      widplote, ev, /screen
    end

    'OP' : begin
      ev = {value : 'Print plot to oosik (Portrait)'}
      widplote, ev, /screen
    end

    'DONE' : stoppro = 1
    'STEP' : stoppro = 1
    'QUIT' : stoppro = 1
    'EXIT' : stoppro = 1
    'Q' : stoppro = 1
    'E' : stoppro = 1

    'FILE' : begin
      stoppro = 0
      getfiles = 1
    end 

    'PTYPE' : begin
      stoppro = 0
      getplottype = 1
    end

    'SELECT' : begin
      ev = {value : 'SELECT'}
      widplote, ev, /screen
    end 

    'MEAN' : begin
      ev = {value : 'MEAN'}
      widplote, ev, /screen
    end 

    'PUBLISH' : begin
      ev = {value : 'PUBLISH'}
      widplote, ev, /screen
    end 

    'LOGX' : begin
      ev = {value : 'LOGX'}
      widplote, ev, /screen
    end 

    'LOGY' : begin
      ev = {value : 'LOGY'}
      widplote, ev, /screen
    end 

    'CLR' : begin
      ev = {value : 'CLR'}
      widplote, ev, /screen
    end 

    'STATION' : begin
      ev = {value : 'STATION'}
      widplote, ev, /screen
    end 

    'INVAR' : begin
      ev = {value : 'INVAR'}
      widplote, ev, /screen
    end 

    'HELP' : d_menu_bat
    '?' : d_menu_bat

    'DISPLAY' : d_vars_bat, filename
    'DISP' : d_vars_bat, filename

; Things to do :

    'VAR' : begin
      get_vars_bat,titles, varcount, numofvar
      change.data = 1
      change.place = 1
    end

    'TIME' : begin
      get_time_bat, time
      change.data = 1
    end

    'RANGE' : begin
      print, ' '
      print, 'APSR  -  All plots same range'
      print, 'SPSR  -  Some plots same range'
      print, 'APSS  -  All plots same scale'
      print, 'SPSS  -  Some plots same scale'
      print, 'NONE  -  None of the above'
      print, 'Enter choice :'
      que = ''
      read, que

      bs = byte(que)
      loc = where((bs ge 97) and (bs lt 122), count)
      if count gt 0 then bs(loc) = bs(loc)-byte(32)
      que = string(bs)

      case (que) of

	'APSR' : begin
	  maxminque = 1
	  print, 'All plots will now be printed will same range'
	end

	'APSS' : begin
	  maxminque = 3
	  print, 'All plots will now be printed will same scale'
	end

	'SPSR' : begin
	  maxminque = 2
	  change.range = 1
	  print, 'Some plots will now be printed will same range'
	end

	'SPSS' : begin
	  maxminque = 4
	  change.range = 1
	  print, 'Some plots will now be printed will same range'
	end

	'NONE' : begin
	  maxminque = 5
	  print, 'All plots will have their individual ranges'
	end

	else : print, 'Nothing done'

      endcase

    end

    'XYRAN' : begin

      if strlen(inminx) gt 0 then print, 'Old min for X is ',inminx
      print, 'Enter new min for X (return for none) : '
      read, inminx
      if strlen(inminy) gt 0 then print, 'Old min for Y is ',inminy
      print, 'Enter new min for Y (return for none) : '
      read, inminy
      if strlen(inmaxx) gt 0 then print, 'Old max for X is ',inmaxx
      print, 'Enter new max for X (return for none) : '
      read, inmaxx
      if strlen(inmaxy) gt 0 then print, 'Old max for Y is ',inmaxy
      print, 'Enter new max for Y (return for none) : '
      read, inmaxy

    end

    'SYMBOL' : new_symbol, /screen

    'MAXR' : begin

      print, 'Old maximum range = ',maxran
      print, 'Enter new maximum range : '
      read, maxran
      print, 'Maximum range now set to ',maxran

    end

    'PPP' : begin

      print, 'Old number of plots per page = ',nb
      print, 'Enter number of plots per page : '
      read, nb
      print, 'Plots per page now set to ',nb

    end

    'DG' : begin

      if unconpts gt 0 then 					$
	print, 'Old number of seconds for a data gap = ',unconpts
      print, 'Enter number of seconds for a data gap : '
      read, unconpts
      print, 'Minimum data gap now set to ',tostr(unconpts),' seconds.'

    end

    'SM' : begin

      if offmul ne 1.0 then 					$
	print, 'Old offset multiplication factor = ',offmul
      print, 'Enter new offset multiplication factor : '
      read, offmul
      print, 'Offset multiplication factor set to ',offmul

    end

    'PTS' : begin

      if step ne 1 then 					$
	print, 'Old number of points to skip in between plots = ',step
      print, 'Enter new number of points to skip in between plots : '
      read, step
      print, 'I will now skip ',tostr(step),' time steps when plotting.'

    end

    'VEL' : begin

      if velocity ne 0.0 then 					$
	print, 'Old phase velocity = ',velocity
      print, 'Enter new phase velocity : '
      read, velocity
      print, 'Phase velocity now changed to ',velocity

    end

    'BEAM' : begin

      if beam ne 361.0 then 					$
	print, 'Old beam angle for imitation dwell = ',beam
      print, 'Enter new beam angle for imitation dwell : '
      read, beam
      print, 'Beam angle now set for ',beam

    end

    'COL' : begin

      if columns ne 1 then 					$
	print, 'Old number of columns = ',columns
      print, 'Enter number of columns : '
      read, columns
      print, 'Number of columns now set to ',columns
    end

    'TITLE' : begin

      if strlen(titleofplot) eq 0 then 				$
	print, 'Old title = ',titleofplot
      print, 'Enter new title : '
      read, titleofplot
      print, 'Title of plot is now "',titleofplot,'"'

    end

    'ROW' : begin

      if inarow ne 1 then 					$
	print, 'Old number of plots in a row = ',inarow
      print, 'Enter number of plots in a row : '
      read, inarow
      print, 'Number of plots in a row now set to ',inarow

    end

    'PERCENT' : begin

      if percent ne 1.0 then 					$
	print, 'Old reduction factor on PS file = ',percent
      print, 'Enter new reduction factor on PS file : '
      read, percent
      print, 'PS file will now be shrunk to ',percent

    end

    'PERFONT' : begin

      if change.fontsize ne 1.0 then 					$
	print, 'Old factor on font size = ',change.fontsize
      print, 'Enter new factor on font size : '
      read, change.fontsize
      print, 'New factor for font size ',change.fontsize

    end

    'PERTITLE' : begin

      if change.titlesize ne 1.0 then 					$
	print, 'Old factor on title size = ',change.titlesize
      print, 'Enter new factor on title size : '
      read, change.titlesize
      print, 'New factor for title size : ',change.titlesize

    end

    'PSN' : begin

      if strlen(psfile) eq 0 then 				$
	print, 'Old PS file name = ',psfile
      print, 'Enter new PS file name : '
      read, psfile
      print, 'The new postscript file name is ',psfile

    end



    else: begin

      print, 'Not an Option. Please type ? or HELP for a options list'

    endelse

  endcase

  return

end

