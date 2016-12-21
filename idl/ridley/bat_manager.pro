pro bat_manager

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

  done = 0

  while done eq 0 do begin

    stoppro = 0
    print, ' '
    print, 'Enter command :'
    command = ''
    read, command
    bat_event, command

    if stoppro eq 1 then done = 1
    if getfiles eq 1 then done = 1 
    if getplottype eq 1 then done = 1 

  endwhile

  return

end
