
; -------------------------------------------------------------------
;  MAIN WIDGET
; -------------------------------------------------------------------


pro d_vars_bat, datafile

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

    for i=0,nfile-1 do begin

      j = 0
      while (strmid(datafile(i),j,1) ne ']') and 			$
	    (j lt strlen(datafile(i))) do j=j+1
      if j lt strlen(datafile(i)) then j=j+1 else j=0
      fn = strmid(filename(i),j,strlen(datafile(i))-j)
      print, ' '
      print, 'File ',tostr(i+1),'     : ',fn
      print, 'start time : ',starttime(i)
      print, 'end time   : ',endtime(i)
      print, ' '
      for j = 0, varcount(i),2 do begin
	str1 = tostr(i+1)+', '+tostr(j+2)+'.  '+titles(i,j)
	for k=1,40 do str1 = str1 + ' '
	str1 = strmid(str1,0,40)
	if j+1 le varcount(i) then begin
	  str2 = tostr(i+1)+', '+tostr(j+3)+'.  '+titles(i,j+1)
	  for k=1,40 do str2 = str2 + ' '
	  str2 = strmid(str2,0,40)
	endif else str2=''
	print, str1+str2
      endfor
    endfor

;
;   type_8
;

    if plottype eq 8 then begin

      print, 'Need 5 Variables:'
      print, 'Range, Azimuth, Velocity, Elevation, Event'

    endif

;
;   type_13
;

    if plottype eq 13 then begin

      print, 'Need 6 Variables:'
      print, 'Range, Azimuth, Radial Vel, Angular Vel, Elevation, Event'

    endif

;
;   type_9
;

    if plottype eq 9 then begin

      print, 'Need 4 Variables:'
      print, 'Range, Azimuth, Velocity, Elevation'

    endif

  RETURN

END


