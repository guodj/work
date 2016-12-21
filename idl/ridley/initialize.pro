pro initialize

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

  mon='JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'

  displays = { dshape : intarr(200), dswitch : intarr(200), 		$
	       dstring : strarr(200,3), dchars : fltarr(200,3)+1.0,	$
	       ddate : 0 }

  columns = 1
  maxminque=5
  unconpts = 10000000.0
  titleofplot=''
  plotwhere = 2
  symtype = 0
  titleofplot = ''
  psfile = ''

  if n_elements(change) gt 0 then begin
    ysize = change.y_win_size
  endif else begin
    ysize = 600
  endelse

  change = {data : 1, place : 1, invar : 1, range : 0, 		$
		fontsize : 1.0, titlesize : 1.0,		$
		y_win_size : ysize, flat : intarr(6),		$
		cdf : intarr(6), location : intarr(6,150,3),	$
		rangescale:strarr(100,2), minvel : 0}
  numofvar = -1
  offdiv = 1.0
  offmul = 1.0
  psfont = 28
  checkinvar = 1
  checkstat = 1
  nb = 4
  publish = 0
  remmean = 0
  beam = 361
  velocity = 0
  step = 0
  !P.font=0
  crb = 0
  clr = 1
  draw = 0
  arrow = 0
  str_x = -1
  checkselect = 0
  percent = 1.0
  inarow = 1
  flags = {todo, maxx : -1.0e-30, maxy : -1.0e-30, 			      $
		 minx : -1.0e-30, miny : -1.0e-30, 			      $
		 logx : 0,    logy : 0}
  ext = {opos : fltarr(4), line_x : fltarr(100), line_y : fltarr(100)}
  ext.line_x(0) = -1
  ext.line_y(0) = -1

  return

end
