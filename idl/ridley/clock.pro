
;------------------------------------------------------------------------------
;
;               Clock Dial Plots of Magnetometer Data
;
;------------------------------------------------------------------------------

PRO clock, single, t, times, moncheck, basetime, titlenam,    	      $
                  numofpts, placement, mullen, plotwhere, lat, mln,  	      $
		  numofvar, stat, datafile, titleofplot, step, publish,       $
		  nb, percent, opos, ani, fontsize, titlesize

  if n_elements(percent) eq 0 then percent=1.0
  cscale = percent*fontsize
  tscale = percent*titlesize

  opos = [1,1,0,0]

  relxlab = strarr(30)
  relxtl  = ' '
  btr       = 0.0
  etr       = 0.0
  curvar    = 0
  nl        = 0
  factor    = 100.0/3.0

;  Compute_Axis, times, placement, btr, etr, curvar, relxlab, relxtl,     $
;		     basetime, moncheck, nl, ticktime, nminor

  if (etr gt max(t)) and (max(t) ne 0.0) then etr=max(t)

  r0 = 90.0-lat
  range = 10.0*fix(max(r0)/10.0)+10.0

  a = findgen(16) * (!pi*2/16.)
  usersym, 0.3*cos(a), 0.3*sin(a), /fill

  theta0 = fltarr(numofvar/3+1)

  for i=0, numofvar/3 do begin

    if (numofpts(i*3) gt 0) and (numofpts(i*3+1) gt 0) then begin

      theta0(i) = mln(i)*!pi/12.0
      st_time = fltarr(6)
      b_time = basetime(i)
      c_r_to_a, st_time, b_time
      t0 = st_time(3)+(st_time(4)/60.0)+(st_time(5)/3600.0)
      theta0(i) = t0*!pi/12.0 - theta0(i) + !pi

    endif

  endfor

  nplots = numofpts(0)/((step+1)*nb) + 1

  if ani eq 1 then begin
    base = widget_base(title = 'Animation')
    animator = cw_animate(base,400,400,nplots)
    widget_control, base, /realize
    winnum = !d.window
  endif 
  frame = 0

  space = 0.04
  pos_space, nb, space, sizes

  if numofpts(0) gt 0 then begin

    zerotonb = 0

    for j=0, numofpts(0)-1, (step+1) do begin

;      if (nb gt 1) or (j eq 0) then begin
      if nb gt 0 then begin

        starttm = double(basetime(0)) + double(t(0,j))

        c_r_to_a, timearray, starttm

        times(0,0,0) = timearray(0)
        times(0,0,1) = timearray(1)
        times(0,0,2) = timearray(2)
        times(0,0,3) = timearray(3)
        times(0,0,4) = timearray(4)
        times(0,0,5) = timearray(5)

        if zerotonb eq 0 then datetl  = 				      $
	  strmid(moncheck,(times(0,0,1)-1)*3,3)+' '+			      $
	  strcompress(string(times(0,0,2)),/remove_all)+ ', '+ 		      $
	  strcompress(string(1900+times(0,0,0)),/remove_all)

        if times(0,0,3) gt 9 then 					      $
	  hour = strcompress(string(times(0,0,3)),/remove_all)		      $	
        else								      $
	  hour = strcompress('0'+string(times(0,0,3)),/remove_all)

        if times(0,0,4) gt 9 then 					      $
	  min = strcompress(string(times(0,0,4)),/remove_all)		      $	
        else								      $
	  min = strcompress('0'+string(times(0,0,4)),/remove_all)

        if times(0,0,5) gt 9 then 					      $
	  sec = strcompress(string(times(0,0,5)),/remove_all)		      $	
        else								      $
	  sec = strcompress('0'+string(times(0,0,5)),/remove_all)

        inserttl = hour+ ':'+ min + ':'+ sec + ' UT'

        if zerotonb eq 0 then plot,[0,1],xstyle=4,ystyle=4,/nodata

        get_position, nb, space, sizes, zerotonb,  $
		      pos

        plot, [-range-1,range+1], [-range-1, range+1],              	      $
              xstyle = 5,                      				      $
              ystyle = 5,      		       				      $
              title = inserttl,						      $
	      position = pos,						      $
              charsize = cscale,					      $
	      /nodata,							      $
	      /noerase

	if pos(0) lt opos(0) then opos(0) = pos(0)
	if pos(1) lt opos(1) then opos(1) = pos(1)
	if pos(2) gt opos(2) then opos(2) = pos(2)
	if pos(3) gt opos(3) then opos(3) = pos(3)

        a = findgen(361) * (!pi*2/360.)
        oplot, range*sin(a),range*cos(a)

        for i=0,range-10,10 do 						      $
          oplot, float(i)*sin(a),float(i)*cos(a), linestyle=1

        zerotonb = (zerotonb+1) mod nb

      endif

      for i=0, numofvar/3 do begin

        r1 = r0(i)

        if (single(i*3,j) ne -1e10) and (single(i*3+1,j) ne -1e10) then begin

          theta1 = theta0(i) + !pi*2.0*t(i,j)/86400.0
          theta2 = theta1 + !pi/2.0

          x1 = r1*sin(theta1)
          y1 = -r1*cos(theta1)

          x2dum = single(i*3,j) / factor
          y2dum = single(i*3+1,j) / factor

	  x2 = x1 + (x2dum * cos(theta2) - y2dum * sin(theta2))*mullen
	  y2 = y1 + (x2dum * sin(theta2) + y2dum * cos(theta2))*mullen

	  oplot, [x1], [y1], psym=8
	  oplot, [x1,x2], [y1,y2]

        endif

      endfor
                           
;      if (nb gt 1) or 							      $
      if ((j+step+1 gt numofpts(0)-1) or (zerotonb eq 0)) then begin

        if plotwhere eq 2 then begin

	  oplot, [-range],[-range], psym = 8
	  oplot, [-range,-range],[-range,-range+200.0/factor]
          xyouts, -range+range/200.0, -range,		      $
		  strcompress(string(fix(200.0/mullen+0.5))+' nT'),	      $
		  charsize=cscale

        endif

        if (publish eq 0) and 						      $
	   (plotwhere eq 1) and 					      $
	   ((j+step+1 gt numofpts(0)-1) or (zerotonb eq 0)) then begin

          chartime=''
          spawn,'show time',chartime

          side  = 'date printed : '+chartime+' from file(s) '+datafile(0)
          xyouts, 0, 0, side, orientation=90, charsize=0.5*cscale, /device

          xyouts, 0.5, 1.01, titleofplot, charsize=2.0*tscale, /norm, 	      $
	    alignment=0.5

        endif

        if (plotwhere eq 1) and 					      $
	   ((j+step+1 gt numofpts(0)-1) or (zerotonb eq 0)) then begin

          xd = pos(2) - pos(0)
          yd = pos(3) - pos(1)
          pos = [0.0, 0.0, xd, yd]

          plot, [-range-1,range+1], [-range-1,range+1],			      $
              xstyle    = 5,						      $
              ystyle    = 5,						      $
              pos=pos,							      $
	      /noerase,							      $
	      /nodata

	  oplot, [-range],[-range], psym = 8
	  oplot, [-range,-range+200.0/factor],[-range,-range]
          xyouts, 0.0, -0.05, 						      $
		  strcompress(string(fix(200.0/mullen+0.5))+' nT'),	      $
		  /norm, alignment = 0, charsize = 1.5*cscale

          xyouts, 0.5, -0.05, datetl, charsize=1.5*cscale, /norm, 	      $
		  alignment=0.5
          xyouts, 1.0, -0.05, 'Maximum Invariant Latitude : '+		      $
		  strcompress(string(round(90-range)),/remove_all)+'!Eo',     $
	          charsize=1.5*cscale, /norm, alignment=1.0

        endif

;	if (nb gt 1) and (zerotonb eq 0) and (plotwhere eq 2) then begin
	if (zerotonb eq 0) and (plotwhere eq 2) then begin

          xyouts, 0.5, 0.95, datetl, charsize=cscale, /norm, alignment=0.5
	  if (j+step+1 lt numofpts(0)-1) and (ani ne 1) then prompt_for_next

	endif

	if (ani eq 1) and 				$
	   ((zerotonb eq 0) or (j+step+1 gt numofpts(0)-1)) then begin
	  cw_animate_load, animator, window=winnum, 	$
		frame = frame
	  frame = frame + 1
	endif 

      endif

    endfor

    if ani eq 1 then begin
      cw_animate_run, animator
      xmanager, 'cw_animate thing', base, event_handler = 'animate_event'
    endif

  endif

  RETURN

END
