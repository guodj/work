;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;              Electric Field Plots
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO radar_ef, single, t, times, moncheck, seconds, basetime, titlenam,    $
		  unitnam, numofpts, placement, mullen, maxran, plotwhere,    $
		  survey, datafile, titleofplot, givenaz, rotation, clr,      $
		  nb, publish, dshape, dswitch, dstring, dchars, percent,     $
		  opos, fontsize, titlesize, ddate

  if n_elements(dshape) eq 0 then dshape = intarr(200)
  if n_elements(dswitch) eq 0 then dswitch = intarr(200)
  if n_elements(dstring) eq 0 then dstring = strarr(200,3)
  if n_elements(dchars) eq 0 then dchars = intarr(200,3)+1.0

  if n_elements(clr) eq 0 then clr = 0
  if n_elements(nb) eq 0 then nb = 1
  if n_elements(publish) eq 0 then publish = 0
  if n_elements(percent) eq 0 then percent = 1.0
  cscale = percent*fontsize
  tscale = percent*titlesize

  opos = [1.0,1.0,0.0,0.0]

  red     = intarr(256)
  blue    = intarr(256)
  green   = intarr(256)  
  relxlab = strarr(30)

  checkplace = where(placement ne 0, count)

  curvar    = 0

  if count eq 0 then begin

    altitude  = 0
    azimuth   = 1
    er        = 2
    et        = 3
    elevation = 4
    event     = 5

  endif else begin

    altitude  = placement(0,0)
    azimuth   = placement(0,1)
    er        = placement(0,2)
    et        = placement(0,3)
    elevation = placement(0,4)
    event     = placement(0,5)

  endelse

  if clr eq 1 then begin

    red(0:2)   = [0,0,255]
    green(0:2) = [0,0,0]
    blue(0:2)  = [0,255,0]

  endif else begin

    if plotwhere eq 2 then begin

      red(0:2)   = [0,200,100]
      green(0:2) = [0,200,100]
      blue(0:2)  = [0,200,100]

    endif else begin

      red(0:2)   = [0,50,180]
      green(0:2) = [0,50,180]
      blue(0:2)  = [0,50,180]

    endelse

  endelse

  red(3:255) = 255
  green(3:255) = 255
  blue(3:255) = 255
  tvlct, red, green, blue

  sizeoftype = cscale
  space = 0.04
  pos_space, nb, space, sizes

  minalt = 0.0
  maxalt = maxran

  ytlr =  ' (' + strcompress(unitnam(altitude), /remove_all) + ')'

  ytl  = strcompress(titlenam(altitude), /remove_all) + ' ' + ytlr

  numofeve = 0

  if numofpts(event) gt 1 then begin

    for i = 1, numofpts(event)-1 do 					      $
      if single(event,i) ne single(event,i-1) then numofeve = numofeve + 1

  endif

  startpt = 0
  dummy   = 0
  starttm = double(0.0)

  aaa = findgen(361) * (!pi*2/360.)

  a = findgen(16) * (!pi*2/16.)
  usersym, 0.15*cos(a), 0.15*sin(a), /fill

  for cureve = 0, numofeve do begin

    onetofour = cureve mod nb

    starttm = double(basetime(event)) + double(t(event,startpt))

    c_r_to_a, timearray, starttm

    times(0,0,0) = timearray(0)
    times(0,0,1) = timearray(1)
    times(0,0,2) = timearray(2)
    times(0,0,3) = timearray(3)
    times(0,0,4) = timearray(4)
    times(0,0,5) = timearray(5)

    datetl  = strmid(moncheck,(times(0,0,1)-1)*3,3)+' '+		      $
	       strcompress(seconds(1+times(0,0,2))+ ', 19'+    $
               strcompress(string(times(0,0,0)),/remove_all))

    inserttl = strcompress(seconds(1+times(0,0,3)) +			      $
                           ':'+seconds(1+times(0,0,4)),/remove_all)

    case dswitch(cureve mod nb) of

      0 : begin

	    inserttl = inserttl + ':'+		      $
		strcompress(seconds(1+times(0,0,5)),/remove_all)+ ' UT'

	  end

      1 : begin

	    inserttl = inserttl + ' UT'

	  end

      2 : begin

	    inserttl = dstring(cureve mod nb, 1)

	  end

    endcase

    elev = round(single(elevation,startpt))
    insertel = 'Elevation ='+strcompress(string(elev),/remove_all)+'!Eo'

    get_position, nb, space, sizes, onetofour, pos

    if onetofour eq 0 then begin

      plot, [-maxalt,-maxalt], [-maxalt,-maxalt],			$
            charsize  = sizeoftype, 					$
            xra       = [-maxalt, maxalt],				$
            yra       = [-maxalt, maxalt],				$
            ytitle    = ytl,						$
            xstyle    = 5,						$
            ystyle    = 5,						$
            pos = pos

    endif else begin

      plot, [-maxalt,-maxalt], [-maxalt,-maxalt],			$
            charsize  = sizeoftype, 					$
            xra       = [-maxalt, maxalt],				$
            yra       = [-maxalt, maxalt],				$
            ytitle    = ytl,						$
            xstyle    = 5,						$
            ystyle    = 5,						$
            pos       = pos,						$
	    /noerase

    endelse

    if pos(0) lt opos(0) then opos(0) = pos(0)
    if pos(1) lt opos(1) then opos(1) = pos(1)
    if pos(2) gt opos(2) then opos(2) = pos(2)
    if pos(3) gt opos(3) then opos(3) = pos(3)

    case dshape(cureve mod nb) of

      0 : begin

	    plots, maxalt*sin(aaa), maxalt*cos(aaa)

	  end

      1 : begin

	    plots, maxalt*[-1.0,-1.0,1.0,1.0,-1.0],			$
		   maxalt*[-1.0,1.0,1.0,-1.0,-1.0]

	  end

      2 : begin

	    plots, maxalt*sin(aaa), maxalt*cos(aaa)

	    plots, [-maxalt,maxalt],[0,0]
	    plots, [0,0],[-maxalt,maxalt]

	  end

    endcase

    if numofpts(er)-startpt gt 0 then begin
  
      facvel = (maxalt-minalt) / 5000.0

      curpt  = startpt
      eveold = single(event,startpt)

      dalt = single(altitude,startpt)/2.0
      dazi = single(azimuth,startpt)/2.0
      dazi = 0.0
      dalt = 0.0

      while (curpt le numofpts(er) - 1) and 				      $
	    (eveold eq single(event,curpt)) do begin
                                                 
        testazm = single(azimuth,curpt) - givenaz
	if single(azimuth,curpt) ge 180.0 then testazm = testazm - 180.0

;        if (givenaz gt 360.0) or 
;	((testazm lt 4.9) and (testazm gt 0.0)) then begin

	  if single(er,curpt) ge 0.0 then colval = 2 else colval = 1

          r   = single(altitude,curpt) - dalt
          t1  = (single(azimuth,curpt) - dazi) * 2.0 * !pi / 360.0
          efr = single(er,curpt)*facvel*mullen
          eft = single(et,curpt)*facvel*mullen

          x1 = r * sin(t1)
          x2 = x1 + eft*cos(t1) + efr*sin(t1)

          y1 = r * cos(t1)
          y2 = y1 - eft*sin(t1) + efr*cos(t1)

          if (single(er,curpt) ne -1e10) and				      $
	     (single(altitude,curpt) le maxalt) then begin

	    plots, [x1], [y1], psym = 8, color = colval
	    plots, [x1,x2],[y1,y2], color = colval

	  endif

;	endif

        eveold = single(event,curpt)
        curpt  = curpt + 1

      endwhile

      if plotwhere eq 2 then begin

        x1 = -maxalt + maxalt/20.0
        x2 = x1 + 1000.0*facvel/5.0
        y1 = x1
        y2 = y1 + 1000.0*facvel

        oplot, [x1], [y1], psym = 8

        oplot, [x1,x2], [y1,y2]

        xyouts, x2 + 0.05*facvel/4.0, y1, 				      $
	  strcompress(string(1000.0/mullen),/remove_all)+' '+unitnam(er),     $
	  charsize = 0.75*cscale

      endif

      xyouts, pos(0), pos(3)+0.015, 				$
	dstring(cureve mod nb,0), alignment=0.0, /norm, charsize=cscale
      xyouts, pos(2), pos(3)+0.015, 				$
	dstring(cureve mod nb,2), alignment=1.0, /norm, charsize=cscale
      xyouts, pos(0)+(pos(2)-pos(0))/2.0, pos(3)+0.015, 		$
	inserttl, alignment=0.5, /norm, charsize=cscale

;      disfac = 1.0

;      Put_Stations, disfac, elev

    endif

    startpt = curpt

    if (cureve eq nb-1) or (cureve eq numofeve) then begin

      if plotwhere eq 2 then begin
	if not ddate then						$
	  xyouts, 0.5, 0.95, datetl, charsize=1.5*cscale, /norm, 	      $
	  alignment = 0.5
        xyouts, 0.0, 0.00, insertel, charsize=1.5*cscale, /norm
        xyouts, 1.0, 0.00,'Maximum Range = '+				      $
          strcompress(string(round(maxalt)),/remove_all)+' Km', 	      $
          alignment=1.0, /norm, charsize = 1.5*cscale
      endif else begin
	if not ddate then						$
	xyouts, 0.5, -0.05, datetl, charsize=1.5*cscale, /norm, 	      $
	  alignment = 0.5
        xyouts, 0.0, -0.05, insertel, charsize=1.5*cscale, /norm
        xyouts, 1.0, -0.05,'Maximum Range = '+				      $
          strcompress(string(round(maxalt)),/remove_all)+' Km', 	      $
          alignment=1.0, /norm, charsize = 1.5*cscale

	xd = pos(2) - pos(0)
	yd = pos(3) - pos(1)
	pos = [0.0,0.0,xd,yd]

        plot, [-maxalt-1,maxalt+1], [-maxalt-1,maxalt+1],		$
            xstyle    = 5,						$
            ystyle    = 5,						$
            pos       = pos,						$
	    /nodata,							$
	    /noerase

        x1 = -maxalt + maxalt/20.0
        x2 = x1 + 1000.0*facvel/5.0
        y1 = x1
        y2 = y1 + 1000.0*facvel

        oplot, [x1], [y1], psym = 8

        oplot, [x1,x2], [y1,y2]

        xyouts, x2 + 0.05*facvel/4.0, y1, 				      $
	  strcompress(string(1000.0/mullen),/remove_all)+' '+unitnam(er),     $
	  charsize = 0.75*cscale

      endelse

    endif

    if (publish eq 0) and 						      $
       ((cureve eq nb-1) or (cureve eq numofeve)) then begin

      chartime=''
      spawn,'show time',chartime

      side  = 'date printed : '+chartime+' from file(s) '+datafile(0)+	      $
	'   '+datafile(1)
      side1 = datafile(2)+'   '+datafile(3)
      side2 = datafile(4)+'   '+datafile(5)
      xyouts, 0, 0, side, orientation=90, charsize=0.5*cscale, /device
      xyouts, 200, 2000, side1, orientation=90, charsize=0.5*cscale, /device
      xyouts, 400, 2000, side2, orientation=90, charsize=0.5*cscale, /device

    endif

    if (cureve eq nb-1) or (cureve eq numofeve) then begin

      xyouts, 0.5, 1.01, titleofplot, charsize=2.0*tscale, /norm, 	      $
	alignment=0.5

    endif

    if (plotwhere ne 1) and 						      $
       (cureve eq nb-1) and						      $
       (cureve ne numofeve) then prompt_for_next

  endfor

  !p.multi = 0

  RETURN

END

