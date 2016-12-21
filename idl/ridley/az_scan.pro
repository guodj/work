;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;               Polar Vector Survey Plots
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO az_scan, single, t, times, moncheck, seconds, basetime, titlenam,     $
		 unitnam, numofpts, placement, mullen, maxran, plotwhere,     $
		 survey, datafile, titleofplot, checkinvar, nb, publish,      $
		 checkstat, crb, clr, dshape, dswitch, dstring, dchars,	      $
		 percent, opos, fontsize, titlesize, ddate, minvel

  if n_elements(dshape) eq 0 then dshape = intarr(200)
  if n_elements(dswitch) eq 0 then dswitch = intarr(200)
  if n_elements(dstring) eq 0 then dstring = strarr(200,3)
  if n_elements(dchars) eq 0 then dchars= intarr(200,3)+1.0
  if n_elements(percent) eq 0 then percent = 1.0
  if n_elements(minvel) eq 0 then minvel = 20.0
  cscale = percent*fontsize
  tscale = percent*titlesize

  opos = [1.0,1.0,0.0,0.0]

  red     = intarr(256)
  blue    = intarr(256)
  green   = intarr(256)  
  relxlab = strarr(30)

  checkplace = where(placement ne 0, count)

  if count eq 0 then begin

    altitude  = 0
    azimuth   = 1
    velocity  = 2
    elevation = 3
    event     = 4

  endif else begin

    altitude  = placement(0,0)
    azimuth   = placement(0,1)
    velocity  = placement(0,2)
    elevation = placement(0,3)
    event     = placement(0,4)

  endelse

  curvar    = 0

  sizeoftype = 1.5

  minalt = 0.0
  maxalt = maxran

  ytlr =  ' (' + strcompress(unitnam(altitude), /remove_all) + ')'

  ytl  = strcompress(titlenam(altitude), /remove_all) + ' ' + ytlr

  if clr eq 1 then begin

    red(0:2)      = [0,0,255]
    green(0:2)    = [0,0,0]
    blue(0:2)     = [0,255,0]

  endif else begin

    if plotwhere eq 2 then begin

      red(0:2)      = [0,200,100]
      green(0:2)    = [0,200,100]
      blue(0:2)     = [0,200,100]

    endif else begin

      red(0:2)      = [0,50,120]
      green(0:2)    = [0,50,120]
      blue(0:2)     = [0,50,120]

    endelse

  endelse

  red(3:255)    = 255
  green(3:255)  = 255
  blue(3:255)   = 255
  tvlct, red, green, blue

  space = 0.05
  pos_space, nb, space, sizes

  numofeve = 0

  if numofpts(event) gt 1 then begin

    eventdum = single(event, 1:numofpts(event)-1) - 			      $
	       single(event, 0:numofpts(event)-2)

    eventdum2 = where(eventdum ne 0, numofeve)

  endif

  startpt = 0
  dummy   = 0
  starttm = double(0.0)
  onetofour = 0

  aaa = findgen(361) * (!pi/180.)

  for cureve = 0, numofeve do begin

    if crb eq 1 then begin
      crb_xlc = fltarr(1)
      crb_ylc = fltarr(1)
    endif

    found = 0

    starttm = double(basetime(event)) + double(t(event,startpt))

    c_r_to_a, timearray, starttm

    times(0,0,0) = timearray(0)
    times(0,0,1) = timearray(1)
    times(0,0,2) = timearray(2)
    times(0,0,3) = timearray(3)
    times(0,0,4) = timearray(4)
    times(0,0,5) = timearray(5)

    if (onetofour eq 0) then datetl  =	 				$
	strmid(moncheck,(times(0,0,1)-1)*3,3)+' '+			$
	strcompress(string(times(0,0,2)),/remove_all)+ ', '+ 		$
	strcompress(string(1900+times(0,0,0)),/remove_all)

    inserttl = strcompress(seconds(1+times(0,0,3)),/remove_all)+ ':'+ 	      $
	       strcompress(seconds(1+times(0,0,4)),/remove_all)+ ':'+ 	      $
               strcompress(seconds(1+times(0,0,5)),/remove_all)+ ' UT'

    elev = single(elevation,startpt)

    insertel = 'Elevation : '+						      $
	       strcompress(string(round(single(elevation,startpt))),	      $
	       /remove_all)+'!Eo'

    mulfac = 0.95
    ystart = 0.04

    get_position, nb, space, sizes, onetofour, pos

    a = findgen(16) * (!pi*2./16.)

    usersym, 0.25*cos(a)*cscale, 0.25*sin(a)*cscale, /fill

    if onetofour eq 0 then plot,[0,1],xstyle=5,ystyle=5,/nodata

    plot, [-maxalt-1,maxalt+1], [-maxalt-1,maxalt+1],			      $
          charsize  = cscale,						      $
          ytitle    = ytl,						      $
          xstyle    = 5,						      $
          ystyle    = 5,						      $
          pos=pos,							      $
	  title = inserttl,						      $
	  /noerase,							      $
	  /nodata

    if pos(0) lt opos(0) then opos(0) = pos(0)
    if pos(1) lt opos(1) then opos(1) = pos(1)
    if pos(2) gt opos(2) then opos(2) = pos(2)
    if pos(3) gt opos(3) then opos(3) = pos(3)

    case dshape(cureve mod nb) of 

      0 : begin

            plots, maxalt*sin(aaa),maxalt*cos(aaa)

	  end

      1 : begin

	    plots, maxalt*[-1.0,-1.0,1.0,1.0,-1.0],			$
		   maxalt*[-1.0,1.0,1.0,-1.0,-1.0]

	  end

      2 : begin

            plots, maxalt*sin(aaa),maxalt*cos(aaa)

	    oplot, [-maxalt,maxalt],[0,0]
	    oplot, [0,0],[-maxalt,maxalt]

	  end

    endcase

    onetofour = (onetofour + 1) mod nb

    if numofpts(velocity)-startpt gt 0 then begin
  
      facvel = mullen * (maxalt-minalt) / 5000.0
      factim = 1.0 / 2000.0

      curpt  = startpt
      eveold = single(event,startpt)

      while (curpt lt numofpts(velocity) - 1) and 			      $
	    (eveold eq single(event,curpt)) do begin
                                                 
        if (single(velocity,curpt) ne -1e10) and 			$
	   (single(altitude,curpt) le maxalt) and			$
	   (abs(single(velocity,curpt)) lt 2000.0) and			$
	   (abs(single(velocity,curpt)) gt minvel) then begin

	  if single(velocity,curpt) ge 0.0 then colval=2 else colval=1

	  r   = single(altitude,curpt)
          t1  = single(azimuth,curpt) * 2.0 * !pi / 360.0
          v   = single(velocity,curpt)*facvel
          f   = 15.0 * 2.0 * !pi / 360.0
      
          s  = v * sin(f)
          p  = r + v * cos(f)
          sp = sqrt(s*s + p*p)
          if sp ne 0.0 then t2 = asin(s/sp) else t2=0.0
	  if p lt 0.0 then t2 = !pi - t2
          if cos(t2) ne 0.0 then e  = p / cos(t2) else e=p

          x1 = r * sin(t1)
          x2 = e * sin(t1 + t2)

          y1 = r * cos(t1)
          y2 = e * cos(t1 + t2)

	  plots, [x1], [y1], psym = 8, color = colval

          plots, [x1,x2],[y1,y2], color = colval

          if (found eq 0) and (crb eq 1) and (curpt gt 0) then begin

	    sign = single(velocity, curpt)*single(velocity, curpt-1)

	    if (sign lt 0) and 						      $
	      (t(altitude, curpt) eq t(altitude, curpt-1)) then begin

	      found = 1

	      ypos = single(altitude,curpt)*				      $
		      sin(!pi/2.0 - !pi*single(azimuth,curpt)/180.0)
	      xpos = single(altitude,curpt)*				      $
		      cos(!pi/2.0 - !pi*single(azimuth,curpt)/180.0)
	      if crb_xlc(0) eq 0 then begin
	        crb_ylc(0) = ypos
	        crb_xlc(0) = xpos
	      endif else begin 
	        crb_ylc = [crb_ylc,ypos]
	        crb_xlc = [crb_xlc,xpos]
              endelse

            endif

	  endif

          if (curpt gt 0) and (found eq 1) and (crb eq 1) then		      $ 
            if (t(altitude, curpt) ne t(altitude, curpt-1)) then 	      $
	      found = 0

        endif

        eveold = single(event,curpt)
        curpt  = curpt + 1

      endwhile

      if plotwhere eq 2 then begin

        x1 = -maxalt
        x2 = x1 + 1000.0*facvel/5.0
        y1 = x1
        y2 = y1 + 1000.0*facvel

        plots, [x1], [y1], psym = 8

        plots, [x1,x2], [y1,y2]

        xyouts, x2 + 1000.0*facvel/4.0, y1, '1.0 Km/s', charsize=0.9

      endif

      if checkinvar eq 2 then begin

        if titlenam(0) eq 'Range' then 					      $
	  maxrange = maxalt*cos(min(					      $
	    single(elevation,0:numofpts(elevation)-1))*!pi/180.0) 	      $
        else maxrange = maxalt
        data_array = single(5:7,*)
        basic_contour, data_array, maxrange, pos

      endif

      if checkstat eq 2 then begin

        disfac = 1.0

        elev = 30.0

        Put_Stations, disfac, elev, cscale

      endif

      if (crb eq 1) and (n_elements(crb_xlc) gt 0) then begin
	oplot, crb_xlc, crb_ylc
      endif

    endif

    startpt = curpt

    if (numofeve eq cureve) or (onetofour eq 0) then begin

      if plotwhere eq 2 then begin

	if not ddate then 						$
        xyouts, 0.5, 0.95, datetl, charsize=1.2*cscale, /norm, alignment=0.5
        xyouts, 1.0, 0.0, 'Maximum Range : '+				      $
		strcompress(string(round(maxalt)),/remove_all)+' Km',         $
	        charsize=1.2*cscale, /norm, alignment=1.0
        xyouts, 0.0, 0.0, insertel, charsize=1.2*cscale, /norm

      endif else begin

        if publish eq 0 then begin

          chartime=''
          spawn,'show time',chartime

          side  = 'date printed : '+chartime+' from file(s) '+datafile(0)+    $
	          '   '+datafile(1)
          xyouts, -0.01, 0.05, side, orientation=90, charsize=0.5*cscale,    $
		/norm

	endif

	ymin = opos(1)-space
	xmin = opos(0)-space/2.0+0.005
	xmax = opos(2)+space/2.0-0.005
	xmid = (xmin+xmax)/2.0
	if not ddate then 						$
        xyouts, xmid, ymin+0.01, datetl, charsize=1.2*cscale, /norm, 	$
		alignment=0.5
        xyouts, xmax, ymin+0.01, 'Maximum Range : '+			$
		strcompress(string(round(maxalt)),/remove_all)+' Km',   $
	        charsize=1.2*cscale, /norm, alignment=1.0
        xyouts, xmin, ymin+0.01, insertel, charsize=1.2*cscale, /norm

        xd = pos(2) - pos(0)
        yd = pos(3) - pos(1)
        pos = [xmin, ymin+space/2.0, xmin+xd, ymin+yd+space/2.0]

        plot, [-maxalt-1,maxalt+1], [-maxalt-1,maxalt+1],		      $
              xstyle    = 5,						      $
              ystyle    = 5,						      $
              pos=pos,							      $
	      /noerase,							      $
	      /nodata

        x1 = -maxalt
        x2 = x1 + 1000.0*facvel/5.0
        y1 = x1 + 0.12*maxalt
        y2 = y1 + 1000.0*facvel

        plots, [x1], [y1], psym = 8

        plots, [x1,x2], [y1,y2]

        xyouts, x2 + 1000.0*facvel/6.0, y1, '1.0 Km/s', charsize=0.9*cscale

      endelse

      xyouts, (opos(0)+opos(2))/2.0, opos(3), titleofplot, 		$
	charsize=2.0*tscale, /norm, alignment=0.5

      if (plotwhere ne 1) and (numofeve ne cureve) then prompt_for_next

    endif

  endfor

  !p.multi = 0

  RETURN

END

