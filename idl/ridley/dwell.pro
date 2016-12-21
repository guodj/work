;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;               Dwell Vector Plots
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO dwell, single, t, times, moncheck, seconds, basetime, titlenam,     $
                 unitnam, numofpts, placement, mullen, maxran, plotwhere,     $
		 datafile, titleofplot, final, crb, clr, percent, opos,	$
		 fontsize, titlesize

  if n_elements(percent) eq 0 then percent = 1.0
  cscale = percent*fontsize
  tscale = percent*titlesize

  red     = intarr(256)
  blue    = intarr(256)
  green   = intarr(256)  
  relxlab = strarr(30)
  relxtl  = ' '
  if crb eq 1 then begin
    crb_loc = fltarr(1)
    crb_tim = fltarr(1)
  endif

  btr       = 0.0
  etr       = 0.0
  curvar    = 0
  nl        = 0

  checkplace = where(placement ne 0, count)

  if count eq 0 then begin

    altitude  = 0
    azimuth   = 1
    velocity  = 2
    elevation = 3

  endif else begin

    altitude  = placement(0,0)
    azimuth   = placement(0,1)
    velocity  = placement(0,2)
    elevation = placement(0,3)

  endelse

  Compute_Axis, times, placement, btr, etr, curvar, relxlab, relxtl,          $
	basetime, moncheck, nl, ticktime, nminor

  sizeoftype = 1.25*cscale

  minalt = min(single(altitude,0:numofpts(altitude)-1))
  minalt = minalt - 0.25*(maxran - minalt)
  maxalt = maxran

  divfac = 0.25*(maxalt - minalt)

  ytlr =  ' (' + strcompress(unitnam(altitude), /remove_all) + ')'

  ytl  = strcompress(titlenam(altitude), /remove_all) + ' ' + ytlr

  if clr eq 1 then begin

    red(1:99)		= 0.0
    red(100:199)	= 100.0+ 155.0 * sqrt(findgen(100)/100.0)
    green(0:199)	= 0.0
    blue(1:99)		= 255.0 - 155.0*sqrt(findgen(99)/100.0)
    blue(100:199)	= 0.0

  endif else begin

    if plotwhere eq 2 then begin

      red(1:99)      = 200
      red(100:199)   = 100
      green(1:99)    = 200
      green(100:199) = 100
      blue(1:99)     = 200
      blue(100:199)  = 100

    endif else begin

      red(1:99)      = 50
      red(100:199)   = 180
      green(1:99)    = 50
      green(100:199) = 180
      blue(1:99)     = 50
      blue(100:199)  = 180

    endelse

  endelse

  red(0) = 0
  green(0) = 0
  blue(0) = 0
  red(200:255)   = 255
  green(200:255) = 255
  blue(200:255)  = 255

  tvlct, red, green, blue
  
  plot, [btr,etr], [minalt,minalt],                                $
        charsize  = sizeoftype,                                    $
        xticks    = nl,                                            $
        xtickname = relxlab,                                       $
        xra       = [btr,etr],                                     $
        yra       = [minalt, maxalt],				   $
        xtitle    = relxtl,                                        $
        ytitle    = ytl,                                           $
        xminor    = nminor,                                        $
        yminor    = 5,                                             $
        xstyle    = 1,                                             $
        ystyle    = 1,						   $
        xtickv=ticktime
                        
  ppos, opos

  if numofpts(velocity) gt 0 then begin

    facvel = mullen * (maxalt-minalt) / 10000.0
    timadd = (etr - btr)/20.0
    factim = 1.0 / 2000.0

    azdiv = 5.0
    oldaz = 0.0
    oldel = 0.0

    a = findgen(16) * (!pi*2/16.)

    usersym, 0.5*cos(a), 0.5*sin(a), /fill

    found = 0

    for curpt = 0, numofpts(velocity) - 1 do begin

      if abs(single(velocity,curpt)) le 1000.0 then 			      $
	colval= 100.0 + 100.0*single(velocity,curpt)/1000.0 		      $
      else								      $
	colval = 100.0+99.0*single(velocity,curpt)/abs(single(velocity,curpt))

      if single(velocity,curpt) ne -1e10 then 				      $
	oplot, [t(altitude,curpt)], [single(altitude, curpt)], 		      $
	       psym = 8, color = colval

      if (single(velocity,curpt) ne -1e10) and 				      $
	 (single(altitude,curpt) lt maxalt) then begin

        if (found eq 0) and (crb eq 1) and (curpt gt 0) then begin

	  sign = single(velocity, curpt)*single(velocity, curpt-1)

	  if (sign lt 0) and 						      $
	     (single(altitude, curpt) ne single(altitude, curpt-1)) then begin

	    found = 1

	    if crb_loc(0) eq 0 then begin
	      crb_tim(0) = t(altitude, curpt)
	      crb_loc(0) = single(altitude, curpt)
	    endif else begin 
	      crb_tim = [crb_tim,t(altitude, curpt)]
	      crb_loc = [crb_loc,single(altitude, curpt)]
            endelse

          endif

	endif

        if (curpt gt 0) and (found eq 1) and (crb eq 1) then		      $ 
          if (t(altitude, curpt) ne t(altitude, curpt-1)) then		      $
	    found = 0

        if (curpt gt 0) and (crb eq 1) then				      $
          if (t(altitude, curpt) ne t(altitude, curpt-1)) then begin
	    sign = single(velocity, curpt)*single(velocity, 0)
	    if sign lt 0 then begin

	      found = 1

	      if crb_loc(0) eq 0 then begin
	        crb_tim(0) = t(altitude, curpt)
	        crb_loc(0) = single(altitude, curpt)
	      endif else begin 
	        crb_tim = [crb_tim,t(altitude, curpt)]
	        crb_loc = [crb_loc,single(altitude, curpt)]
              endelse

            endif

          endif

        oplot, [t(altitude,curpt), 					      $
		t(altitude,curpt)+timadd*single(velocity, curpt)*factim],     $
               [single(altitude, curpt),				      $
		single(altitude, curpt)+single(velocity, curpt)*facvel],      $
               color = colval

        if oldaz ne round(single(azimuth,curpt)) then begin

	  oplot, [t(altitude,curpt),t(altitude,curpt)], [minalt,maxalt]
	  xyouts, t(altitude,curpt)+10.0, minalt+divfac/azdiv, 'Az = ' +      $
		  strcompress(string(round(single(azimuth,curpt))),	      $
			      /remove_all) + '!Eo',  			      $
	          charsize=cscale
	  oldaz=round(single(azimuth,curpt))

	  if azdiv eq 5.0 then azdiv = 10.0 else azdiv = 5.0

        endif

        if oldel ne round(single(elevation,curpt)) then begin

	  oplot, [t(altitude,curpt),t(altitude,curpt)], [minalt,maxalt]
	  xyouts, t(altitude,curpt)+10.0, minalt+divfac/2.0, 'El = '+	      $
	    	  strcompress(string(round(single(elevation,curpt))),	      $
			      /remove_all) + '!Eo', $
		  charsize=cscale

	  oldel=round(single(elevation,curpt))

        endif

      endif

    endfor

    mid = btr + 1.0*(etr - btr) / 8.0

    oplot, [mid+timadd], [minalt+divfac/6.0], psym = 8

    oplot, [mid+timadd, mid+timadd+timadd*1000.0*factim], 		      $
	   [minalt+divfac/6.0,minalt+divfac/6.0+1000.0*facvel]

    xyouts, mid+timadd+timadd*2000.0*factim,				      $
	    minalt+divfac/6.0+500.0*facvel, '1.0 Km/s', charsize=cscale

    if (crb eq 1) then if (n_elements(crb_loc) gt 1) then begin		      $
      nadd = 5
      rave, crb_loc, nadd
      if plotwhere eq 1 then oplot, crb_tim, crb_loc, thick = 5               $
      else oplot, crb_tim, crb_loc
    endif

  endif

  if final eq 1 then begin

      chartime=''
      spawn,'show time',chartime

      side  = 'date printed : '+chartime+' from file(s) '+datafile(0)
      xyouts, 0, 0, side, orientation=90, charsize=0.5*cscale, /device

      xyouts, 0.5,1.01, titleofplot, charsize=2.0*tscale, /norm, 	      $
	alignment=0.5

  endif  

  RETURN

END

