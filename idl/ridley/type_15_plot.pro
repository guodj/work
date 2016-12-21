;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;	       Special Polar Contours for Sondrestrom
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO Type_15_Plot, single, t, dt, times, moncheck, seconds, same, basetime,         $
                  rows, columns, numofgra, titlenam, unitnam, placement, maxminque, numofpts, plotwhere, datafile, titleofplot

  !p.multi=[0,4,2]

  rad = placement(0,0)
  azm = placement(0,1)
  val = placement(0,2)
  ele = placement(0,3)
  eve = placement(0,4)

  evecha = intarr(100)

  latitude = fltarr(41)
  longitude = fltarr(73)

  evevalues = single(eve, 0:numofpts(eve)-1)

  numofeve = 1
  evecha(numofeve-1) = -1
  numofeve = 2

  for i=1, numofpts(eve)-1 do if (evevalues(i) ne evevalues(i-1)) then begin

    evecha(numofeve-1) = i-1
    numofeve = numofeve + 1

  endif

  evecha(numofeve-1) = numofpts(eve)-1

  numofeve = numofeve - 1

  place = 1

  for cureve = 0, numofeve-1 do begin

    polear = fltarr(73,41)
    xcor   = fltarr(73,41)
    ycor   = fltarr(73,41)
    latts  = fltarr(73,41)

    radvalues = single(rad, evecha(cureve)+1:evecha(cureve+1))
    azmvalues = single(azm, evecha(cureve)+1:evecha(cureve+1))
    valvalues = single(val, evecha(cureve)+1:evecha(cureve+1))

    maxrad = max(radvalues)

    raddiv = 40.0 / maxrad
    azmdiv = 72.0 / 360.0

    lastrad = 40

    for curnum = 0, n_elements(valvalues)-1 do begin

      if (valvalues(curnum) ne -1e10) then begin

        catch = 0

	azmvalues(curnum) = azmvalues(curnum) + 197.0
	if azmvalues(curnum) gt 360.0 then azmvalues(curnum) = azmvalues(curnum) - 360.0
        polear(azmvalues(curnum)*azmdiv,(radvalues(curnum))*raddiv) = valvalues(curnum)

      endif else begin

        if catch eq 0 then begin

          catch = 1
          lastrad = round(radvalues(curnum)*raddiv) - 1

        endif

      endelse

    endfor

    polear(72,*) = polear(0,*)

    for j=0,lastrad do latitude(j) = -1*float(j)*50.0/40.0 + 90
  
    for j=0,72 do longitude(j) = (float(j)*2*!pi/72.0)

    for i=0,lastrad do begin

       latts(*,i) = (90.0 - latitude(i)) * 0.8
       xcor(*,i) = latts(*,i)*sin(longitude)
       ycor(*,i) = latts(*,i)*cos(longitude)

    endfor

    xcor = xcor(*,0:lastrad)
    ycor = ycor(*,0:lastrad)
    polear = polear(*,0:lastrad)

    xstart = 0.01
    dx     = (1.0 - 2.0*xstart)/4.1
    ystart = 0.075

    if place eq 1 then begin

      posx1 = xstart
      posy1 = 0.45 + ystart

    endif

    if place eq 2 then begin

      posx1 = xstart + dx + 0.01
      posy1 = 0.45 + ystart

    endif

    if place eq 3 then begin

      posx1 = xstart + 2.0*dx + 0.02
      posy1 = 0.45 + ystart

    endif

    if place eq 4 then begin

      posx1 = xstart + 3.0*dx + 0.03
      posy1 = 0.45 + ystart

    endif

    if place eq 5 then begin

      posx1 = xstart
      posy1 = ystart

    endif

    if place eq 6 then begin

      posx1 = xstart + dx + 0.01
      posy1 = ystart

    endif

    if place eq 7 then begin

      posx1 = xstart + 2.0*dx + 0.02 
      posy1 = ystart

    endif

    if place eq 8 then begin

      posx1 = xstart + 3.0*dx + 0.03
      posy1 = ystart
      place = 0

    endif

    place = place + 1
    posx2 = posx1 + dx
    posy2 = posy1 + dx * (11.0/8.5)

    tickname = [' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']

    divlev = round((max(polear) - min(polear))/29.0)

    lev = fltarr(21)

    for i=0,20 do lev(i) = -60 + 6*i

    starttm = double(basetime(rad)) + double(t(rad,evecha(cureve)+1))

    c_r_to_a, timearray, starttm

    datetl  = strmid(moncheck,(timearray(1)-1)*3,3)+' '+strcompress(seconds(1+timearray(2))+ ', 19'+    $
               strcompress(string(timearray(0)),/remove_all))

    inserttl = strcompress(seconds(1+timearray(3)) +                                         $
                           ' :'+seconds(1+timearray(4)) + ' :'+seconds(1+timearray(5))+ ' UT')

    plot, [-40],[-40], xrange=[-40.0,40.0], yrange=[-40.0,40.0], title = inserttl,           $
          xtickname = tickname, ytickname = tickname, xstyle = 5, ystyle = 5,                $
          pos = [posx1,posy1,posx2,posy2], charsize = 2.0

    contour, polear, xcor, ycor, nlevels=21, /follow, levels = lev, /overplot, $
             charsize = 1.25, c_labels = [0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0] 

    x = findgen(361)*!pi/180.0 

    oplot, 40.0*sin(x), 40.0*cos(x)

    oplot, [0,0], [-40.0,40.0]
    oplot, [-40.0,40.0], [0,0]

    if (cureve eq numofeve-1) or (place eq 1) then begin

      xyouts, 0.5, 0.95, datetl, charsize=1.5, /norm, alignment=0.5
      xyouts, 0.5, 0.0, 'Maximum Range : '+strcompress(string(maxrad),/remove_all)+' Km',   $
		        charsize=1.5, /norm, alignment=0.5

      chartime=''
      spawn,'show time',chartime

      side  = 'date printed : '+chartime+' from file(s) '+datafile(0)+'   '+datafile(1)
      side1 = datafile(2)+'   '+datafile(3)
      side2 = datafile(4)+'   '+datafile(5)
      xyouts, 0, 0, side, orientation=90, charsize=0.5, /device
      xyouts, 200, 2000, side1, orientation=90, charsize=0.5, /device
      xyouts, 400, 2000, side2, orientation=90, charsize=0.5, /device

      xyouts, 0.5, 1.01, titleofplot, charsize=2.0, /norm, alignment=0.5

      if (plotwhere ne 1) and (cureve ne numofeve-1) then prompt_for_next

    endif

  endfor

  !p.multi = 0

  !p.charsize = 1.0

  RETURN

END


