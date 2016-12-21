
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;	       Polar Contour
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PRO Type_5_Plot, single, t, dt, times, moncheck, seconds, same, basetime, rows, columns, numofgra,    $
		 titlenam, unitnam, placement, maxminque, numofpts, plotwhere, datafile, titleofplot, $
		 rotation

  !p.multi=0

  rad = placement(0,0)
  azm = placement(0,1)
  val = placement(0,2)
  ele = placement(0,3)
  eve = placement(0,4)

  evecha = intarr(100)

  x = where(single(rad,*) eq max(single(rad,*)))
  y = where(single(azm,*) eq max(single(azm,*)))

  nran = x(0) + 1
  naz = y(0) / nran + 1

  latitude = fltarr(nran + 1)
  longitude = fltarr(naz + 1)

  evevalues = single(eve, 0:numofpts(eve)-1)

  numofeve = 1
  evecha(numofeve-1) = 0
  numofeve = 2

  for i=1, numofpts(eve)-2 do if (evevalues(i) ne evevalues(i-1)) then begin

    evecha(numofeve-1) = i
    numofeve = numofeve + 1

  endif

  evecha(numofeve-1) = numofpts(eve)-1

  numofeve = numofeve - 1

  for cureve = 0, numofeve-1 do begin

    polear = fltarr(naz+1,nran)
    xcor   = fltarr(naz+1,nran)
    ycor   = fltarr(naz+1,nran)
    latts  = fltarr(naz+1,nran)

    radvalues = single(rad, evecha(cureve):evecha(cureve+1)-1)
    azmvalues = single(azm, evecha(cureve):evecha(cureve+1)-1)
    valvalues = single(val, evecha(cureve):evecha(cureve+1)-1)

    maxrad = max(radvalues)

    raddiv = nran / maxrad
    azmdiv = naz / 360.0

    lastrad = nran

    for curnum = 0, n_elements(valvalues)-1 do begin

      if (valvalues(curnum) ne -1.0e10) then begin

        catch = 0
        azdum = floor((azmvalues(curnum) + rotation)*azmdiv)
        if azdum gt naz - 1 then azdum = azdum - naz
	raddum = round(radvalues(curnum)*raddiv)
        polear(azdum,raddum) = valvalues(curnum) 

      endif else begin

        if catch eq 0 then begin

          catch = 1
          lastrad = round(radvalues(curnum)*raddiv) - 1

        endif

      endelse

    endfor

    polear(naz,*) = polear(0,*)

    for j=0,naz do begin

      longitude(j) = (float(j)*2.0*!pi/float(naz))

    endfor

    for i=0,lastrad do begin

       latts(*,i) = float(i)
       xcor(*,i) = latts(*,i)*sin(longitude)
       ycor(*,i) = latts(*,i)*cos(longitude)

    endfor

    xcor = xcor(*,0:lastrad)
    ycor = ycor(*,0:lastrad)
    polear = polear(*,0:lastrad)

    posx1 = 0.15
    posy1 = 0.0
    posx2 = posx1 + (8.5/11.0) * 8.1/8.9
    posy2 = posy1 + 1.0

    tickname = [' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']

    divlev = round((max(polear) - min(polear))/29.0)

    lev = fltarr(30)

    for i=0,29 do lev(i) = -75 + i * 5

    starttm = double(basetime(rad)) + double(t(rad,evecha(cureve)+1))

    c_r_to_a, timearray, starttm

    times(0,0,0) = timearray(0)
    times(0,0,1) = timearray(1)
    times(0,0,2) = timearray(2)
    times(0,0,3) = timearray(3)
    times(0,0,4) = timearray(4)
    times(0,0,5) = timearray(5)

    datetl  = strcompress(seconds(1+times(0,0,1))+' /'+seconds(1+times(0,0,2))+ ' /'+string(times(0,0,0)))

    inserttl = strcompress(seconds(1+times(0,0,3)) +                                         $
                           ' :'+seconds(1+times(0,0,4)) + ' :'+seconds(1+times(0,0,5))+ ' UT')

    contour, polear, xcor, ycor, nlevels=30, /follow, pos = [posx1,posy1,posx2,posy2],       $
           xtickname = tickname, ytickname = tickname, c_linestyle = (lev lt 0.0),           $
           levels = lev, xrange=[-40.0,40.0], yrange=[-40.0,40.0]

    x = findgen(360)*!pi/180.0 

    oplot, 40.0*sin(x), 40.0*cos(x)

    oplot, [0,0], [-40.0,40.0]
    oplot, [-40.0,40.0], [0,0]

    xyouts, -38.0, 35.0, datetl
    xyouts, 20.0, 35.0, inserttl
    xyouts, -38.0, -38.0, 'Maximum '+unitnam(rad)+' '+strcompress(string(round(maxrad)),/remove_all)+' Km'

    disfac = 40.0/500.0

    elev = 30.0

;    Put_Stations, disfac, elev

    chartime=''
    spawn,'show time',chartime

    side  = 'date printed : '+chartime+' from file(s) '+datafile(0)+'   '+datafile(1)
    side1 = datafile(2)+'   '+datafile(3)
    side2 = datafile(4)+'   '+datafile(5)
    xyouts, 0, 0, side, orientation=90, charsize=0.5, /device
    xyouts, 200, 2000, side1, orientation=90, charsize=0.5, /device
    xyouts, 400, 2000, side2, orientation=90, charsize=0.5, /device

    xyouts, 0.5, 1.01, titleofplot, charsize=2.0, /norm, alignment=0.5

    if (cureve ne numofeve-1) and (plotwhere ne 1) then prompt_for_next

  endfor

  RETURN

END
