;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;	       Eight Polar Line of Sight Velocity Plots
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO Type_14_Plot, single, t, times, moncheck, seconds, basetime, titlenam, unitnam, numofpts, placement,   $
                  mullen, maxran, plotwhere, survey, datafile, titleofplot, checkinvar

  red     = intarr(256)
  blue    = intarr(256)
  green   = intarr(256)  
  relxlab = strarr(30)

  curvar    = 0
  altitude  = 0
  azimuth   = 1
  velocity  = 2
  elevation = 3
  event     = 4

  sizeoftype = 1.5

  minalt = 0.0
  maxalt = maxran

  ytlr =  ' (' + strcompress(unitnam(altitude), /remove_all) + ')'

  ytl  = strcompress(titlenam(altitude), /remove_all) + ' ' + ytlr

  !p.multi=[0,4,2]

    red(0:2)      = [0,0,255]
    red(3:255)    = 255
    green(0:2)    = [0,0,0]
    green(3:255)  = 255
    blue(0:2)     = [0,255,0]
    blue(3:255)   = 255

  tvlct, red, green, blue

  numofeve = 0

  if numofpts(event) gt 1 then begin

    for i = 1, numofpts(event)-1 do if single(event,i) ne single(event,i-1) then numofeve = numofeve + 1

  endif

  startpt = 0
  dummy   = 0
  starttm = double(0.0)
  place = 1

  for cureve = 0, numofeve do begin

    starttm = double(basetime(event)) + double(t(event,startpt))

    c_r_to_a, timearray, starttm

    if place eq 1 then datetl  = strmid(moncheck,(timearray(1)-1)*3,3)+' '+strcompress(seconds(1+timearray(2))+ ', 19'+    $
               strcompress(string(timearray(0)),/remove_all))

    inserttl = strcompress(seconds(1+timearray(3)) + ' :'+seconds(1+timearray(4)) + ' UT')

    elev = single(elevation,startpt)

    if elev gt 59.0 then elev=59.0

    insertel = strcompress('E='+seconds(1+round(elev)))

    xstart = 0.02
    dx     = (1.0 - 2.0*xstart)/4.0
    ystart = 0.1

    if place eq 1 then begin

      posx1 = xstart
      posy1 = 0.4 + ystart

    endif

    if place eq 2 then begin

      posx1 = xstart + dx
      posy1 = 0.4 + ystart

    endif

    if place eq 3 then begin

      posx1 = xstart + 2.0*dx
      posy1 = 0.4 + ystart

    endif

    if place eq 4 then begin

      posx1 = xstart + 3.0*dx 
      posy1 = 0.4 + ystart

    endif

    if place eq 5 then begin

      posx1 = xstart
      posy1 = ystart

    endif

    if place eq 6 then begin

      posx1 = xstart + dx 
      posy1 = ystart

    endif

    if place eq 7 then begin

      posx1 = xstart + 2.0*dx 
      posy1 = ystart

    endif

    if place eq 8 then begin

      posx1 = xstart + 3.0*dx 
      posy1 = ystart
      place = 0

    endif

    place = place + 1
    posx2 = posx1 + dx
    posy2 = posy1 + dx * (11.0/8.5) * 8.9/8.1

    a = findgen(16) * (!pi*2/16.)

    gentick = [' ',' ',' ']
    lowtick = [round(-1*maxalt),0,round(maxalt)]

    usersym, 0.3*cos(a), 0.3*sin(a), /fill

      plot, [-maxalt,-maxalt], [-maxalt,maxalt],                     $
            charsize  = sizeoftype,                                  $
            xra       = [-maxalt, maxalt],                           $
            yra       = [-maxalt, maxalt],		             $
	    xticks    = 2,                                           $
	    yticks    = 2,                                           $
            xtickname = gentick,                                     $
	    ytickname = gentick,				     $
            xminor    = maxalt/100,                                  $
            yminor    = maxalt/100,                                  $
            xstyle    = 1,                                           $
            ystyle    = 1,					     $
            title     = inserttl,				     $
            pos=[posx1, posy1, posx2, posy2]
                                             
    if numofpts(velocity)-startpt gt 0 then begin
  
      facvel = mullen * (maxalt-minalt) / 5000.0
      factim = 1.0 / 2000.0

      curpt  = startpt
      eveold = single(event,startpt)

      while (curpt le numofpts(velocity) - 1) and (eveold eq single(event,curpt)) do begin
                                                 
        if single(velocity,curpt) ge 0.0 then colval=2 else colval=1

        r   = single(altitude,curpt)
        t1  = single(azimuth,curpt) * 2.0 * !pi / 360.0
        v   = single(velocity,curpt)*facvel
        f   = 15.0 * 2.0 * !pi / 360.0
      
        s  = v * sin(f)
        p  = r + v * cos(f)
        sp = sqrt(s*s + p*p)
        if sp ne 0.0 then t2 = asin(s/sp) else t2=1.0
        if p lt 0.0 then t2 = !pi-t2
        if cos(t2) ne 0.0 then e  = p / cos(t2) else e=p

        x1 = r * sin(t1)
        x2 = e * sin(t1 + t2)

        y1 = r * cos(t1)
        y2 = e * cos(t1 + t2)

        if (single(velocity,curpt) ne -1e10) then  $
           oplot, [x1], [y1], psym = 8, color = colval

        if (single(velocity,curpt) ne -1e10) then  $
           oplot, [x1,x2],[y1,y2], color = colval

        eveold = single(event,curpt)
        curpt  = curpt + 1

      endwhile

      x1 = -maxalt + maxalt/20.0
      x2 = x1 + 1000.0*facvel/5.0
      y1 = x1
      y2 = y1 + 1000.0*facvel

      oplot, [x1], [y1], psym = 8

      oplot, [x1,x2], [y1,y2]

      xyouts, x2 + 1000.0*facvel/4.0, y1+500.0*facvel, '1.0 Km/s'
      xyouts, maxalt*4.0/7.0, y1+500.0*facvel, insertel

      if checkinvar eq 2 then begin

        posarr = [posx1, posy1, posx2, posy2]
        if titlenam(0) eq 'Range' then maxrange = maxalt*cos(single(elevation,0)*!pi/180.0) else maxrange = maxalt
        data_array = single(5:7,*)
        basic_contour, data_array, maxrange, posarr

      endif

      disfac = 1.0

      Put_Stations, disfac, elev

    endif

    startpt = curpt

    if (numofeve eq cureve) or (place eq 1) then begin

      xyouts, 0.5, 0.95, datetl, charsize=1.5, /norm, alignment=0.5
      xyouts, 0.5, 0.0, 'Maximum Range : '+strcompress(string(round(maxalt)),/remove_all)+' Km',   $
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
    
      if (plotwhere ne 1) and (numofeve ne cureve) then prompt_for_next

    endif

  endfor

  !p.multi = 0

  RETURN

END

