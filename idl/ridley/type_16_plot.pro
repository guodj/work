;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;	       DMSP Horizontal Velocities
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO Type_16_Plot, single, t, dt, times, moncheck, seconds, same, basetime, titlenam, unitnam,    $
                  placement, numofpts, plotwhere, datafile, titleofplot, maxran, mullen

  !p.multi=[0,3,2]

  !p.charsize = 2.0

  a = findgen(30) * (!pi*2/29.)

  cirsiz = ((1200.0/(6380.0*2.0*!pi))*20.0)*360.0/(90.0-maxran)

  usersym, 0.5*cos(a), 0.5*sin(a),/fill

  rad = placement(0,0)
  azm = placement(0,1)
  val = placement(0,2)
  rad2 = placement(0,3)

  minlat = maxran

  place = 1

  curpt = 1

  tickname = [' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']

  if curpt lt numofpts(rad) then begin

    if single(rad2,curpt) ne 0.0 then sgn = abs(single(rad2,curpt))/single(rad2,curpt) else sgn = 1.0

    if sgn eq 1.0 then place = 1 else place = 2

  endif

  xstart = 0.02
  dx     = (1.0 - 2.0*xstart)/3.5
  ystart = 0.05

  if place eq 2 then begin

    posx1 = xstart
    posy1 = 0.45 + ystart
    posx2 = posx1 + dx
    posy2 = posy1 + dx * (11.0/8.5) * 8.9/8.1

    plot, [-45.0],[-45.0], pos = [posx1,posy1,posx2,posy2], xrange=[-45.0,45.0], yrange=[-45.0,45.0], $
          xtickname = tickname, ytickname = tickname, xstyle = 5, ystyle = 5, title = 'empty', charsize=1.5

    x = findgen(361)*!pi/180.0 

    oplot, 40.0*sin(x), 40.0*cos(x)

    if (maxran eq 60) then begin

      oplot, (2./3.)*40.*sin(x),(2./3.)*40.*cos(x), linestyle = 1, thick = 0.5
      oplot, (1./3.)*40.*sin(x),(1./3.)*40.*cos(x), linestyle = 1, thick = 0.5

    endif

    oplot, [0,0], [-40.0,40.0]
    oplot, [-40.0,40.0], [0,0]

    pole = 'North Pole' 
    ml = '60!Uo'

    xyouts, 47.0, 46.0, pole, alignment=1.0
    xyouts, 30.0, 30.0, ml
    xyouts, -3.0,-45, '00'
    xyouts, -3.0, 41.0, '12'
    xyouts, 41.0, -2.0, '06'
    xyouts,-47.0, -2.0, '18'

  endif

  if max(single(azm,*)) lt 25.0 then begin

    for i=0,numofpts(azm)-1 do begin

      single(azm,i) = 360.0 * single(azm,i) / 24.0 + 180.0

      if single(azm,i) gt 360.0 then single(azm,i) = single(azm,i) - 360.

      single(azm,i) = 360.0 - single(azm,i)

    endfor

  endif

  while (numofpts(val)-1 gt curpt) do begin

    if place eq 1 then begin

      posx1 = xstart
      posy1 = 0.45 + ystart

    endif

    if place eq 3 then begin

      posx1 = xstart + dx + 0.07
      posy1 = 0.45 + ystart

    endif

    if place eq 5 then begin

      posx1 = xstart + 2.0*dx + 0.14
      posy1 = 0.45 + ystart

    endif

    if place eq 2 then begin

      posx1 = xstart
      posy1 = ystart

    endif

    if place eq 4 then begin

      posx1 = xstart + dx + 0.07
      posy1 = ystart

    endif

    if place eq 6 then begin

      posx1 = xstart + 2.0*dx + 0.14
      posy1 = ystart
      place = 0

    endif

    place = place + 1
    posx2 = posx1 + dx
    posy2 = posy1 + dx * (11.0/8.5) * 8.9/8.1

    mulfac = 8.0 * mullen

    if single(rad2,curpt) ne 0.0 then sgn = abs(single(rad2,curpt))/single(rad2,curpt) else sgn = 1.0

    sgnold = sgn

    if (place le 3) and (place ne 1) then    $
       datetl  = strmid(moncheck,(times(0,0,1)-1)*3,3)+' '+strcompress(seconds(1+times(0,0,2))+ ', 19'+    $
                 strcompress(string(times(0,0,0)),/remove_all))

    plot, [-45.0],[-45.0], pos = [posx1,posy1,posx2,posy2], xrange=[-45.0,45.0], yrange=[-45.0,45.0], $
          xtickname = tickname, ytickname = tickname, xstyle = 5, ystyle = 5, charsize=1.5

    x = findgen(361)*!pi/180.0 

    oplot, 40.0*sin(x), 40.0*cos(x)

    if maxran eq 60 then begin

      oplot, (2./3.)*40.*sin(x),(2./3.)*40.*cos(x), linestyle = 1, thick=0.5
      oplot, (1./3.)*40.*sin(x),(1./3.)*40.*cos(x), linestyle = 1, thick=0.5

    endif

    oplot, [0,0], [-40.0,40.0]
    oplot, [-40.0,40.0], [0,0]

    if single(rad2,curpt) ge 0.0 then begin

      pole = 'North Pole' 
      ml = '60!Uo'

    endif else begin

      pole = 'South Pole'
      ml = '-60!Uo'

    endelse

    xyouts, 47.0, 46.0, pole, alignment=1.0
    xyouts, 30.0, 30.0, ml
    xyouts, 0.0,-45, '00 MLT', alignment=0.5
    xyouts, 0.0, 41.0, '12', alignment=0.5
    xyouts, 47.0, -1.75, '06', alignment=1.0
    xyouts,-47.0, -1.75, '18'

    firstterm = 1
    distold = 1.0e32

    while (sgnold eq sgn) and (curpt lt numofpts(val)-2) do begin

      starttm = double(basetime(rad)) + double(t(rad,curpt-1))

      c_r_to_a, timearray, starttm

      plon = float(timearray(3)) - 2.2 + float(timearray(4))/60.0 + float(timearray(5))/3600.0
      plon = (plon / 24.) * 2.0 * !pi + !pi
      if plon gt 2.*!pi then plon = plon - 2.*!pi
      plon = 2.0*!pi - plon
      px = (90.0 - 73.6)*sin(plon)
      py = (90.0 - 73.6)*cos(plon)
      nlat = single(rad2,curpt)
      nlon = single(azm,curpt) * !pi / 180.0
      nx = (90.0 - abs(nlat)) * sin(nlon)
      ny = (90.0 - abs(nlat)) * cos(nlon)

      dist = sqrt((nx-px)^2.0+(ny-py)^2.0)

      if (firstterm eq 1) and (dist gt distold) then begin

        inserttl = strcompress(seconds(1+timearray(3)) +                                         $
                   ' :'+seconds(1+timearray(4)) + ' :'+seconds(1+timearray(5))+ ' UT')

	xyouts, 0.0, 46.0, inserttl, alignment=0.5

        if single(rad2,curpt) ge 0.0 then linsty = 0 else linsty = 1

        plat = 73.6
        prad = 40.0 * (90.0 - abs(plat)) / (90.0-minlat)
        px = prad * sin(plon)
        py = prad * cos(plon)

        oplot, [px], [py], psym = 8

	oplot, px+[cirsiz*cos(a)], py+[cirsiz*sin(a)], linestyle = linsty

        xyouts, px+1.0, py+1.0, 'STF', charsize = 0.5, alignment=0.5	

        firstterm = 0

      endif

      if (single(rad,curpt) gt 50.0) then distold = dist

      plat = single(rad,curpt-1)
      plon = single(azm,curpt-1) * !pi / 180.0
      prad = 40.0 * (90.0 - abs(plat)) / (90.0-minlat)
      px = prad * sin(plon)
      py = prad * cos(plon)

      qlat = single(rad,curpt)
      qlon = single(azm,curpt) * !pi / 180.0
      qrad = 40.0 * (90.0 - abs(qlat)) / (90.0-minlat)
      qx = qrad * sin(qlon)
      qy = qrad * cos(qlon)

      pqx = qx - px
      pqy = qy - py
      dpq = (pqx^2.0 + pqy^2.0)^0.5
      if dpq eq 0.0 then dpq = 0.0001

      nx = px - sgn * single(val,curpt-1) * mulfac * pqy / dpq
      ny = py + sgn * single(val,curpt-1) * mulfac * pqx / dpq

      if (abs(ny) le 40.0) and (abs(nx) le 40.0) and (abs(plat) ge (minlat)) then oplot, [px,nx], [py,ny]

      if abs(plat) ge (minlat) then oplot, [px,qx], [py,qy]

      curpt = curpt + 1

      if single(rad2,curpt) ne 0.0 then sgn = abs(single(rad2,curpt))/single(rad2,curpt) else sgn = 1.0

    endwhile

    px = -40.0
    py = -40.0

    nx = px
    ny = py + 1.0 * mulfac

    oplot, [px,nx],[py,ny]

    xyouts, px + 2.0, py, '1.0 Km/s'

    curpt = curpt + 1

    if (curpt eq numofpts(val)-1) or (place eq 1) then begin

      xyouts, 0.5, 0.95, datetl, charsize=1.5, /norm, alignment = 0.5
      xyouts, 0.5, 0, 'Minimum Latitude : '+strcompress(string(round(minlat)),/remove_all)+' Degrees', $
                      charsize=1.5, /norm, alignment=0.5

      chartime=''
      spawn,'show time',chartime

      side  = 'date printed : '+chartime+' from file(s) '+datafile(0)+'   '+datafile(1)
      side1 = datafile(2)+'   '+datafile(3)
      side2 = datafile(4)+'   '+datafile(5)
      xyouts, 0, 0, side, orientation=90, charsize=0.5, /norm
      xyouts, 0.02, 0.02, side1, orientation=90, charsize=0.5, /norm
      xyouts, 0.04, 0.02, side2, orientation=90, charsize=0.5, /norm

      xyouts, 0.5, 1.01, titleofplot, charsize=2.0, /norm, alignment=0.5

      if (plotwhere ne 1) and (numofpts(val)-1 gt curpt) then prompt_for_next

    endif

  endwhile

  !p.multi = 0

  !p.charsize = 1.0

  RETURN

END


