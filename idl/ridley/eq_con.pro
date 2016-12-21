;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;               Magnetometer Vector Plots
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO eq_con, single, t, times, moncheck, seconds, basetime, titlenam,    $
                  unitnam, numofpts, placement, mullen, plotwhere, lat, lon,  $
		  velocity, numofvar, stat, datafile, titleofplot, 	      $
		  alpha, beam, percent, opos, fontsize, titlesize

  if n_elements(percent) eq 0 then percent = 1.0
  cscale = percent*fontsize
  tscale = percent*titlesize
  relxlab = strarr(30)
  relxtl  = ' '
  btr       = 0.0
  etr       = 0.0
  curvar    = 0
  nl        = 0

  Compute_Axis, times, placement, btr, etr, curvar, relxlab, relxtl, 	      $
	basetime, moncheck, nl, ticktime, nminor

  if (etr gt max(t)) and (max(t) ne 0.0) then etr=max(t)

  if velocity ne 0.0 then tdis = lon / abs(velocity) else tdis = lon * 0.0

  relxtl  = relxtl + ' ' + 	      $
	' ;  Velocity = '+strcompress(string(velocity),/remove_all)

  relxlab = [' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',    $
             ' ',' ',' ',' ',' ',' ',' ']

  if min(tdis) lt 0.0 then btr = btr + min(tdis)
  if max(tdis) gt 0.0 then etr = etr + max(tdis)

  factor = (etr - btr) / 50.0

  btr = btr - 2*factor
  etr = etr + 2*factor

  relytl = 'Northward Displacement in Km'

  ymin = min(lat)
  ymax = max(lat)

  ymin = 50
  ymax = 90

  xoff = 0.0
  yoff = 0.0

  if ymin lt -200 then scalept = -375
  if ymin lt -500 then scalept = -750

  xpos1 = 0.1
  xpos2 = 0.95
  ypos1 = 0.1
  ypos2 = 0.95

  xyratio = 10.5 / 7.5

  yrange = ymax - ymin

  alp2 = -1.0 * (-1.0*alpha + 90.0) * !pi / 180.0
;  alp2 = -1.0 * (-1.0*alpha) * !pi / 180.0

  plot, [btr, btr], [ymin, ymin],              			$
        xrange = [btr, etr],                   			$
        yrange = [ymin, ymax],					$
        xticks = nl,                           			$
        xstyle = 1,                            			$
        ystyle = 1,                            			$
        xtickname = relxlab,					$
        xtitle = relxtl,					$
        ytitle = relytl,					$
	position = [xpos1,ypos1,xpos2,ypos2],			$
        charsize = 1.5*cscale,xtickv=ticktime

  ppos, perpos

  xyratio = (perpos(2)-perpos(0))/(perpos(3)-perpos(1))
  xdenorm = (etr - btr) / 13.3*1.25
  ydenorm = (yrange / (etr - btr)) * xdenorm * xyratio

  opos = [xpos1,ypos1,xpos2,ypos2] 

  a = findgen(16) * (!pi*2/16.)

  usersym, 0.25*cos(a), 0.25*sin(a), /fill

  if abs(beam) lt 360.0 then azim = (beam+17.0)*!pi/180.0

  for i=0, numofvar/3 do begin

    if (numofpts(i*3) gt 0) and (numofpts(i*3+1) gt 0) then begin

      for j=0, numofpts(i*3)-1 do begin

        y1 = lat(i)
        if velocity lt 0.0 then x1 = t(i*3,j) + tdis(i) else x1 = t(i*3,numofpts(i*3)-j-1) + tdis(i)

        x2dum = single(i*3,j) / 200.0
        y2dum = single(i*3+1,j) / 200.0

        if abs(beam) gt 360.0 then begin

	  x2 = x1 + 				$
	    (x2dum * cos(alp2(i)) + 		$
	     y2dum * sin(alp2(i)) - yoff/200.0)*mullen*xdenorm
	  y2 = y1 + 				$
	    (-1.0 * x2dum * sin(alp2(i)) + 	$
	    y2dum * cos(alp2(i)) + xoff/200.0)*mullen*ydenorm

        endif else begin

          x2 = (x2dum * cos(alp2(i)) + y2dum * sin(alp2(i)) - yoff/200.0)
	  y2 = (-1.0 * x2dum * sin(alp2(i)) + y2dum * cos(alp2(i)) + xoff/200.0)

          r2 = sqrt(x2^2.0 + y2^2.0)
          if x2 ne 0.0 then the2 = atan(y2/x2) else the2 = 90.0
          if x2 lt 0.0 then the2 = !pi + the2

	  y2 = y1 + r2*cos(the2 + azim - !pi/2.0)*mullen*ydenorm 
	  x2 = x1 + r2*cos(the2 + azim - !pi/2.0)*mullen*xdenorm*sin(10.0*!pi/180.0)

	endelse

        if (single(i*3,j) ne -1e10) and (single(i*3+1,j) ne -1e10) then begin
	  oplot, [x1], [y1], psym=8
	  oplot, [x1,x2], [y1,y2]
        endif

        if (j eq 0) and (i eq 0) then begin

	  xyouts, x1+factor/2.0, ymin - abs((ymax-ymin)/20.0),		      $
		'First', charsize=cscale

	endif

        if (j eq numofpts(i*3)-1) and (i eq 0) then begin

	  xyouts, x1+factor/2.0, ymin - abs((ymax-ymin)/20.0), 		      $
		'Last', charsize=cscale

	endif

      endfor

      xyouts, btr+factor, lat(i), stat(i), charsize=cscale

    endif

  endfor

   y1 = ymin + yrange/10.0
   x1 = btr+10.0*factor

   x2 = x1 + xdenorm
   y2 = y1

   oplot, [x1], [y1], psym=8
   oplot, [x1,x2], [y1,y2]

   xyouts, x2+2.0*factor, y1, strcompress(string(200.0/mullen)+' nT'),	      $
	charsize=cscale

  chartime=''
  spawn,'date',chartime

      side  = 'date printed : '+chartime+' from file(s) '+datafile(0)+'   '+datafile(1)
      side1 = datafile(2)+'   '+datafile(3)
      side2 = datafile(4)+'   '+datafile(5)
      xyouts, 0, 0, side, orientation=90, charsize=0.5*cscale, /device
      xyouts, 200, 2000, side1, orientation=90, charsize=0.5*cscale, /device
      xyouts, 400, 2000, side2, orientation=90, charsize=0.5*cscale, /device

  xyouts, 0.5, 1.01, titleofplot, charsize=2.0*tscale, /norm, alignment=0.5

  RETURN

END

