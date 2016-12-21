;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;               Magnetometer Latitude Compairison
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO Type_18_Plot, single, t, times, moncheck, seconds, basetime, titlenam, unitnam, numofpts, placement,  $
		  plotwhere, lat, lon, numofvar, stat, datafile, titleofplot, wvoff, mullen, alpha, rotation

  relxlab = strarr(30)
  relxtl  = ' '
  btr       = 0.0
  etr       = 0.0
  curvar    = 0
  nl        = 0

  Compute_Axis, times, placement, btr, etr, curvar, relxlab, relxtl, basetime, moncheck, nl, ticktime, nminor

  if (etr gt max(t)) and (max(t) ne 0.0) then etr=max(t) 

  al2 = (-1.0*alpha+rotation)*!pi/180.0

  factor = (etr - btr) / 50.0

  btr = btr - 5*factor
  etr = etr + 2*factor

  ymin = min(lat)
  ymax = max(lat)

  factory = (ymax - ymin) / 50.0

  xpos1 = 0.1
  xpos2 = 0.95
  ypos1 = 0.1
  ypos2 = 0.95

  relytl = 'Northward Distance from STF'

  plot, [btr, btr], [ymin, ymin],              $
        xrange = [btr, etr],                   $
        yrange = [ymin - 3.0*factory, ymax + 3.0*factory],                 $
        xticks = nl,                           $
        xstyle = 1,                            $
        ystyle = 1,                            $
        xtickname = relxlab,                   $
        xtitle = relxtl,                       $
        ytitle = relytl,                       $
	position = [xpos1,ypos1,xpos2,ypos2],  $
        charsize = 1.5,                        $
        xtickv=ticktime,                       $
	xminor = nminor

  a = findgen(16) * (!pi*2/16.)

  usersym, 0.25*cos(a), 0.25*sin(a), /fill

  maxpts = max(numofpts) - 1

  for i = 0, maxpts do begin

    oplot, [t(0,i), t(0,i)], [ymin, ymax], linestyle=1

    jold = -1

    for j=0,numofvar/3 do begin

      if (numofpts(j*3+wvoff-1)-1 ge i) and   $
         (abs(single(j*3+wvoff-1,i)) lt 2000.0) then begin

;	if wvoff eq 3 then begin

; 	  value1 = single(j*3+2,i)*factor*mullen/100.0
;	  if (jold gt -1) then value2 = single(jold*3+2,i)*factor*mullen/100.0
 	  value1 = single(j*3+wvoff-1,i)*factor*mullen/100.0
	  if (jold gt -1) then value2 = single(jold*3+wvoff-1,i)*factor*mullen/100.0

;	endif else begin

;	   xval1 = single(j*3,i)*factor*mullen/100.0
;	   yval1 = single(j*3+1,i)*factor*mullen/100.0

;	   if (jold gt -1) then begin

;	     xval2 = single(jold*3,i)*factor*mullen/100.0
;	     yval2 = single(jold*3+1,i)*factor*mullen/100.0

;	   endif

;	   if wvoff eq 1 then begin

;	     value1 = xval1 * cos(al2(j)) + yval1 * sin(al2(j))
;	     if (jold gt -1) then value2 = xval2 * cos(al2(jold)) + yval2 * sin(al2(jold))

;           endif else begin

;	     value1 = -1.0 * xval1 * sin(al2(j)) + yval1 * cos(al2(j))
;	     if (jold gt -1) then value2 = -1.0 * xval2 * sin(al2(jold)) + yval2 * cos(al2(jold))

;	   endelse

;	endelse

	if (jold eq -1) and (j gt 0) then oplot, [t(0,i) + value1, t(0,i)], [lat(j), ymax]

	if (jold gt -1) then oplot, [t(0,i) + value1, t(0,i) + value2], [lat(j), lat(jold)]

	jold = j

        if (i eq 0) then xyouts, -4.0*factor, lat(j), stat(j)

      endif

    endfor   

  endfor

  y1 = ymin - factory*1.5
  x1 = btr + 5.0*factor

  x2 = x1 + 2.0*factor
  y2 = y1

  oplot, [x1], [y1], psym=8
  oplot, [x1,x2], [y1,y2]

  if wvoff eq 1 then addon = 'East / West Variation'
  if wvoff eq 2 then addon = 'North / South Variation'
  if wvoff eq 3 then addon = 'Vertical Variation'

  xyouts, x2+factor, y1, strcompress(string(200.0/mullen)+' nT ')+addon

  chartime=''
  spawn,'show time',chartime

  side  = 'date printed : '+chartime+' from file(s) '+datafile(0)+'   '+datafile(1)
  xyouts, -200, 0, side, orientation=90, charsize=0.5, /device

  xyouts, 0.5, 1.01, titleofplot, charsize=2.0, /norm, alignment=0.5

  RETURN

END

