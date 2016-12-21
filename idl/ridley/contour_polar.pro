;
; procedure polar_contour
;
;  contours a 2 dimensional polar array in the specified position
;
;  data_array - fltarr(3,*) 	1 - azimuth values
;				2 - range values
;				3 - values to contour
;  maxrange - should be pretty close max(data_array(1,*))
;  posarr   - contains array specifying position to plot contour in /norm
;

PRO contour_polar, data_array, maxrange, posarr, levels=levels, 	$
	plotmaxmin = plotmaxmin, white = white, color=color

  if n_elements(plotmaxmin) eq 0 then pmm = 0				$
  else pmm = 1

  if n_elements(white) eq 0 then begin
    if !d.name eq 'PS' then pw = 0 else pw = 255
  endif else begin
    if !d.name eq 'PS' then pw = 255 else pw = 0
  endelse

  if n_elements(color) eq 0 then pcolor = 0                             $
  else pcolor = 1

  azm = 0
  rad = 1
  val = 2

  nele = n_elements(data_array(0,*))

  radvalues = data_array(rad, 0:nele-1)
  azmvalues = data_array(azm, 0:nele-1)
  valvalues = data_array(val, 0:nele-1)

  loc = where(radvalues eq radvalues(0), nazm)
  loc = where(azmvalues eq azmvalues(0), nrad)

  locr = where(radvalues(loc) le maxrange, nrad)

  maxrad = max([max(radvalues),maxrange])+1.0
  locr = where(radvalues eq 0, count)
  if count eq 0 then begin
    nrad = nrad + 1
    zero = 0
  endif else zero = 1

  polear = fltarr(nazm+1,nrad)
  xcor   = fltarr(nazm+1,nrad)
  ycor   = fltarr(nazm+1,nrad)

  for i = 0, nrad-1 do begin

    loc = where(radvalues eq min(radvalues))
    if (zero eq 0) and (i eq 0) then begin
      ave = mean(valvalues(loc))
      polear(0:nazm-1,i) = ave
      xcor(0:nazm-1,i) = 0.0
      ycor(0:nazm-1,i) = 0.0
    endif else begin
      if i eq 0 then begin
        ave = mean(valvalues(loc))
	polear(0:nazm-1,i) = ave
      endif else polear(0:nazm-1,i) = valvalues(loc)
      xcor(0:nazm-1,i) = radvalues(loc)*cos(azmvalues(loc)*!pi/180.0)
      ycor(0:nazm-1,i) = radvalues(loc)*sin(azmvalues(loc)*!pi/180.0)
      radvalues(loc) = maxrad
    endelse

  endfor
    
  xcor(nazm,*) = xcor(0,*)
  ycor(nazm,*) = ycor(0,*)
  polear(nazm,*) = polear(0,*)
        
  minval = round(min(polear))
  maxval = round(max(polear))

  tickname = [' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']

  if n_elements(levels) eq 0 then begin

    lev = fltarr(30)

    if (maxval lt 100.0) and (maxval-minval gt 20.0) then begin

      minval = -75
      divlev = 5

      for i=0,29 do lev(i) = minval + i * divlev

    endif else begin

      divlev = 2*round((maxval - minval)/29.0)
      for i=0,29 do lev(i) = minval + i * divlev

      if divlev eq 0 then begin

        divlev = 2.0
        for i=0,29 do lev(i) = minval + (float(i) - 0.5) * divlev

      endif

    endelse

  endif else begin

    if n_elements(levels) le 30 then lev = levels 			$
    else lev = levels(0:29)

  endelse

  nlev = n_elements(lev)

  plot, [-maxrange-1,maxrange+1],[-maxrange-1,maxrange+1], 		$
	/nodata, pos = posarr, xstyle=5, ystyle=5, /noerase

  contour, polear, xcor, ycor, nlevels=nlev, /follow, 			$
	c_linestyle = 3.0*(lev lt 0.0),		 			$
	levels = lev, /overplot, 					$
	charsize=1.0, max_value=maxval+1.0, color = pw

  if pmm eq 1 then begin
    clip = float(!p.clip)
    yoff = 0.5*maxrange*float(!d.y_ch_size)/(clip(3)-clip(1))
    loc = where(polear eq min(polear))
    xyouts, xcor(loc(0)), ycor(loc(0))-yoff, '-', alignment=0.5
    loc = where(polear eq max(polear))
    xyouts, xcor(loc(0)), ycor(loc(0))-yoff, '+', alignment=0.5
  endif

  RETURN

END

