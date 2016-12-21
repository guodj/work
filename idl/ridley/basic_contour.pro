;
; procedure basic_contour
;
;  contours a 2 dimensional polar array in the specified position
;
;  data_array - fltarr(3,*) 	1 - azimuth values
;				2 - range values
;				3 - values to contour
;  maxrange - should be pretty close max(data_array(1,*))
;  posarr   - contains array specifying position to plot contour in /norm
;

PRO basic_contour, data_array, maxrange, posarr, levels

  azm = 0
  rad = 1
  val = 2

  azmdum = where(data_array(azm,*) eq max(data_array(azm,*)), azmcount)
  raddum = where(data_array(rad,*) eq max(data_array(rad,*)), radcount)

  if raddum(0) lt azmdum(0) then begin

    nrad = raddum(0)+1
    nazm = (azmdum(azmcount-1)+1)/nrad

  endif else begin

    nazm = azmdum(0)+1
    nrad = (raddum(radcount-1)+1)/nazm

  endelse

  nele = where(data_array(rad,*) eq 0, count)

  if count eq 0 then nele = n_elements(data_array(rad,*)) else nele = nele(0)

  radvalues = data_array(rad, 0:nele-1)
  azmvalues = data_array(azm, 0:nele-1)
  valvalues = data_array(val, 0:nele-1)

  maxrad = max(radvalues)

  if maxrange ne maxrad then begin

    raddum = where(radvalues gt maxrange)

    if raddum(0) ne -1 then begin

       maxrad = radvalues(raddum(0)-1)
       nrad   = raddum(0)

    endif

  endif

  raddiv = float(nrad) / maxrad
  azmdiv = float(nazm) / 360.0

  latitude = fltarr(nrad+1)
  longitude = fltarr(nazm+1)
  polear = fltarr(nazm+1,nrad+1)
  xcor   = fltarr(nazm+1,nrad+1)
  ycor   = fltarr(nazm+1,nrad+1)
  latts  = fltarr(nazm+1,nrad+1)

  lastrad = nrad

  firstval = 0

  maxval = round(max(valvalues))

  for curnum = 0, n_elements(valvalues)-1 do begin

    if (valvalues(curnum) ne -1e10) and (radvalues(curnum) le maxrad) then begin

      catch = 0
      azdum = floor(azmvalues(curnum)*azmdiv)
      raddum = round(radvalues(curnum)*raddiv)
      polear(azdum,raddum) = valvalues(curnum) 

    endif else begin

      if catch eq 0 then begin

        catch = 1
        lastrad = round(radvalues(curnum)*raddiv) - 1

      endif

    endelse

  endfor

  polear(nazm,*) = polear(0,*)

  for j=0,nazm do begin

    longitude(j) = (float(j)*2.0*!pi/float(nazm))

  endfor

  for i=0,lastrad do begin

    latts(*,i) = (float(i))*maxrad/float(lastrad)
    xcor(*,i)  = latts(*,i)*sin(longitude)
    ycor(*,i)  = latts(*,i)*cos(longitude)

  endfor

  xcor = xcor(*,0:lastrad)
  ycor = ycor(*,0:lastrad)
  polear = polear(*,0:lastrad)

  sum = 0.0

  for i=0,nazm-1 do sum = sum + polear(i,1)

  polear(*,0) = sum / float(nazm)

  minval = round(min(polear))

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
	c_linestyle = (lev lt 0.0), levels = lev, /overplot, 		$
	charsize=1.0, max_value=maxval+1.0

  RETURN

END

