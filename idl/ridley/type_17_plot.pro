;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;              Density Plots
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO Type_17_Plot, single, t, times, moncheck, seconds, basetime, titlenam,         $
                  unitnam, numofpts, maxran, plotwhere, datafile, titleofplot, minden, maxden
    
  red     = intarr(201)
  blue    = intarr(201)
  green   = intarr(201)  
  relxlab = strarr(30)

  curvar    = 0
  altitude  = 0
  azimuth   = 1
  density   = 2
  elevation = 3
  event     = 4

  sizeoftype = 1.5

  if (numofpts(altitude) gt 0) then minalt = min(single(altitude,0:numofpts(altitude)-1)) else minalt = 0.0
  if (minalt gt 25.0) then minalt = minalt - 25.0
  maxalt = maxran

  ytlr =  ' (' + strcompress(unitnam(altitude), /remove_all) + ')'

  ytl  = strcompress(titlenam(altitude), /remove_all) + ' ' + ytlr

  tickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']

  numofeve = 0

  if numofpts(event) gt 1 then begin

    for i = 1, numofpts(event)-1 do if single(event,i) ne single(event,i-1) then numofeve = numofeve + 1

  endif

  startpt = 0
  dummy   = 0
  starttm = double(0.0)

;  minden = 1e4
;  if (numofpts(density) gt 0) then maxden = max(single(density,0:numofpts(density)-1)) else maxden = 1.0e10
;  maxden = 0.8 * maxden
;  if (numofpts(density) gt 0) then minden = min(single(density,0:numofpts(density)-1)) else minden = 0.0
;  minden = 1.2 * minden
;  maxden = 1.0e6
;  maxden = 1.2e12
;  minden = 2.0e10
  dden = 200.0 / alog(maxden/minden)

  for i=1,99 do begin

    red(i) = 0.0
    green(i) = 255.0*(float(i)-1.0)/99.0
    blue(i) = 255.0 - 255.0*(float(i)-1.0)/99.0

  endfor

  for i=100,199 do begin

    red(i) =  255.0*(float(i)-100.0)/99.0
    green(i) = 255.0 - 255.0*(float(i)-100.0)/99.0
    blue(i) = 0.0

  endfor

  red(200) = 255
  green(200) = 255
  blue(200) = 255

  tvlct, red, green, blue

  for cureve = 0, numofeve do begin

    starttm = double(basetime(event)) + double(t(event,startpt))

    c_r_to_a, timearray, starttm

    datetl  = strcompress(strmid(moncheck,(timearray(1)-1)*3,3)+' '+string(timearray(2))+ $
              ', 19'+strcompress(string(timearray(0)),/remove_all))

    inserttl = strcompress(seconds(1+timearray(3)) +                                         $
                           ' :'+seconds(1+timearray(4)) + ' :'+seconds(1+timearray(5))+ ' UT')

    search = startpt+1
    found = 0

    while (search lt numofpts(elevation)-1) and (found eq 0) do begin

      if (round(single(elevation,startpt)) ne round(single(elevation,search))) then begin

        change = elevation 
        insertran = 'Range ='+strcompress(string(round(maxalt)))+' Km'
        insertsca = 'Elevation Scan'
        insertang = strcompress('Azimuth ='+string(round(single(azimuth,startpt)))+'!Eo')
	found = 1

      endif else begin

        if (round(single(azimuth,startpt)) ne round(single(azimuth,search))) then begin

          insertran = 'Range ='+strcompress(string(round(maxalt)))+ 'Km'
          insertsca = 'Azimuth Scan'
          insertang = strcompress('Elevation ='+string(round(single(elevation,startpt)))+'!Eo')
          change = azimuth
	  found = 1

	endif else search = search + 1

      endelse

    endwhile


    dr = maxalt/20.0
    dr2 = 2.0*750.0/maxalt

    if (change eq azimuth) then begin

      ymax = maxalt
      ymin = -maxalt
      posx1 = 0.05
      posy1 = 0.0
      posx2 = posx1 + 0.95*(7.5/10.0)
      posy2 = posy1 + 0.95
      plot, [-maxalt], [ymin],						$
            charsize  = sizeoftype,					$
            xra       = [-maxalt, maxalt],				$
            yra       = [ymin,ymax],					$
            ytitle    = ytl,						$
            xstyle    = 5,						$
            ystyle    = 5,						$
	    xtickname = tickname,					$
	    ytickname = tickname,					$
            pos=[posx1, posy1, posx2, posy2]

    endif else begin

      ymax = maxalt
      ymin = 0.0
      posx1 = 0.05
      posy1 = 0.35
      posx2 = posx1 + (7.5/10.0)
      posy2 = posy1 + 0.5
      plot, [-maxalt], [ymin],						$
            charsize  = sizeoftype,					$
            xra       = [-maxalt, maxalt],				$
            yra       = [ymin,ymax],					$
            ytitle    = ytl,						$
            xstyle    = 5,						$
            ystyle    = 5,						$
	    xtickname = tickname,					$
	    ytickname = tickname,					$
            pos=[posx1, posy1, posx2, posy2]

    endelse

    if numofpts(density)-startpt gt 0 then begin
  
      curpt  = startpt
      eveold = single(event,startpt)
      rold = minalt

      while (curpt le numofpts(density) - 1) and (eveold eq single(event,curpt)) do begin
                                                 
        colval = 0.0

        if (single(density,curpt) gt 0.0) then colval = alog(single(density,curpt)/minden)*dden else colval = 0

	if colval gt 200 then colval = 200
	if colval lt 0 then colval = 0

        r2   = single(altitude,curpt)
	if rold gt r2 then r1 = minalt else r1   = rold

	if change eq elevation then azm1 = single(change,curpt) - 90.0 else $
	  if single(elevation,curpt) gt 90.0 then azm1 = single(change,curpt) + 180.0 else $
	    azm1 = single(change,curpt)

	if azm1 lt 0.0 then azm1 = azm1 + 360.0
	if azm1 gt 360.0 then azm1 = azm1 - 360.0

	found = 0
	search = curpt+1

	while found eq 0 do begin

	  if numofpts(change) - 1 lt search then begin

	    azm2 = azm1 + 10.0

	    found = 1

	  endif else begin

	    if (eveold ne single(event,search)) then begin

	      if da lt 0.0 then begin

		azm2 = azm1 - 10.0
	        found = 1 

	      endif else begin

	       azm2 = azm1 + 10.0
	       found = 1	 

	      endelse

	    endif else begin

              if single(change,search) eq single(change,curpt) then search = search + 1 else begin

	        found = 1

	        if change eq elevation then azm2 = single(change,search) - 90.0 else $
	          if single(elevation,curpt) gt 90.0 then azm2 = single(change,search) + 180.0 else  $
	            azm2 = single(change,search)

	        if azm2 lt 0.0 then azm2 = azm2 + 360.0
	        if azm2 gt 360.0 then azm2 = azm2 - 360.0

	      endelse

	    endelse

	  endelse

	endwhile

	da = azm2 - azm1

	if abs(da) gt 300.0 then if azm2 gt azm1 then da = da - 360.0 else da = da + 360.0

	da = da / 10.0

        a  = (azm1 + (da * findgen(11))) * (!pi*2.0/360.0)
        a2 = (azm1 + (da * (10-findgen(11)))) * (!pi*2.0/360.0)

        if (change eq elevation) and (single(azimuth,curpt) gt 63.0) and (single(azimuth,curpt) lt 243.0) then begin
	  a = 2.0*!pi - a
	  a2 = 2.0*!pi - a2
	endif

        if (single(density,curpt) gt 0.0) and (r2 lt maxalt) then      $
	   polyfill, [r2*sin(a),r1*sin(a2)],[r2*cos(a),r1*cos(a2)], color = colval

	rold = r2

        eveold = single(event,curpt)
        curpt  = curpt + 1

      endwhile

    endif

    startpt = curpt

      posy1 = 0.05
      posy2 = 0.95

      dposy=(posy2-posy1)/4.0
      py1 = posy1+dposy
      py2 = posy2-dposy

      x=findgen(361)*2.0*!pi/360.0

      for i=maxalt/3.0,maxalt, maxalt/3. do oplot, float(i)*sin(x), float(i)*cos(x)
      oplot, [-maxalt,maxalt],[0,0]
      oplot, [0,0],[ymin,ymax]

      values = fltarr(4)

      values(0) = minden
      values(1) = (maxden - minden) * 0.25
      values(2) = (maxden - minden) * 0.75
      values(3) = maxden

      plot, [dr,dr],[minden,minden],xrange=[0,dr],yrange=[minden,maxden],xstyle=5,ystyle=5,     $
            position=[posx2+0.03,py1,posx2+0.08,py2],/noerase, /ytype, yticks=3, ytickv = values

      axis, yaxis=1, ytitle = 'Density in '+unitnam(density),xrange=[0,dr],yrange=[minden,maxden],xstyle=5,ystyle=5,  $
	    yticks = 3,/ytype, ytickv = values

      step = (maxden - minden) / 400.0	
      for curcol = minden,maxden,step do begin
         polyfill,[0,0,dr,dr],[curcol,curcol+step,curcol+step,curcol],color=alog(curcol/minden)*dden
      endfor

      xyouts, posx2+0.03, py1-0.05, insertsca,charsize=1.5,/norm	
      xyouts, posx2+0.03, py1-0.10, insertran,charsize=1.5,/norm	
      xyouts, posx2+0.03, py1-0.15, insertang,charsize=1.5,/norm	
      xyouts, posx2+0.03, py2+0.1, inserttl, charsize=1.5, /norm
      xyouts, posx2+0.03, py2+0.15, datetl, charsize=1.5, /norm
    
      chartime=''
      spawn,'show time',chartime

      side  = 'date printed : '+chartime+' from file(s) '+datafile(0)+'   '+datafile(1)
      side1 = datafile(2)+'   '+datafile(3)
      side2 = datafile(4)+'   '+datafile(5)
      xyouts, 0, 0, side, orientation=90, charsize=0.5, /norm
      xyouts, 0.02, 0.2, side1, orientation=90, charsize=0.5, /norm
      xyouts, 0.04, 0.2, side2, orientation=90, charsize=0.5, /norm

      xyouts, 0.5, 1.01, titleofplot, charsize=2.0, /norm, alignment=0.5

      if (plotwhere ne 1) and (numofeve ne cureve) then prompt_for_next

  endfor

  !p.multi = 0

  RETURN

END

