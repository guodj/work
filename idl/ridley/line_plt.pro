

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;              One or More Standard X vs. T Graph(s)
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO line_plt, single, t, dt, times, moncheck, seconds, same, basetime,     $
		 rows, columns, numofvar, titlenam, unitnam, placement,       $
		 maxminque, numofpts, unconpts, rangepl, numofrang,	      $
                 datafile, titleofplot, plotwhere, symtype, publish,	      $
		 percent, opos, fontsize, titlesize, rangescale

  if n_elements(publish) le 0 then publish = 0
  if n_elements(percent) le 0 then percent = 1.0

  cscale = fontsize*percent
  tscale = titlesize*percent
  ylabs = strarr(20)
  if fontsize gt 1.0 then for i=0,19,2 do ylabs(i) = ' '

  opos = [1.0,1.0,0.0,0.0]

  comxlab = strarr(30)+' '
  relxlab = strarr(30)+' '
  comxtl  = ' '
  relxtl  = ' '

  btr = 0.0
  etr = 0.0
  nl  = 0

  if same eq 1 then begin

    curvar = 0

    Compute_Axis, times, placement, btr, etr, curvar, relxlab, relxtl,   $
		       basetime, moncheck, nl, ticktime, nminor

  endif

  sizeoftype = rows * 1.0

  if rows eq 2 then sizeoftype = rows*0.75
  if rows eq 1 then sizeoftype = rows*1.25
  if rows gt 3 then sizeoftype = 2.0

  bottom = 0

  maxar = fltarr(numofrang+1)
  minar = fltarr(numofrang+1)

  if (maxminque eq 1) or (maxminque eq 2) then 				      $
    for currang = 0, numofrang do begin

     maxar(currang) = min(single)
     minar(currang) = max(single)

     ranar = where(rangepl eq currang, count)

     for i=0, count-1 do begin

       nop = numofpts(placement(ranar(i),0))-1

       if nop lt 0 then nop = 0
       if max(single(ranar(i),0:nop)) gt maxar(currang) then     $
	 maxar(currang) = max(single(ranar(i),0:nop))
       if min(single(ranar(i),0:nop)) lt minar(currang) then     $
	 minar(currang) = min(single(ranar(i),0:nop))

     endfor

     if strlen(rangescale(currang,0)) gt 0 then 			$
	minar(currang) = float(rangescale(currang,0))
     if strlen(rangescale(currang,1)) gt 0 then 			$
	maxar(currang) = float(rangescale(currang,1))
 
  endfor

  if (maxminque eq 3) or (maxminque eq 4) then 				      $
    for currang = 0, numofrang do begin

     maxar(currang) = 0.0

     ranar = where(rangepl eq currang, count)

     for i=0, count-1 do begin

       nop = numofpts(placement(ranar(i),0))-1

       if nop lt 0 then nop = 0
       if abs(max(single(ranar(i),0:nop)) - 		      $
          min(single(ranar(i),0:nop))) gt maxar(currang) then    $
          maxar(currang) = 1.2*abs(max(single(ranar(i),0:nop)) -     $
			   min(single(ranar(i),0:nop)))

     endfor
	 
     if strlen(rangescale(currang,0)) gt 0 then 			$
	maxar(currang) = float(rangescale(currang,0))

  endfor

  for curvar=0, numofvar do begin

    if (same ne 1) or (curvar gt numofvar-columns) then bottom = 1

    dc = 3.0/float(rows)

    if plotwhere eq 2 then begin

      curcol = (curvar)/columns
      ymarbeg=dc * float(curcol+1)
      ymarend=-dc* float(curcol)

    endif else begin

      ymarbeg=0.0
      ymarend=0.0

    endelse

    ytlr =  ' (' + strcompress(unitnam(placement(curvar,0)), /remove_all) + ')'

    ytl  = strcompress(titlenam(placement(curvar,0)), /remove_all)

    if same eq 2 then begin

      Compute_Axis, times, placement, btr, etr, curvar, relxlab, 	      $
		         relxtl, basetime, moncheck, nl, ticktime, nminor

    endif

    if bottom eq 1 then begin

      comxlab = relxlab
      comxtl  = relxtl  

    endif

    if numofpts(placement(curvar,0)) gt 0 then begin

      sin2 = single(placement(curvar,0),0:numofpts(placement(curvar,0))-1)
      t2   = t(placement(curvar,0),0:numofpts(placement(curvar,0))-1)
      dt2  = dt(placement(curvar,0),0:numofpts(placement(curvar,0))-2) 

    endif else begin

      sin2 = [0,0]
      t2   = [0,1]
      dt2  = [0,0]

    endelse		

    if maxminque eq 5 then begin

     sinmax = max(sin2)
     sinmin = min(sin2)
     if strlen(rangescale(placement(curvar,0),0)) gt 0 then 		$
	sinmin = float(rangescale(placement(curvar,0),0))
     if strlen(rangescale(placement(curvar,0),1)) gt 0 then 		$
	sinmax = float(rangescale(placement(curvar,0),1))
     ys = 0

    endif else begin

      if maxminque le 2 then begin

        sinmax = maxar(rangepl(placement(curvar,0)))
        sinmin = minar(rangepl(placement(curvar,0)))
	ys = 0

      endif else begin

        drange = (maxar(rangepl(placement(curvar,0))) - abs(max(sin2) -  $
		  min(sin2))) / 2.0
        sinmax = max(sin2) + drange
        sinmin = min(sin2) - drange
        ys = 1

      endelse

    endelse

    spt = 0

    misdata = where(dt2 GT unconpts, count)

    for i=0,count do begin

      if count eq 0 then ept = numofpts(placement(curvar,0))-1 else begin

	if i gt 0 then spt = misdata(i-1)+1
        if i ne count then 						      $
	  ept = misdata(i) 						      $
	else								      $
	  ept = numofpts(placement(curvar,0))-1

      endelse

      if numofpts(placement(curvar,0)) gt 0 then begin

        sin3 = sin2(spt:ept)
        t3   = t2(spt:ept)
    
      endif else begin

        sin3 = sin2
        t3   = t2

      endelse

      if sinmin eq sinmax then sinmax = sinmin + 1.0 

      if i gt 0 then oplot, t3, sin3, psym=symtype else begin
 
        if bottom ne 0 then begin
 
          if same eq 1 then plot, t3, sin3, xticks = nl, xtickname = relxlab, $
                xtitle= relxtl, ytitle=ytl, xstyle=1, xra=[btr,etr],          $
	        yra=[sinmin,sinmax], ymargin=[ymarbeg,ymarend],		      $
                charsize=cscale*sizeoftype, xminor=nminor, yminor=5, 	      $
		ystyle = ys, ytickname = ylabs,   $
		xtickv=ticktime, psym = symtype 

          if same ne 1 then plot, t3, sin3, xticks = nl, xtickname = relxlab, $
                xtitle= relxtl, ytitle=ytl, xstyle=1, xra=[btr,etr], 	      $
		yra=[sinmin,sinmax], ymargin=[2,2], 			      $
		charsize=cscale*sizeoftype, ytickname = ylabs,     $
		xminor=nminor, yminor=5, ystyle =ys, xtickv=ticktime,	      $
		psym= symtype

        endif else begin

          if ys eq 1 then						      $ 
		plot, t3, sin3, xticks = nl, xtickname = comxlab,	      $
		xtitle= comxtl, ytitle=ytl, xstyle=1, 			      $
		ymargin=[ymarbeg,ymarend], xra=[btr,etr], yra=[sinmin,sinmax],$
		charsize=cscale*sizeoftype, xminor=nminor, yminor=5, 	      $
		ystyle =ys, ytickname = ylabs,    $
		xtickv=ticktime, psym=symtype       			      $
	  else								$
            plot, t3, sin3, xticks = nl, xtickname = comxlab, 		      $
		ytickname=ylabs, xtitle= comxtl, ytitle=ytl, xstyle=1,	      $
		ymargin=[ymarbeg,ymarend], xra=[btr,etr], yra=[sinmin,sinmax],$
                charsize=cscale*sizeoftype, xminor=nminor, yminor=5, 	      $
		ystyle =ys,    $
		xtickv=ticktime, psym=symtype
        endelse

	ppos, pos
	if pos(0) lt opos(0) then opos(0) = pos(0)
	if pos(1) lt opos(1) then opos(1) = pos(1)
	if pos(2) gt opos(2) then opos(2) = pos(2)
	if pos(3) gt opos(3) then opos(3) = pos(3)

      endelse

    endfor

    axis, yaxis=1, ytitle = ytlr, yra=[sinmin,sinmax], yminor=5,	      $
          charsize= cscale*sizeoftype, 				      $
	  ytickname =[' ',' ',' ',' ',' ',' ',' ',' ',  $
	  ' ',' '], ystyle = ys 

    if (sinmin lt 0) and (sinmax gt 0) then 				      $
      oplot, [btr,etr],[0.0,0.0],linestyle=5

  endfor

  if publish eq 0 then begin

    chartime=''
    side  = 'file(s) '+datafile(0)+'   '+  $
	    datafile(1)
    side1 = datafile(2)+'   '+datafile(3)
    side2 = datafile(4)+'   '+datafile(5)
    xyouts, 0, 0, side, orientation=90, charsize=0.5*cscale, /device
    xyouts, 200, 2000, side1, orientation=90, charsize=0.5*cscale, /device
    xyouts, 400, 2000, side2, orientation=90, charsize=0.5*cscale, /device

  endif

  xyouts, 0.5, 1.01, titleofplot, charsize=2.0*tscale, /norm, alignment = 0.5

  RETURN

END

