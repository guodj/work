;-----------------------------------------------------------------------------------------------------------------
;
;	       Multiple Variables on single Plots
;
;-----------------------------------------------------------------------------------------------------------------

PRO multi_v, single, t, dt, times, moncheck, seconds, same, basetime,     $
		 rows, columns, numofgra, titlenam, unitnam, placement,       $
		 maxminque, numofpts, unconpts, datafile, titleofplot,        $
		 plotwhere, numofrang, rangepl, symtype, percent, opos,	$
		 fontsize, titlesize

  if n_elements(percent) eq 0 then percent = 1.0
  cscale = percent*fontsize
  tscale = percent*titlesize

  opos = [1.0,1.0,0.0,0.0]

  comxlab = strarr(30)
  comxlab = [' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',  $
	     ' ',' ',' ',' ',' ',' ']
  relxlab = strarr(30)
  comxtl  = ' '
  relxtl  = ' '

  btr = 0.0
  etr = 0.0
  nl  = 0

  if same eq 1 then begin

    curvar = 0

    Compute_Axis, times, placement, btr, etr, curvar, relxlab, relxtl, basetime, moncheck, nl, ticktime, nminor

  endif

  sizeoftype = cscale*rows * 1.0

  if rows eq 2 then sizeoftype = cscale*rows*0.75
  if rows eq 1 then sizeoftype = cscale*rows*1.25
  if rows gt 3 then sizeoftype = cscale*2.0

  bottom = 0

  maxar = fltarr(numofrang+1)
  minar = fltarr(numofrang+1)

  if (maxminque eq 1) or (maxminque eq 2) then 				$
    for currang = 0, numofrang do begin

     maxar(currang) = min(single)
     minar(currang) = max(single)

     ranar = where(rangepl eq currang, count)

     for curgra=0, count-1 do begin

       i = 1

       while placement(ranar(curgra),i+1) ne 0 do begin

         if numofpts(placement(ranar(curgra),i)) gt 0 then begin

           nop = numofpts(placement(ranar(curgra),i))-1

           if nop lt 0 then nop = 0
           if max(single(placement(ranar(curgra),i),0:nop)) gt 		$
	     maxar(currang) then 					$
	     maxar(currang) = max(single(placement(ranar(curgra),i),0:nop))
           if min(single(placement(ranar(curgra),i),0:nop)) lt 		$
	     minar(currang) then 					$
	     minar(currang) = min(single(placement(ranar(curgra),i),0:nop))

         endif

         i=i+1

       endwhile

     endfor
	 
  endfor

  if (maxminque eq 3) or (maxminque eq 4) then 			$
    for currang = 0, numofrang do begin

     maxar(currang) = 0.0

     ranar = where(rangepl eq currang, count)

     for curgra=0, count-1 do begin

       i = 1

       while placement(ranar(curgra),i+1) ne 0 do begin

         if numofpts(placement(ranar(curgra),i)) gt 0 then begin

           nop = numofpts(placement(ranar(curgra),i))-1

           if nop lt 0 then nop = 0
           if abs(max(single(placement(ranar(curgra),i),0:nop)) - 	$
	          min(single(placement(ranar(curgra),i),0:nop))) gt 	$
	     maxar(currang) then  					$
             maxar(currang) = 						$
	       abs(max(single(placement(ranar(curgra),i),0:nop)) - 	$
		   min(single(placement(ranar(curgra),i),0:nop)))

         endif

         i=i+1

       endwhile

     endfor
	 
  endfor

  for curgra = 0, numofgra do begin

    if (same ne 1) or (curgra gt numofgra-columns) then bottom = 1

    dc = 3.0/float(rows)

    if plotwhere eq 2 then begin

      curcol = (curgra)/columns
      ymarbeg =  dc * float(curcol+1)
      ymarend = -dc * float(curcol)

    endif else begin

      ymarbeg=0.0
      ymarend=0.0

    endelse

    ytl =  ' (' + strcompress(unitnam(placement(curgra,0)), /remove_all) + ')'

    if same eq 2 then begin

      Compute_Axis, times, placement, btr, etr, curgra, relxlab, 	$
	relxtl, basetime, moncheck, nl, ticktime, nminor

    endif

    if bottom eq 1 then begin

      comxlab = relxlab
      comxtl  = relxtl  

    endif

    if numofpts(placement(curgra,0)) ge 2 then begin

      sin2 = single(placement(curgra,0),0:numofpts(placement(curgra,0))-1)
      t2   = t(placement(curgra,0),0:numofpts(placement(curgra,0))-1)
      dt2  = dt(placement(curgra,0),0:numofpts(placement(curgra,0))-2) 

    endif else begin

      sin2 = [0,0]
      t2   = [0,1]
      dt2  = [0,0]

    endelse

    if maxminque eq 5 then begin

     sinmax = max(sin2)
     sinmin = min(sin2)
     ys = 0

    endif else begin

      if maxminque le 2 then begin

;        sinmax = maxar(rangepl(placement(curgra,0)))
;        sinmin = minar(rangepl(placement(curgra,0)))
	sinmax = maxar(0)
	sinmin = minar(0)
	ys = 0

      endif else begin

        i = 1

	sinmax = max(sin2)
	sinmin = min(sin2)

        while placement(curgra,i+1) ne 0 do begin

	  if numofpts(placement(curgra,i)) gt 0 then begin

	    sin4 = single(placement(curgra,i),0:numofpts(placement(curgra,i))-1)
	    if sinmax lt max(sin4) then sinmax=max(sin4)
	    if sinmin gt min(sin4) then sinmin=min(sin4)

	  endif

	  i=i+1

	endwhile

        drange = (maxar(0) - 	$
	  abs(sinmax - sinmin)) / 2.0
        sinmax = sinmax + drange
        sinmin = sinmin - drange
        ys = 1

      endelse

    endelse

    if sinmin eq sinmax then sinmax = sinmin + 1.0 

    i = 1

    while placement(curgra,i+1) ne 0 do begin

      if numofpts(placement(curgra,i)) gt 0 then begin

        sin4 = single(placement(curgra,i),0:numofpts(placement(curgra,i))-1)
        if sinmax lt max(sin4) then sinmax=max(sin4)
        if sinmin gt min(sin4) then sinmin=min(sin4)

      endif

      i=i+1

    endwhile

    num = i+1

    spt = 0
    
    misdata = where(dt2 GT unconpts, count)

    for i=0,count do begin

      if count eq 0 then ept = numofpts(placement(curgra,0))-1 else begin

	if i gt 0 then spt = misdata(i-1)+1
        if i ne count then ept = misdata(i) else ept = numofpts(placement(curgra,0))-1

      endelse

      if numofpts(placement(curgra,0)) gt 0 then begin

        sin3 = sin2(spt:ept)
        t3   = t2(spt:ept)
    
      endif else begin

        sin3 = sin2
        t3   = t2

      endelse 

      if i gt 0 then oplot, t3, sin3 else begin
 
        if bottom ne 0 then begin
 
          if same eq 1 then plot, t3, sin3, xticks = nl, 		$
		xtickname = relxlab, xtitle= relxtl, ytitle=ytl, 	$
		xstyle=1, xra=[btr,etr], yra=[sinmin,sinmax], 		$
		ymargin=[ymarbeg,ymarend], charsize=sizeoftype, 	$
		xminor=nminor, yminor=5,xtickv=ticktime, ystyle = ys

          if same ne 1 then plot, t3, sin3, xticks = nl, 		$
		xtickname = relxlab, xtitle= relxtl, ytitle=ytl, 	$
		xstyle=1, xra=[btr,etr], yra=[sinmin,sinmax], 		$
		ymargin=[2,2], charsize=sizeoftype, xminor=nminor, 	$
		yminor=5,xtickv=ticktime,ystyle = ys

        endif else begin

          plot, t3, sin3, xticks = nl, xtickname = comxlab, 		$
		ytickname = [' '], xtitle= comxtl, ytitle=ytl, 		$
		xstyle=1, ymargin=[ymarbeg,ymarend], xra=[btr,etr], 	$
		yra=[sinmin,sinmax], charsize=sizeoftype, 		$
		xminor=nminor, yminor=5,xtickv=ticktime,ystyle = ys

        endelse

        ppos, pos
	if pos(0) lt opos(0) then opos(0) = pos(0)
	if pos(1) lt opos(1) then opos(1) = pos(1)
	if pos(2) gt opos(2) then opos(2) = pos(2)
	if pos(3) gt opos(3) then opos(3) = pos(3)

      endelse

    endfor

    if (sinmin lt 0) and (sinmax gt 0) then 				$
	oplot, [btr,etr],[0.0,0.0],linestyle=5

    i      = 1

    xplac  = (etr - btr)/20
    xplac2 = (etr - btr)/8
    xplac3 = xplac2 + (etr - btr)/100

    if (numofpts(placement(curgra,0)) eq 0) and 			$
       (maxminque eq 5) then yplac = 0.8

    ysiz = !y.crange
    yplac = ysiz(0) + (ysiz(1)-ysiz(0))/25.0

    oplot, [xplac,xplac2], [yplac, yplac], linestyle = 0
    xyouts, xplac3, yplac, titlenam(placement(curgra,0)), charsize=cscale

    while placement(curgra,i+1) ne 0 do begin

      if numofpts(placement(curgra,i)) ge 2 then begin

        sin2 = single(placement(curgra,i),0:numofpts(placement(curgra,i))-1)
        t2  = t(placement(curgra,i),0:numofpts(placement(curgra,i))-1)
        dt2 = dt(placement(curgra,i),0:numofpts(placement(curgra,i))-2)

      endif else begin

        sin2 = [0,0]
	t2   = [0,1]
        dt2  = [0,0]

      endelse

      spt = 0
    
      misdata = where(dt2 GT unconpts, count)

      for j=0,count do begin

        if count eq 0 then ept = numofpts(placement(curgra,i))-1 else begin

	  if j gt 0 then spt = misdata(j-1)+1
          if j ne count then ept = misdata(j) else ept = numofpts(placement(curgra,i))-1

        endelse

        if numofpts(placement(curgra,i)) gt 0 then begin

          sin3 = fltarr(numofpts(placement(curgra,i)))
          t3   = fltarr(numofpts(placement(curgra,i)))

          sin3 = sin2(spt:ept)
          t3   = t2(spt:ept)
    
        endif else begin

          sin3 = sin2
          t3   = t2

        endelse 

        oplot, t3, sin3, linestyle = i

      endfor

      xplac  = xplac + (etr - btr)/num
      xplac2 = xplac2 + (etr - btr)/num
      xplac3 = xplac3 + (etr - btr)/num
      
      if (numofpts(placement(curgra,0)) eq 0) and 			$
	 (maxminque eq 5) then yplac = 0.8

      oplot, [xplac,xplac2], [yplac, yplac], linestyle = i
      xyouts, xplac3, yplac, titlenam(placement(curgra,i)), charsize=cscale

      i=i+1

    endwhile

  endfor

  chartime=''

  side  = 'file(s) '+datafile(0)+'   '+datafile(1)
  side1 = datafile(2)+'   '+datafile(3)
  side2 = datafile(4)+'   '+datafile(5)
  xyouts, 0, 0, side, orientation=90, charsize=0.5*cscale, /device
  xyouts, 200, 2000, side1, orientation=90, charsize=0.5*cscale, /device
  xyouts, 400, 2000, side2, orientation=90, charsize=0.5*cscale, /device

  xyouts, 0.5, 1.01, titleofplot, charsize=2.0*tscale, /norm, alignment = 0.5

  RETURN

END

