;------------------------------------------------------------------------------
;
;               Over Lapping Stack Plot (Offset Plots)
;
;------------------------------------------------------------------------------

PRO stack_pl, single, t, dt, times, moncheck, seconds, same, basetime,     $
                 rows, columns, numofvar, titlenam, unitnam, placement,       $
		 maxminque, numofpts, unconpts, offset, datafile, 	      $
	         titleofplot, plotwhere, offdiv, offmul, symtype, 	      $
		 publish, percent, opos, fontsize, titlesize, remmean

  if n_elements(publish) eq 0 then publish = 0
  if n_elements(percent) eq 0 then percent = 1.0
  cscale = percent*fontsize
  tscale = percent*titlesize

  opos = [1.0,1.0,0.0,0.0]

  relxlab = strarr(30)
  relxtl  = ' '

  ytickname = [' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',   $
        ' ',' ',' ',' ',' ',' ',' ',' ']

  btr = 0.0
  etr = 0.0
  nl  = 0

;  if plotwhere eq 1 then ymar = [0,0] else ymar = [3,0]
  ymar = [3,0]

  title = strarr(numofvar+1)

  for i=0,numofvar do 							      $
    title(i) = strcompress(titlenam(placement(0,i)),/remove_all)

  curvar = 0

  Compute_Axis, times, placement, btr, etr, curvar, relxlab, relxtl,     $
		     basetime, moncheck, nl, ticktime, nminor

  if remmean then begin

    s_save = single

    for i=0,numofvar do begin

      ii = placement(0,i)

      if numofpts(ii) gt 0 then begin
	single(ii,0:numofpts(ii)) = 				$
	  single(ii,0:numofpts(ii)) - 			$
	  mean(single(ii,0:numofpts(ii)))
      endif

    endfor

  endif

  sizeoftype = 1.5*cscale

  if offset eq 0.0 then begin

    for i=0,numofvar do begin

      if numofpts(placement(0,i)) gt 0 then begin

        range = abs(max(single(placement(0,i),0:numofpts(placement(0,i))-1))- $
		    min(single(placement(0,i),0:numofpts(placement(0,i))-1)))
        if range gt offset then offset=range

      endif

    endfor

    offset = floor((floor(offset/100.0)*100.0+100.0)*offdiv*offmul)

  endif

  if (numofpts(placement(0,0)) eq 0) then begin

    curvar = 1

    while (numofpts(placement(0,curvar)) eq 0) or (curvar ne numofvar+1) do   $
      curvar=curvar+1

    if curvar eq numofvar+1 then begin

      sin2 = [0,0]
      t2   = [0,1]
      dt2  = [0,0]

    endif else begin

      temp = placement(0,0)
      placement(0,0) = placement(0,curvar)
      placement(0,curvar) = temp

      sin2 = fltarr(numofpts(placement(0,0)))
      t2   = fltarr(numofpts(placement(0,0)))
      dt2  = fltarr(numofpts(placement(0,0)))

    endelse

  endif 

  sinmax = (numofvar+2)*offset
  sinmin = 0.0

  sin2 = sinmax - offset + single(placement(0,0),0:numofpts(placement(0,0))-1)
  t2   = t(placement(0,0),0:numofpts(placement(0,0))-1)
  dt2  = dt(placement(0,0),0:numofpts(placement(0,0))-2) 

  spt = 0

  misdata = where(dt2 GT unconpts, count)

  for i=0,count do begin

    if count eq 0 then ept = numofpts(placement(0,0))-1 else begin

      if i gt 0 then spt = misdata(i-1)+1
      if i ne count then ept = misdata(i) else ept = numofpts(placement(0,0))-1

    endelse

    if numofpts(placement(0,0)) gt 0 then begin

      sin3 = fltarr(numofpts(placement(0,0)))
      t3   = fltarr(numofpts(placement(0,0)))

      sin3 = sin2(spt:ept)
      t3   = t2(spt:ept)

    endif else begin

      sin3 = sin2
      t3   = t2

    endelse 

    if i gt 0 then oplot, t3, sin3, psym=symtype else begin

      plot, t3, sin3, xticks = nl, xtickname = relxlab, xtitle= relxtl,       $
	    ytitle=' ', xstyle=1, xra=[btr,etr], yra=[sinmin,sinmax], 	      $
	    ymargin=ymar, charsize=sizeoftype, yticks = numofvar+2,           $
            xminor=nminor, yminor=5, ytickname=ytickname, ystyle=1,  	      $
	    xtickv=ticktime, psym=symtype

      xyouts, btr-(etr-btr)/100.0, sinmax-offset, title(0), alignment=1.0, $
	    charsize=sizeoftype

      ppos, opos

    endelse

  endfor

  oplot, [btr,etr], [sinmax - offset, sinmax - offset],linestyle=5           

  axis, yaxis=1, ytitle = 'Division Between Base Lines = '+		      $
	strcompress(string(offset),/remove_all)+' '+unitnam(placement(0,0)),  $
        yticks=numofvar+2, ytickname = ytickname, yminor = 5, 		      $
	charsize= sizeoftype

  for curvar=1, numofvar do begin

    if numofpts(placement(0,curvar)) gt 0 then begin

      sin2 = sinmax - offset*(curvar+1) + single(placement(0,curvar),0:numofpts(placement(0,curvar))-1)
      t2   = t(placement(0,curvar),0:numofpts(placement(0,curvar))-1)
      dt2  = dt(placement(0,curvar),0:numofpts(placement(0,curvar))-2) 

    endif else begin

      sin2 = [0,0]
      t2   = [0,1]
      dt2  = [0,0]

    endelse		

    spt = 0
    
    misdata = where(dt2 GT unconpts, count)

    for i=0,count do begin

      if count eq 0 then ept = numofpts(placement(0,curvar))-1 else begin

	if i gt 0 then spt = misdata(i-1)+1
        if i ne count then ept = misdata(i) else ept = numofpts(placement(0,curvar))-1

      endelse

      if numofpts(placement(0,curvar)) gt 0 then begin

        sin3 = sin2(spt:ept)
        t3   = t2(spt:ept)

      endif else begin

	sin3 = sin2
	t3   = t2

      endelse 

      oplot, t3, sin3, psym=symtype

    endfor

    oplot, [btr,etr], 							      $
	   [sinmax-(curvar+1)*offset, sinmax-(curvar+1)*offset],linestyle=5
 
    xyouts, btr-(etr-btr)/100.0, sinmax-(curvar+1)*offset, title(curvar),  $
	    alignment=1.0, charsize=sizeoftype

  endfor

  if publish eq 0 then begin

    side  = 'file(s) '+datafile(0)+'   '+datafile(1)
    side1 = datafile(2)+'   '+datafile(3)
    side2 = datafile(4)+'   '+datafile(5)
    xyouts, -200, 0, side, orientation=90, charsize=0.5*cscale, /device
    xyouts, 0, 2000, side1, orientation=90, charsize=0.5*cscale, /device
    xyouts, 200, 2000, side2, orientation=90, charsize=0.5*cscale, /device

  endif

  xyouts, 0.5, 1.01, titleofplot, charsize=2.0*tscale, /norm, alignment= 0.5

  if remmean then single = s_save

  RETURN

END

