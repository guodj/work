;----------------------------------------------------------------------------
;
;		Magnetometer Quick Views
;
;----------------------------------------------------------------------------

PRO mag_qick, single, t, dt, times, moncheck, seconds, same, basetime,    $
	          numofvar, titlenam, unitnam, placement, maxminque, numofpts,$
		  unconpts, offset, datafile, titleofplot, stat, firstnum,    $
		  plotwhere, offdiv, offmul, symtype, remmean, percent, opos, $
		  fontsize, titlesize
  if n_elements(percent) eq 0 then percent = 1.0
  cscale = percent*fontsize
  tscale = percent*titlesize

  opos = [1.0,1.0,0.0,0.0]

  btr = 0.0
  etr = 0.0
  nl  = 0.0

  order = [1,0,2]

  !P.Multi=[0,3,1]

  mf = 1.0/3.0

  comp = ['H          H','E          E','Z          Z']

  offset = 0.0

  max_val = 1000.0
  missing = -1.0e32

  if remmean eq 1 then begin

    s_save = single

    for i=0,numofvar do begin

      if numofpts(i) gt 0 then begin
        sin2 = single(placement(0,i),0:numofpts(i))	
        mn = mean(sin2)
        sin2 = sin2 - mn
        single(placement(0,i),0:numofpts(i)) = sin2
      endif

    endfor

  endif

  for i=0,numofvar-1 do begin

    if numofpts(placement(0,i)) gt 0 then begin

      range = abs(max(single(placement(0,i),0:numofpts(i)-1))-min(single(placement(0,i),0:numofpts(i)-1)))
      if range gt offset then offset=range

    endif

  endfor

  offset = float(floor(offset/100.0)*100.0+100.0)*offdiv*offmul

  print, 'Offset = ',offset

  nv = numofvar/3

  if (nv gt firstnum) and (nv gt 10) then numplot = 1 else begin

    numplot = 0
    firstnum=nv+1

  endelse

  for twice = 0,numplot do begin

    if twice eq 0 then begin
      first = 0
      last = firstnum-1
    endif else begin
      first = firstnum
      last = nv
    endelse

    swtitle = strarr(last-first+1)

    for i = 0, last-first do 					$
      swtitle(i) = strcompress(titlenam(last-i),/remove_all)

    for gord = 0,2 do begin

      graph = order(gord)

      curvar = first*3 + graph

      relxlab=strarr(30)
      relxtl = ' '

      Compute_Axis, times, placement, btr, etr, curvar, relxlab, relxtl,      $
		         basetime, moncheck, nl, ticktime, nminor

      sizeoftype = 1.5*cscale

      if (numofpts(curvar) le 0) then begin

        sin2 = [0,0]
        t2   = [0,1]
        dt2  = [0,0]

      endif else begin

        sin2 = single(curvar,0:numofpts(curvar)-1)
        t2   = t(curvar,0:numofpts(curvar)-1)
        dt2  = dt(curvar,0:numofpts(curvar)-2) 

      endelse

      sinmax = ((last-first)+2)*offset
      sinmin = 0.0
      sin2 = sinmax - offset + sin2

      spt = 0

      misdata = where(dt2 GT unconpts, count)

      for i=0,count do begin

        if count eq 0 then ept = numofpts(curvar)-1 else begin

          if i gt 0 then spt = misdata(i-1)+1
          if i ne count then ept = misdata(i) else ept = numofpts(curvar)-1

        endelse

        if numofpts(curvar) gt 0 then begin

          sin3 = sin2(spt:ept)
          t3   = t2(spt:ept)

        endif else begin

          sin3 = sin2
          t3   = t2

        endelse 

        if i gt 0 then oplot, t3, sin3, psym = symtype else begin
 
  	  if gord ne 0 then begin

	    for j = 0, (last-first) do swtitle(j) = ' '

	  endif

          plot, t3, sin3, 						      $
	        xticks = nl, 						      $
	        xtickname = relxlab, 					      $
	        xtitle= relxtl, 					      $
	        ytitle=' ', 						      $
	        xstyle=1,                          			      $
                xra=[btr,etr], 						      $
                xmargin=[2.0-float(gord)*mf,-(2.0/3.0)+float(gord)*mf],       $
	        yra=[sinmin,sinmax], 					      $
	        ymargin=[1.0,0.0], 					      $
	        charsize=sizeoftype,					      $
	        yticks = last-first+2,         	       			      $
                xminor=nminor,						      $
	        yminor=5,						      $
	        ytickname=[' ',swtitle,' '], 				      $
	        ystyle=1,						      $
	        xtickv=ticktime,					      $
	        title = comp(gord),					      $
		psym = symtype

	  ppos, pos
	  if pos(0) lt opos(0) then opos(0) = pos(0)
	  if pos(1) lt opos(1) then opos(1) = pos(1)
	  if pos(2) gt opos(2) then opos(2) = pos(2)
	  if pos(3) gt opos(3) then opos(3) = pos(3)

        endelse

      endfor

      oplot, [btr,etr], [sinmax - offset, sinmax - offset],linestyle=5                                                

      if gord eq 2 then begin

        axis, yaxis=1, 												$
              ytitle = 'Division Between Base Lines = '+strcompress(string(offset),/remove_all)+ ' nT',		$
              yticks=last-first+2, 										$
	      ytickname = [' ',swtitle,' '], 									$
	      yminor = 5, 											$
	      charsize= sizeoftype

      endif

      for curvar=first+1, last do begin

        varnum = curvar*3+graph

        if numofpts(varnum) gt 2 then begin

          sin2 = single(varnum,0:numofpts(varnum)-1)
          t2   = t(varnum,0:numofpts(varnum)-1)
          dt2  = dt(varnum,0:numofpts(varnum)-2)

	  sin2 = sinmax - offset*(curvar-first+1) + sin2

        endif else begin

          sin2 = [0,0]
          t2   = [0,1]
          dt2  = [0,0]

        endelse		

        spt = 0
    
        misdata = where(dt2 GT unconpts, count)

        for i=0,count do begin

          if count eq 0 then ept = numofpts(varnum)-1 else begin

	    if i gt 0 then spt = misdata(i-1)+1
            if i ne count then ept = misdata(i) else ept = numofpts(varnum)-1

          endelse

          if numofpts(varnum) gt 0 then begin

            sin3 = sin2(spt:ept)
            t3   = t2(spt:ept)

          endif else begin

	    sin3 = sin2
	    t3   = t2

          endelse 
          oplot, t3, sin3, psym=symtype

        endfor

        oplot,[btr,etr], 						      $
	      [sinmax-(curvar-first+1)*offset,sinmax-(curvar-first+1)*offset],$
	      linestyle=5 

      endfor

    endfor

    chartime=''
    spawn,'show time',chartime

    side  = 'date printed : '+chartime+' from file(s) '+datafile(0)+'   '+    $
	    datafile(1)
    xyouts, 0.0, -0.05, side, charsize=0.5*cscale, /norm

    xyouts, 0.5, 1.01, titleofplot, charsize=2.0*tscale, /norm, 	      $
	alignment= 0.5

    if (plotwhere ne 1) and (twice ne numplot) then prompt_for_next

  endfor

  if remmean eq 1 then begin

    single = s_save

  endif

  RETURN

END

