
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;	       One or More Standard X vs. Y Graph(s)
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PRO xy_plt, single, t, dt, times, numofgra, titlenam, unitnam,           $
	         placement, maxminque, numofpts, unconpts,dot, seconds,       $
	         datafile, titleofplot, plotwhere, symtype, maxx, maxy,       $
		 minx, miny, logx, logy, publish, percent, opos,	$
		 fontsize, titlesize, t22 = t22

  if n_elements(maxx) eq 0 then maxx = -1.0e-30
  if n_elements(maxy) eq 0 then maxy = -1.0e-30
  if n_elements(minx) eq 0 then minx = -1.0e-30
  if n_elements(miny) eq 0 then miny = -1.0e-30
  if n_elements(logx) eq 0 then logx = 0
  if n_elements(logy) eq 0 then logy = 0
  if n_elements(percent) eq 0 then percent = 1.0
  if n_elements(t22) eq 0 then t22 = 0
  cscale = percent*fontsize
  tscale = percent*titlesize

  opos = [1.0,1.0,0.0,0.0]

  combin   = fltarr(numofgra+1,2,max(numofpts)+1)
  comtim   = 0.0
  comdt    = fltarr(numofgra+1,max(numofpts)+1)
  compts   = intarr(numofgra+1)

  a = findgen(16) * (!pi*2/16.)
  usersym, 0.5*cos(a), 0.5*sin(a), /fill

  for curgra = 0, numofgra do begin

    dc = 4.0/float(numofgra+1)

    if t22 eq 1 then begin
      ymarbeg = dc*float(curgra+1)
      ymarend = -dc*float(curgra)
    endif else begin
      ymarbeg = 0.0
      ymarend = 0.0
    endelse

    g1 = placement(curgra,0)
    g2 = placement(curgra,1)

    relxtl  = strcompress(seconds(1+times(g1,0,1))+' /'+seconds(1+times(g1,0,2))+' /'+string(times(g1,0,0))+' '+          $
                          seconds(1+times(g1,0,3))+' :'+seconds(1+times(g1,0,4))+' :'+seconds(1+times(g1,0,5))+' to '+       $
                          seconds(1+times(g1,1,1))+' /'+seconds(1+times(g1,1,2))+' /'+string(times(g1,1,0))+' '+          $
                          seconds(1+times(g1,1,3))+' :'+seconds(1+times(g1,1,4))+' :'+seconds(1+times(g1,1,5)))

    xtitle = strcompress(titlenam(g1) + ' ('+unitnam(g1)+')')
    ytitle = strcompress(titlenam(g2) + ' ('+unitnam(g2)+')')

    xtickn = strarr(20)
    ytickn = strarr(20)

    if (t22 eq 1) and (curgra lt numofgra) then begin
      relxtl = ''
      xtitle = ''
      for i=0,19 do xtickn(i) = ' '
      ytickn(0) = ' '
    endif

    if (numofpts(g1) gt 1) and (numofpts(g2) gt 1) then begin

      n1 = 0
      n2 = 0
      curpt = 0

      while (n1 le numofpts(g1)-1) and (n2 le numofpts(g2)-1) do begin

        if round(t(g1, n1)) eq round(t(g2,n2)) then begin

          combin(curgra, 0, curpt) = single(g1,n1)
          combin(curgra, 1, curpt) = single(g2,n2)

          if curpt ne 0 then comdt(curgra, curpt-1) = t(g1,n1) - comtim

          comtim                   = t(g1,n1)

          n1    = n1 + 1
          n2    = n2 + 1
          curpt = curpt + 1

        endif else if round(t(g1, n1)) lt round(t(g2,n2)) then n1 = n1 + 1 else n2 = n2 + 1

      endwhile

      compts(curgra) = curpt - 1

    endif else compts(curgra) = -1

    if dot eq 2 then begin

      if compts(curgra) ne -1 then begin

	if minx eq -1.0e-30 then xmin = min(combin(curgra,0,0:compts(curgra)))    $
	else xmin = minx
	if maxx eq -1.0e-30 then xmax = max(combin(curgra,0,0:compts(curgra)))    $
	else xmax = maxx

	if xmin eq xmax then begin

	  xmin = 0.0
	  xmax = 1.0

	endif 

	if miny eq -1.0e-30 then ymin = min(combin(curgra,1,0:compts(curgra)))    $
	else ymin = miny
	if maxy eq -1.0e-30 then ymax = max(combin(curgra,1,0:compts(curgra)))    $
	else ymax = maxy

	if ymin eq ymax then begin

	  ymin = 0.0
	  ymax = 1.0

	endif 

        if logx eq 1 then begin

	  if xmax le 0 then xmax = 1.0
	  if xmin le 0 then xmin = 10.0^(alog(xmax)-5.0)

	endif

        if logy eq 1 then begin

	  if ymax le 0 then ymax = 1.0
	  if ymin le 0 then ymin = 10.0^(alog(ymax)-5.0)

	endif

	if (logx eq 0) and (logy eq 0) then 				      $
          plot, [combin(curgra,0,0),combin(curgra,0,0)],		      $
                [combin(curgra,1,0),combin(curgra,1,0)],		      $
                xrange = [xmin, xmax], 					      $
                yrange = [ymin, ymax], 					      $
                xtitle = xtitle,    $
                ytitle = ytitle,    $
                xstyle = 1,						      $
                ystyle = 1,						      $
		xtickn = xtickn,	$
		ytickn = ytickn,	$
		ymargin = [ymarbeg, ymarend],	$
                subtitle = relxtl,					      $
		charsize = cscale

	if (logx eq 0) and (logy eq 1) then 				      $
          plot_io, [combin(curgra,0,0),combin(curgra,0,0)],		      $
                [combin(curgra,1,0),combin(curgra,1,0)],		      $
                xrange = [xmin, xmax], 					      $
                yrange = [ymin, ymax], 					      $
                xtitle = xtitle,    $
                ytitle = ytitle,    $
                xstyle = 1,						      $
                ystyle = 1,						      $
		xtickn = xtickn,	$
		ytickn = ytickn,	$
		ymargin = [ymarbeg, ymarend],	$
                subtitle = relxtl,					      $
		charsize = cscale

	if (logx eq 1) and (logy eq 0) then 				      $
          plot_oi, [combin(curgra,0,0),combin(curgra,0,0)],		      $
                [combin(curgra,1,0),combin(curgra,1,0)],		      $
                xrange = [xmin, xmax], 					      $
                yrange = [ymin, ymax], 					      $
                xtitle = xtitle,    $
                ytitle = ytitle,    $
                xstyle = 1,						      $
                ystyle = 1,						      $
		xtickn = xtickn,	$
		ytickn = ytickn,	$
		ymargin = [ymarbeg, ymarend],	$
                subtitle = relxtl,					      $
		charsize = cscale

	if (logx eq 1) and (logy eq 1) then 				      $
          plot_oo, [combin(curgra,0,0),combin(curgra,0,0)],		      $
                [combin(curgra,1,0),combin(curgra,1,0)],		      $
                xrange = [xmin, xmax], 					      $
                yrange = [ymin, ymax], 					      $
                xtitle = xtitle,    $
                ytitle = ytitle,    $
                xstyle = 1,						      $
                ystyle = 1,						      $
		xtickn = xtickn,	$
		ytickn = ytickn,	$
		ymargin = [ymarbeg, ymarend],	$
                subtitle = relxtl,					      $
		charsize = cscale

	  oplot, [combin(curgra,0,0:compts(curgra))],			$
		 [combin(curgra,1,0:compts(curgra))], 	      		$
	  psym = symtype

      endif else begin
 
        plot, [0,1],                                                          $
              [0,0],                                                          $
              yrange = [0,1],                                                 $
              xtitle = xtitle,    $
              ytitle = ytitle,    $
              xstyle = 1,                                                     $
              ystyle = 1,						      $
	      xtickn = xtickn,	$
	      ytickn = ytickn,	$
	      ymargin = [ymarbeg, ymarend],	$
              subtitle = relxtl,					      $
	      charsize = cscale

      endelse

    endif else begin

      if compts(curgra) ne -1 then begin

	if minx eq -1.0e-30 then xmin = min(combin(curgra,0,0:compts(curgra)))    $
	else xmin = minx
	if maxx eq -1.0e-30 then xmax = max(combin(curgra,0,0:compts(curgra)))    $
	else xmax = maxx

	if xmin eq xmax then begin

	  xmin = 0.0
	  xmax = 1.0

	endif 

	if miny eq -1.0e-30 then ymin = min(combin(curgra,1,0:compts(curgra)))    $
	else ymin = miny
	if maxy eq -1.0e-30 then ymax = max(combin(curgra,1,0:compts(curgra)))    $
	else ymax = maxy

	if ymin eq ymax then begin

	  ymin = 0.0
	  ymax = 1.0

	endif 

        if logx eq 1 then begin

	  if xmax le 0 then xmax = 1.0
	  if xmin le 0 then xmin = 10.0^(alog(xmax)-5.0)

	endif

        if logy eq 1 then begin

	  if ymax le 0 then ymax = 1.0
	  if ymin le 0 then ymin = 10.0^(alog(ymax)-5.0)

	endif

	if (logx eq 0) and (logy eq 0) then 				      $
          plot, [combin(curgra,0,0),combin(curgra,0,0)],		      $
                [combin(curgra,1,0),combin(curgra,1,0)],		      $
                xrange = [xmin, xmax], 					      $
                yrange = [ymin, ymax], 					      $
                xtitle = xtitle,    $
                ytitle = ytitle,    $
                xstyle = 1,						      $
                ystyle = 1,						      $
		xtickn = xtickn,	$
		ytickn = ytickn,	$
		ymargin = [ymarbeg, ymarend],	$
                subtitle = relxtl,					      $
		charsize = cscale

	if (logx eq 0) and (logy eq 1) then 				      $
          plot_io, [combin(curgra,0,0),combin(curgra,0,0)],		      $
                [combin(curgra,1,0),combin(curgra,1,0)],		      $
                xrange = [xmin, xmax], 					      $
                yrange = [ymin, ymax], 					      $
                xtitle = xtitle,    $
                ytitle = ytitle,    $
                xstyle = 1,						      $
                ystyle = 1,						      $
		xtickn = xtickn,	$
		ytickn = ytickn,	$
		ymargin = [ymarbeg, ymarend],	$
                subtitle = relxtl,					      $
		charsize = cscale

	if (logx eq 1) and (logy eq 0) then 				      $
          plot_oi, [combin(curgra,0,0),combin(curgra,0,0)],		      $
                [combin(curgra,1,0),combin(curgra,1,0)],		      $
                xrange = [xmin, xmax], 					      $
                yrange = [ymin, ymax], 					      $
                xtitle = xtitle,    $
                ytitle = ytitle,    $
                xstyle = 1,						      $
                ystyle = 1,						      $
		xtickn = xtickn,	$
		ytickn = ytickn,	$
		ymargin = [ymarbeg, ymarend],	$
                subtitle = relxtl,					      $
		charsize = cscale

	if (logx eq 1) and (logy eq 1) then 				      $
          plot_oo, [combin(curgra,0,0),combin(curgra,0,0)],		      $
                [combin(curgra,1,0),combin(curgra,1,0)],		      $
                xrange = [xmin, xmax], 					      $
                yrange = [ymin, ymax], 					      $
                xtitle = xtitle,    $
                ytitle = ytitle,    $
                xstyle = 1,						      $
                ystyle = 1,						      $
		xtickn = xtickn,	$
		ytickn = ytickn,	$
		ymargin = [ymarbeg, ymarend],	$
                subtitle = relxtl,					      $
		charsize = cscale

        spt = 0
    
        misdata = where(comdt(curgra,0:compts(curgra)) GT unconpts, count)

        for i=0,count do begin

          if count eq 0 then ept = compts(curgra) else begin

	    if i gt 0 then spt = misdata(i-1)+1
            if i ne count then ept = misdata(i) else ept = compts(curgra)

          endelse

          sinx = combin(curgra,0,spt:ept)
          siny = combin(curgra,1,spt:ept)

          oplot, sinx, siny

        endfor

	oplot, [combin(curgra,0,0)],[combin(curgra,1,0)],psym=5

      endif else begin
 
        plot, [0,1],                                                          $
              [0,0],                                                          $
              yrange = [0,1],                                                 $
              xtitle = xtitle,    $
              ytitle = ytitle,    $
              xstyle = 1,                                                     $
	      xtickn = xtickn,	$
	      ytickn = ytickn,	$
	      ymargin = [ymarbeg, ymarend],	$
              subtitle = relxtl,					      $
	      charsize = cscale

      endelse

    endelse

    ppos, pos
    if pos(0) lt opos(0) then opos(0) = pos(0)
    if pos(1) lt opos(1) then opos(1) = pos(1)
    if pos(2) gt opos(2) then opos(2) = pos(2)
    if pos(3) gt opos(3) then opos(3) = pos(3)

  endfor

  if publish eq 0 then begin

    side  = 'file(s) '+datafile(0)+'   '+datafile(1)
    side1 = datafile(2)+'   '+datafile(3)
    side2 = datafile(4)+'   '+datafile(5)
    xyouts, 0, 0, side, orientation=90, charsize=0.5*cscale, /device
    xyouts, 200, 2000, side1, orientation=90, charsize=0.5*cscale, /device
    xyouts, 400, 2000, side2, orientation=90, charsize=0.5*cscale, /device

  endif

  xyouts, 0.5, 1.01, titleofplot, charsize=2.0*tscale, /norm, alignment = 0.5

  RETURN

END
