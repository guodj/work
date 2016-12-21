pro plot_data, data, time, nvars, titles, overplots, dt=dt, day = day

  if (n_elements(nVars) eq 0) then begin
     s = size(data)
     print, s
     if (s(0) eq 1) then begin
        nVars = 1
        titles = ['data']
        nPts = s(1)
        data2 = fltarr(1,nPts)
        data2(0,*) = data
        data = data2
     endif
  endif

  if n_elements(overplots) eq 0 then overplots = intarr(n_elements(titles))

  if (n_elements(dt) eq 0) then dt = 1.0e32

  nPts = n_elements(time)
  dtime = time(1:nPts-1) - time(0:nPts-2)
  ltime = where(dtime gt dt, nSkips)
  if (nSkips eq 0) then begin
      ltime = [-1,nPts-1]
  endif else begin
      ltime = [-1,ltime,nPts-1]
  endelse

  stime = min(time)
  c_r_to_a, itime, stime
  if (itime(3) le 1 or n_elements(day) gt 0) then begin
      itime(3) = 0
      itime(4) = 0
      itime(5) = 0
      c_a_to_r, itime, stime
  endif

  etime = max(time)
  c_r_to_a, itime, etime
  if (itime(3) ge 22 or n_elements(day) gt 0) then begin
      itime(3) = 24
      itime(4) = 0
      itime(5) = 0
      c_a_to_r, itime, etime
  endif

  tmptime = time - stime

  time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

  loc = where(overplots eq 0, ppp)
  if ppp eq 0 then ppp = 1
  space = 0.01
  pos_space, ppp, space, sizes, ny = ppp

  plotnum = 0

  for nv = 0,nvars-1 do begin

    if plotnum eq 0 then plotdumb

    get_position, ppp, space, sizes, plotnum, pos, /rect

    pos(0) = pos(0) + 0.1

    if nv lt nvars-1 then begin
      xtn = strarr(60)+' '
      xt  = ' '
    endif else begin
      xtn = xtickname
      xt  = xtitle
    endelse

    loc = where(data(nv,*) gt -1.0e10,count)
    if count gt 0 then range = mm(data(nv,loc)) else range = [-1.0,1.0]

    if (nv lt nvars-1) then begin
        if (overplots(nv+1)) then begin
            loc = where(data(nv+1,*) gt -1.0e10,count)
            if count gt 0 then begin
                range(0) = min([min(data(nv+1,loc)),range(0)])
                range(1) = max([max(data(nv+1,loc)),range(1)])
            endif 
        endif
    endif

    if (range(0) eq range(1)) then begin
      if range(0) eq 0.0 then begin
        range(0) = -1.0
        range(1) = 1.0
      endif else begin
        range(0) = range(0) - 0.1*range(0)
        range(1) = range(1) + 0.1*range(1)
      endelse
    endif else begin
      r = range(1) - range(0)
      if (range(0) gt 0 and range(0)-r*0.1 gt 0) then $
         range(0) = range(0) - r*0.1
      range(1) = range(1) + r*0.1
    endelse

    if (overplots(nv) eq 0) or (nv eq 0) then begin

      plot, tmptime, data(nv,*), xstyle=1,		$
	  ytitle = titles(plotnum),		$
	  xtickname = xtn,			$
	  xtitle = xt,			$
	  xtickv = xtickv,			$
	  xminor = xminor,			$
	  xticks = xtickn,   $
          pos = pos, /noerase, $
          min_val = -1.0e10, yrange = range, $
          thick = 3, /nodata

      for iSkip = 0, nSkips do begin
          oplot, tmptime(ltime(iSkip)+1:ltime(iSkip+1)), $
            data(nv,ltime(iSkip)+1:ltime(iSkip+1)), thick = 3, $
                 min_val = -1.0e10
      endfor

      oplot, [btr,etr],[0.0,0.0],linestyle = 1

      ndays = fix(max(tmptime)/(24.0*3600.0))

      if (ndays lt 10) then $
        for i=1,ndays do oplot, 24.0*3600.0*[float(i),float(i)],   $
        [-max(data(nv,*))*10.0,max(data(nv,*))*10.0], $
        linestyle = 1

      plotnum = plotnum + 1

    endif else begin

        for iSkip = 0, nSkips do begin
            oplot, tmptime(ltime(iSkip)+1:ltime(iSkip+1)), $
              data(nv,ltime(iSkip)+1:ltime(iSkip+1)), $
              linestyle = 2, thick = 3, min_val = -1.0e10
        endfor
        
    endelse

  endfor

  return

end
