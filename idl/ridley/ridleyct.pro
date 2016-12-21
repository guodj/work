pro ridleyct

  window, xsize = 600, ysize = 400, /free

  ncolors = 0
  colors = fltarr(3,256)
  print, 'Enter filename with default color table (return for none):'
  fname = ''
  read,fname
  if strlen(fname) gt 0 then begin
    openr,1,fname
    readf,1,ncolors
    readf,1,colors
    close,1
  endif else begin
    print, 'How many elements would you like on your color table (50-256):'
    read, ncolors
    colors(0,0:ncolors-1) = 255.0*findgen(ncolors)/ncolors
    colors(1,0:ncolors-1) = 255.0*findgen(ncolors)/ncolors
    colors(2,0:ncolors-1) = 255.0*findgen(ncolors)/ncolors
    colors(*,ncolors:255) = 255.0
  endelse

  colors(*,0) = 0.0
  tvlct, colors(0,*),colors(1,*),colors(2,*)

  dy = 0.25
  dx = 0.7
  sx = 0.07
  sy = 0.07
  space = 0.07

  posr = [sx,sy+2.0*(space+dy),sx+dx,sy+2.0*space+3.0*dy]
  posg = [sx,sy+1.0*(space+dy),sx+dx,sy+1.0*space+2.0*dy]
  posb = [sx,sy,sx+dx,sy+dy]

  pos = fltarr(3,4)
  pos(0,*) = posr
  pos(1,*) = posg
  pos(2,*) = posb

  posct = [sx+dx+space,0.1,sx+dx+2.0*space,0.9]

  cpts = 5
  cloc = (ncolors-1)*(findgen(cpts+2))/float(cpts+1)
  cloc(cpts+1) = ncolors-1
  clip = fltarr(3,4)

  rgb = ['red','green','blue']

  plot, colors(0,0:ncolors-1), pos = posr, xstyle=1, ystyle = 1, 	$
	yrange = [0.0,255.0], ytitle = rgb(0)
  clipr = !p.clip
  clip(0,0:3) = clipr(0:3)
  oplot, cloc, colors(0,cloc), color = 255, psym = 4

  plot, colors(1,0:ncolors-1), pos = posg, xstyle=1, ystyle = 1, 	$
	yrange = [0.0,255.0], /noerase, ytitle = rgb(1)
  clipg = !p.clip
  clip(1,0:3) = clipg(0:3)
  oplot, cloc, colors(1,cloc), color = 255, psym = 4

  plot, colors(2,0:ncolors-1), pos = posb, xstyle=1, ystyle = 1, 	$
	yrange = [0.0,255.0], /noerase, ytitle = rgb(2)
  clipb = !p.clip
  clip(2,0:3) = clipb(0:3)
  oplot, cloc, colors(2,cloc), color = 255, psym = 4

  x = [0.0,0.0,1.0,1.0,0.0]
  y = [0.0,1.0,1.0,0.0,0.0]

  plot, [0.0,1.0],[0.0,ncolors-1], /nodata, /noerase, 		$
    xstyle=5, ystyle=1, pos=posct

  for i=0,ncolors-1 do polyfill,x,y+float(i),color=i

  plots, [0.0,0.0],[0.0,ncolors-1]
  plots, [1.0,1.0],[0.0,ncolors-1]
  plots, [0.0,1.0],[0.0,0.0]
  plots, [0.0,1.0],[ncolors-1.0,ncolors-1.0]

  done = 0

  while not done do begin

    cursor,x,y,1,/device
    mouse = !err

    while mouse eq 1 do begin

      a = -1
      for c = 0, 2 do begin
        if (x ge clip(c,0) and x le clip(c,2)) and 			$
	   (y ge clip(c,1) and y le clip(c,3)) then begin
          i = (ncolors-1.0)*float(x - clip(c,0))/float(clip(c,2)-clip(c,0))
          j = 255.0*float(y - clip(c,1))/float(clip(c,3)-clip(c,1))
          a = c
        endif
      endfor

      if (a ge 0) then begin
        plot, colors(a,0:ncolors-1), pos = pos(a,*), 			$
	  xstyle=1, ystyle = 1, color = 0,				$
	  yrange = [0.0,255.0], /noerase
	oplot, cloc, colors(a,cloc), color = 0, psym = 4

	dis = abs(cloc-i)
        loc = where(dis eq min(dis))
        loc = loc(0)

        if loc lt cpts+1 then begin
	  dy = float(colors(a,cloc(loc+1)) - j)/(cloc(loc+1)-cloc(loc))
	  for ii=cloc(loc),cloc(loc+1)-1 do 				$
	    colors(a,ii) = j+dy*float(ii-cloc(loc))
        endif

        if loc gt 0 then begin
	  dy = float(j - colors(a,cloc(loc-1)))/(cloc(loc)-cloc(loc-1))
	  for ii=cloc(loc-1),cloc(loc) do 				$
	    colors(a,ii) = colors(a,cloc(loc-1))+dy*float(ii-cloc(loc-1))
        endif

	colors(a,cloc(loc)) = j

        plot, colors(a,0:ncolors-1), pos = pos(a,*), 			$
	  xstyle=1, ystyle = 1, color = 255,				$
	  yrange = [0.0,255.0], /noerase
	oplot, cloc, colors(a,cloc), color = 255, psym = 4

        colors(*,0) = 0.0
        tvlct, colors(0,*),colors(1,*),colors(2,*)

      endif

      cursor,x,y,1,/device
      mouse = !err

    endwhile

    if mouse eq 2 then begin

      cpts = cpts + 1

      for a=0,2 do begin
        plot, colors(a,0:ncolors-1), pos = pos(a,*), 			$
	  xstyle=1, ystyle = 1, color = 0,				$
	  yrange = [0.0,255.0], /noerase
        oplot, cloc, colors(a,cloc), color = 0, psym = 4
      endfor

      cloc = (ncolors-1.0)*(findgen(cpts+2))/float(cpts+1)
      cloc(cpts+1) = ncolors-1.0

      for a=0,2 do begin
        plot, colors(a,0:ncolors-1), pos = pos(a,*), 			$
	  xstyle=1, ystyle = 1, color = 255,				$
	  yrange = [0.0,255.0], /noerase
        oplot, cloc, colors(a,cloc), color = 255, psym = 4
      endfor

      while mouse eq 2 do begin
        cursor,x,y,0,/device
        mouse = !err
      endwhile

    endif

    if mouse eq 4 then done = 1

  endwhile

  que=''
  print, 'Would you like this color table to be saved in a file (y/n) ?'
  read,que
  if (strmid(que,0,1) eq 'y') then begin
    fname = ''
    print, 'Enter filename :'
    read,fname
    openw,1,fname
    colors(*,0) = 0.0
    printf,1,ncolors
    printf,1,colors
    close,1
  endif

  wdelete

  return

end