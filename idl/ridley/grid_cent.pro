pro grid_cent, nshrooms, grid, cent, clen, n_cent

	common objects, ob

	xsize	= !d.x_size/20
	ysize	= !d.y_size/20

	xpos	= fix(xsize*randomu(s,nshrooms))
	ypos	= fix(((ysize-1)*randomu(s,nshrooms))*.7+.33*ysize)
	here	= where(xpos ge xsize-1)
	if here(0) ne -1 then xpos(here) = xsize-1
	here	= where(ypos ge ysize-2)
	if here(0) ne -1 then ypos(here) = ysize-2

	loc	= sort(randomu(s,nshrooms))
	loc1	= loc(0:nshrooms/2-1)
	loc2	= loc(nshrooms/2:nshrooms-1)

	grid 	= intarr(xsize,ysize)
	grid(xpos(loc1),ypos(loc1)) = 1
	grid(xpos(loc2),ypos(loc2)) = 2

	for i=0,n_elements(loc1)-1 do					$
	  tv, ob.mushroom1, 20*xpos(loc1(i)),20*ypos(loc1(i))
	for i=0,n_elements(loc2)-1 do					$
	  tv, ob.mushroom2, 20*xpos(loc2(i)),20*ypos(loc2(i))

	cent = intarr(clen*n_cent,clen*n_cent,5)
	cent(0,0,0) = 1
	for i=1,clen-1 do cent(i,0,0) = 5+(i mod 2)
	cent(*,0,1) = -1
	cent(*,0,2) = -1
	for i=0,clen-1 do cent(i,0,3) = xsize+i
	cent(*,0,4) = ysize-1

	if n_cent ge 2 then begin

	  for i = 2, n_cent do begin

	    found = 0
	    while not found do begin
	      ypos = fix(((ysize-1)*randomu(s,1))*.7+.33*ysize)
	      ypos = ypos(0)
	      xpos = fix(xsize*randomu(s,1))
	      if xpos(0) gt xsize/2 then begin
		xpos = xsize 
		dir = 1
	      endif else begin
		xpos = 0
		dir = -1
	      endelse

	      if ypos lt ysize-2 then				$
		if (grid(xpos-dir,ypos) eq 0) then begin
		  loc = where(cent(0,*,4) eq ypos,count)
		  if count eq 0 then found = 1
		endif
	    endwhile

	    if dir eq 1 then begin
	      cent(0,1,0) = 1 
	      for i=1,clen-1 do cent(i,1,0) = 7+(i mod 2)
	    endif else begin
	      cent(0,1,0) = 3
	      for i=1,clen-1 do cent(i,1,0) = 5+(i mod 2)
	    endelse
	    cent(*,1,1) = -dir
	    cent(*,1,2) = -1
	    for i=0,clen-1 do cent(i,1,3) = xpos+dir*i
	    cent(*,1,4) = ypos

	  endfor

	endif

	return

end

