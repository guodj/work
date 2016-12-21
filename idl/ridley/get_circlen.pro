pro get_circlen, pos, plotn, eventn

  out = 0

  plotn  = -1
  eventn = -1

  while not out do begin

    cursor,x,y,2,/normal

    mouse = !err

    if mouse eq 1 then begin

      out = 1

      loc = where( (pos(*,0) lt x) and					$
		   (pos(*,1) lt y) and					$
		   (pos(*,2) gt x) and					$
		   (pos(*,3) gt y), count)

      if count gt 0 then begin

	plotn = loc(0)
	eventn = 0

      endif else begin

	dx = pos(0,2) - pos(0,0)

        loc = where( (pos(*,0) lt x) and				$
		     (pos(*,3) lt y) and				$
		     (pos(*,0)+0.2*dx gt x) and				$
		     (pos(*,3)+0.04 gt y), count)

        if count gt 0 then begin

	  plotn = loc(0)
	  eventn = 1

        endif

        loc = where( (pos(*,0)+0.2*dx lt x) and				$
		     (pos(*,3) lt y) and				$
		     (pos(*,0)+0.8*dx gt x) and				$
		     (pos(*,3)+0.04 gt y), count)

        if count gt 0 then begin

	  plotn = loc(0)
	  eventn = 2

        endif

        loc = where( (pos(*,0)+0.8*dx lt x) and				$
		     (pos(*,3) lt y) and				$
		     (pos(*,2) gt x) and				$
		     (pos(*,3)+0.04 gt y), count)

        if count gt 0 then begin

	  plotn = loc(0)
	  eventn = 3

        endif

      endelse

    endif

  endwhile

end

