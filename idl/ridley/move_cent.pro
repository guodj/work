pro move_cent, cent

  common positions, grid
  common objects, ob

  maxx = n_elements(grid(*,0))
  maxy = n_elements(grid(0,*))

  locx = where(cent(0,*,0) ne 0, nc)

  for i=0, nc-1 do begin

    lx = locx(i)

    locy = where(cent(*,lx,0) ne 0,len)

    move = 1
    j = 0
    ly = locy(j)
    vx = cent(ly,lx,1)
    vy = cent(ly,lx,2)
    gx = cent(ly,lx,3)
    gy = cent(ly,lx,4)
    if ((gx+vx lt 0) or (gx+vx ge maxx)) then down = 1	$
    else if grid(gx+vx,gy) ne 0 then down = 1 else down = 0

    if down eq 1 then begin

      if (gy+vy lt 0) or (gy+vy ge maxy) then vy = -1.0*vy

      over = 2
      if grid(gx,gy+vy) eq 0 then over = 0			$
      else begin
        if ((gx+vx ge 0) and (gx+vx lt maxx)) then 		$
	  if grid(gx+vx,gy+vy) eq 0 then over = 1
	if ((gx-vx ge 0) and (gx-vx lt maxx)) then		$
	  if (grid(gx-vx,gy+vy) eq 0) and (over eq 2) then over = -1
      endelse
      if over eq 2 then move = 0

      if (gy+vy lt 0) or (gy+vy ge maxy) then vy = -1.0*vy

      if move eq 1 then begin

        if (len eq 1) and 				$
	   (gx lt maxx) and (gy lt maxy) and		$
	   (gx ge 0) and (gy ge 0) then begin
	  grid(gx,gy) = 0
          tv, ob.dark_shroom, gx*20,gy*20
	endif
        if (over eq 0) or (over eq -1) then vx = -1.0*vx
	gy = gy+vy
	if over ne 0 then gx = gx + vx

        case (cent(ly,lx,0)) of

	  1 : begin
	    cent(ly,lx,0) = 4
	    tv, ob.cent.head2e, gx*20, gy*20
	  end
	  2 : begin
	    cent(ly,lx,0) = 3
	    tv, ob.cent.head1e, gx*20, gy*20
	  end
	  3 : begin
	    cent(ly,lx,0) = 2
	    tv, ob.cent.head2w, gx*20, gy*20
	  end
	  4 : begin
	    cent(ly,lx,0) = 1
	    tv, ob.cent.head1w, gx*20, gy*20
	  end

	endcase

        grid(gx,gy) = 5

      endif 

    endif else begin

      case (cent(ly,lx,0)) of

	1 : begin
	  cent(ly,lx,0) = 2
	  tv, ob.cent.head2w, (gx+vx)*20, gy*20
	end
	2 : begin
	  cent(ly,lx,0) = 1
	  tv, ob.cent.head1w, (gx+vx)*20, gy*20
	end
	3 : begin
	  cent(ly,lx,0) = 4
	  tv, ob.cent.head2e, (gx+vx)*20, gy*20
	end
	4 : begin
	  cent(ly,lx,0) = 3
	  tv, ob.cent.head1e, (gx+vx)*20, gy*20
	end

      endcase

      grid(gx+vx,gy) = 5
      if (len eq 1) and 				$
	 (gx lt maxx) and (gy lt maxy) and		$
	 (gx ge 0) and (gy ge 0) then begin		$
	grid(gx,gy) = 0
	tv, ob.dark_shroom, gx*20,gy*20
      endif
      gx = gx+vx

    endelse

    if move eq 1 then begin

      for j=len-1, 1,-1 do begin

        ly = locy(j)
        if j eq 1 then begin
          case (cent(ly,lx,0)) of
	    5 : cent(ly,lx,0) = 6
	    6 : cent(ly,lx,0) = 5
	    7 : cent(ly,lx,0) = 8
	    8 : cent(ly,lx,0) = 7
          endcase	
        endif else cent(ly,lx,0) = cent(locy(j-1),lx,0)
        cent(ly,lx,1) = cent(locy(j-1),lx,1)
        cent(ly,lx,2) = cent(locy(j-1),lx,2)
        if (j eq len-1) and 						$
	   (cent(ly,lx,3) lt maxx) and (cent(ly,lx,4) lt maxy) and	$
	   (cent(ly,lx,3) ge 0) and (cent(ly,lx,4) ge 0) then begin
	  grid(cent(ly,lx,3),cent(ly,lx,4)) = 0
	  tv, ob.dark_shroom, cent(ly,lx,3)*20,cent(ly,lx,4)*20
        endif
        cent(ly,lx,3) = cent(locy(j-1),lx,3)
        cent(ly,lx,4) = cent(locy(j-1),lx,4)
        gxb = cent(ly,lx,3)
        gyb = cent(ly,lx,4)

        case (cent(ly,lx,0)) of

	  5 : tv, ob.cent.body1w, gxb*20, gyb*20
	  6 : tv, ob.cent.body2w, gxb*20, gyb*20
	  7 : tv, ob.cent.body1e, gxb*20, gyb*20
	  8 : tv, ob.cent.body2e, gxb*20, gyb*20

        endcase

      endfor

      j = 0
      ly = locy(j)
      cent(ly,lx,1) = vx
      cent(ly,lx,2) = vy
      cent(ly,lx,3) = gx
      cent(ly,lx,4) = gy

    endif

  endfor

  return

end

