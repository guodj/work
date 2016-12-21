pro get_mouse_command, pos, command, ret_pos

  out = 0

  npts = n_elements(pos(*,0))

  while not out do begin

    cursor,x,y,2, /norm
    mouse = !err

    if mouse gt 0 then begin

      l = where( (x gt pos(*,0)) and (x lt pos(*,2)) and 		$
		 (y gt pos(*,1)) and (y lt pos(*,3)), count)

      if count gt 0 then begin

	command = l(0) + (mouse-1)*npts
	ret_pos = fltarr(2)
        ret_pos(0) = (x-pos(l(0),0))/(pos(l(0),2)-pos(l(0),0))
        ret_pos(1) = (y-pos(l(0),1))/(pos(l(0),3)-pos(l(0),1))

	out = 1

      endif
  
      clear_mouse

    endif

  endwhile

  return

end
