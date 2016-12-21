pro clear_mouse

  cursor,x,y,0
  mouse = !err

  while mouse ne 0 do begin

    cursor,x,y,0
    mouse = !err

  endwhile

  return

end
