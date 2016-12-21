pro move_man, mouse, gameover

	common objects, ob
	common positions, grid

	cursor, x, y, 0
	mouse = !err

	if (x*!d.x_size ne ob.oldman.x) or (y*!d.y_size ne ob.oldman.y) and $
	  (x gt 0.0) and (y gt 0.0) and (x lt 1.0) and (y lt 1.0) then begin

	  if x lt 0.02 then x = 0.02
	  if x gt 0.98 then x = 0.98
	  if y lt 0.02 then y = 0.02
	  if y gt 0.3 then y = 0.3

	  x = x * !d.x_size
	  y = y * !d.y_size
	  tv, ob.oldman.image, ob.oldman.x-ob.oldman.xoff, 		$
	    ob.oldman.y-ob.oldman.yoff
	  tv, ob.newman.image, x-ob.newman.xoff, y-ob.newman.yoff

	  ob.oldman.xvel = x - ob.oldman.x
	  ob.oldman.yvel = y - ob.oldman.y
	  ob.oldman.x = x
	  ob.oldman.y = y

	  gx = (x+ob.oldman.xoff)/20
	  gy = (y+ob.oldman.yoff)/20
	  if grid(gx,gy) eq 5 then gameover = 1

	  gx = (x-ob.oldman.xoff)/20
	  gy = (y+ob.oldman.yoff)/20
	  if grid(gx,gy) eq 5 then gameover = 1

	  gx = (x+ob.oldman.xoff)/20
	  gy = (y-ob.oldman.yoff)/20
	  if grid(gx,gy) eq 5 then gameover = 1

	  gx = (x-ob.oldman.xoff)/20
	  gy = (y-ob.oldman.yoff)/20
	  if grid(gx,gy) eq 5 then gameover = 1

	endif

	return

END

