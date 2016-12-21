
	common objects, ob
	common positions, grid

	window, 1, xsize = 20*30, ysize = 20*40
	create_ob
	device, cursor_image=bytarr(16)
	mouse = 0

	plot, [0,1], /nodata, xstyle=5, ystyle=5, pos = [0,0,1,1]

	xsize = float(!d.x_size)
	ysize = float(!d.y_size)

	nshrooms = fix(15*randomu(s,1))+20
	nshrooms = nshrooms(0)
	maxbul = 50
	clen = 20
	n_cent = 1
	bullet = -1
	b_speed = 5
	spider_velocity = 10
	c_slow = 20
	grid_cent, nshrooms, grid, cent, clen, n_cent

	timer = 0
	timer2 = -1
	spidertimer = -1
	s_slow = 1

	gameover = 0

	while not gameover do begin

	  timer2 = fix((timer2 + 1) mod c_slow)
	  if timer2 eq 0 then move_cent, cent

	  if (mouse eq 1) and (n_elements(bullet) gt 1) then begin
	    timer = (timer + 1) mod 25
	    if timer ne 0 then mouse = 0
	  endif else begin
	    timer = 1
	    if mouse eq 0 then timer = -1
	  endelse

	  rand = randomu(s,1)
	  move_bullets, bullet, mouse, maxbul, b_speed, cent,rand

	  if (ob.spider.on eq 1) and (timer2 eq 0) then begin
;	    spidertimer = fix((spidertimer + 1) mod s_slow)
;	    if spidertimer eq 0 then
move_spider, xsize, ysize,s,bullet, spider_velocity
	  endif else if rand(0) ge 0.999 then ob.spider.on = 1 

	  move_man, mouse, gameover
	  if mouse eq 4 then gameover=1

	  loc = where(cent(0,*,0) ne 0, count)

	  if count eq 0 then begin

	    plot, [0,1], /nodata, xstyle=5, ystyle=5, pos = [0,0,1,1]
	    nshrooms = nshrooms + 5
	    clen = clen + 5
	    n_cent = clen/20 + 1
	    length = clen mod 20
	    if length lt 10 then length = 10
	    bullet = -1
	    b_speed = 5
	    c_slow = c_slow*0.9
	    grid_cent, nshrooms, grid, cent, length, n_cent

	    timer = 0
	    timer2 = -1

	  endif

	endwhile

	device, /cursor_cross
	wdelete
end

