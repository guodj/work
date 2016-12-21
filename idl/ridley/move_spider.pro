pro move_spider, xs, ys, s, bullet, spider_velocity

  common objects, ob

  ymax = 0.4*ys

  if ob.spider.ypos eq -1 then begin

    rand = randomu(s,2)
    ob.spider.ypos = (ymax-10)*rand(0)+10
    if rand(1) ne 0.5 then sign = (rand(1)-0.5)/abs(rand(1)-0.5) 	$
    else sign = 1.0
    ang = !pi/3.0*rand(0) + !pi/12.0
    if rand(0) lt 0.5 then begin
      ob.spider.xpos = 0.0
      ob.spider.slope = sign*tan(ang)
      ob.spider.velocity = +spider_velocity
      ob.spider.lastx = ob.spider.xpos - ob.spider.xsize/2
      ob.spider.lasty = ob.spider.ypos - ob.spider.ysize/2
    endif else begin
      ob.spider.xpos = xs
      ob.spider.slope = sign*tan(ang)
      ob.spider.velocity = -spider_velocity
      ob.spider.lastx = ob.spider.xpos - ob.spider.xsize/2
      ob.spider.lasty = ob.spider.ypos - ob.spider.ysize/2
    endelse

  endif else tv, ob.spider.overspider, ob.spider.lastx, ob.spider.lasty

  lx = ob.spider.lastx + ob.spider.velocity
  ly = ob.spider.lasty + ob.spider.velocity*ob.spider.slope
  ob.spider.xpos = ob.spider.xpos + ob.spider.velocity
  ob.spider.ypos = ob.spider.ypos + ob.spider.velocity*ob.spider.slope

  sxmin = ob.spider.xpos - ob.spider.xsize
  symin = ob.spider.ypos - ob.spider.ysize
  sxmax = ob.spider.xpos + ob.spider.xsize
  symax = ob.spider.ypos + ob.spider.ysize

  hit = 0
  if bullet(0) ne -1 then begin
    locx = where((bullet(*,0) ge sxmin) and (bullet(*,0) le sxmax), count)
    if count gt 0 then begin
      locy = where((bullet(locx,1) ge symin) and 			$
		   (bullet(locx,1) le symax), county)
      if county gt 0 then begin
        tv, ob.dark_bullet, bullet(locx(locy(0)),0)-2, 			$
			    bullet(locx(locy(0)),1)-4
        bullet(locx(locy(0)),0) = -50
        bullet(locx(locy(0)),1) = -50
        hit = 1
        ob.spider.xpos = -100
      endif
    endif
  endif

  if hit eq 0 then begin

    if (ob.spider.ypos gt ymax) or (ob.spider.ypos lt 5.0) then 	$
	ob.spider.slope = -ob.spider.slope

    if lx lt 0 then lx = 0
    if lx+ob.spider.xsize ge xs then lx = xs-ob.spider.xsize
    if ly lt 0 then ly = 0
    if ly+ob.spider.ysize ge ys then ly = ys-ob.spider.ysize

    ob.spider.lastx = lx
    ob.spider.lasty = ly

    ob.spider.overspider = tvrd(lx,ly,ob.spider.xsize,ob.spider.ysize)
    back = ob.spider.overspider
    case (ob.spider.body) of

      1 : begin
        over = ob.spider.shape1
        loc = ob.spider.loc1
        ob.spider.body = 2
      end

      2 : begin
        over = ob.spider.shape2
        loc = ob.spider.loc2
        ob.spider.body = 1
      end

    endcase

    back(loc) = over(loc)
    tv, back, lx, ly      

  endif

  if (ob.spider.xpos lt 0) or 			$
     (ob.spider.xpos gt xs) then begin

    ob.spider.on = 0
    if hit eq 0 then 				$
	tv, ob.spider.overspider, ob.spider.lastx, ob.spider.lasty
    ob.spider.ypos = -1

  endif

  return

end

