pro mltpoly, image, range, theta, maxran, pos

plot, [-maxran,maxran],[-maxran,maxran], /nodata, /noerase,	$
	xstyle = 5, ystyle = 5, pos = pos

dr = range(1) - range(0)
dt = (theta(1) - theta(0))*!pi/180.0

nxn = n_elements(theta)
locy = where(range gt 90.-maxran, nyn)

if nyn gt 0 then begin

  r = 90.0 - range

  for i=0,nxn-1 do for j=0,nyn-1 do begin

    rl = r(locy(j)) + dr/2.0  
    rh = rl - dr
    if rl gt maxran then rl = maxran

    tl = theta(i)*!pi/180.0 - dt/2.0
    tm = tl + dt/2.0
    th = tl + dt

    x = [rl*cos(th),rh*cos(th),rh*cos(tm),rh*cos(tl),		$
	 rl*cos(tl),rl*cos(tm),rl*cos(th)]

    y = [rl*sin(th),rh*sin(th),rh*sin(tm),rh*sin(tl),		$
	 rl*sin(tl),rl*sin(tm),rl*sin(th)]

    polyfill, x, y, color = image(i,locy(j))

  endfor

  m = min(r)

  if m gt 0.0 then begin

    loc = where(r eq m)
    im  = mean(image(*,loc(0))) 

    t = findgen(37)*2.0*!pi/36.0
    x = m*cos(t)
    y = m*sin(t)

    polyfill, x,y, color = im

  endif

endif

  return

end
