
min_size = 0.25
rsize = 40.0

print,' variables are:'
print,' 0 - x'
print,' 1 - y'
print,' 2 - z'
print,' 3 - density'
print,' 4 - pressure'
print,' 5 - Ux'
print,' 6 - Uy'
print,' 7 - Uz'
print,' 8 - B1x'
print,' 9 - B1y'
print,'10 - B1z'
print,'11 - Bx'
print,'12 - By'
print,'13 - Bz'
print,'14 - jx'
print,'15 - jy'
print,'16 - jz'
print,'17 - ex'
print,'18 - ey'
print,'19 - ez'
print,'20 - t'
print,'21 - '
print,'22 -'
print,'23 - '
varnum = fix(ask('variable to plot','3'))

yz = ask('whether to plot y=0 or z=0 (y/z)','z')
if strlen(yz) gt 1 then yz = strmid(yz,0,1)

filelist = findfile('bIO'+yz+'NEW_n*.save')
nfiles = n_elements(filelist)

for i=0,nfiles-1 do print, tostr(i),'. ',filelist(i)
filenum = fix(ask('file to plot',tostr(nfiles-1)))
if filenum gt nfiles-1 then filenum = nfiles-1
if filenum lt 0 then filenum = 0
file = filelist(filenum)

filelist2 = findfile('../*/'+file)
nfiles2 = n_elements(filelist2)

for i=0,nfiles2-1 do print, tostr(i),'. ',filelist2(i)
filenum2 = fix(ask('file to plot',tostr(nfiles2-1)))
if filenum2 gt nfiles2-1 then filenum2 = nfiles2-1
if filenum2 lt 0 then filenum2 = -1

if filenum2 gt -1 then begin

  file2 = filelist2(filenum2)
  plotall = ask('whether you want to plot difference only (y/n)','y')
  if strpos(plotall,'y') gt -1 then plotall = 0 else plotall = 1
  streamlines = 0

endif else begin

  plotall = 0
  streamlines = ask('whether you want to stream lines (y/n)','y')
  if strpos(streamlines,'y') gt -1 then streamlines = 1 else streamlines = 0

endelse
 
rsize = float(ask('maximum range in y/z (in Re)',string(rsize)))

psfile = ask('Enter ps filename (null for screen)','')

if strlen(psfile) gt 0 then setdevice, psfile, 'p', 4, 0.9

readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
clevels = (ncolors-2)*findgen(21)/20.0 + 1

if plotall then ppp = 6 else ppp = 2
space = 0.02
pos_space, ppp, space, sizes, nx = 2

pn = 0
if plotall then pn = 4

get_position, ppp, space, sizes, pn, pos
pos(2) = pos(2) + 0.5*(pos(2)-pos(0))
diff = 0.5 - mean(pos([0,2]))
pos([0,2]) = pos([0,2]) + diff

restore,file

x = reform(data(0,*))
y = reform(data(1,*))
z = reform(data(2,*))

xtitle = 'X (Re)'
if (max(y) eq 0) then begin
  y = z
  ytitle = 'Z (Re)'
  if streamlines then begin
    xvar = reform(data(11,*))
    yvar = reform(data(13,*))
  endif

  nstreams = 20
  x0 = 20.0 - (findgen(nstreams)+0.5)*3.0*rsize/float(nstreams)
  y0 = x0*0.0 + rsize

endif else begin

  ytitle = 'Y (Re)'
  if streamlines then begin
    xvar = reform(data(5,*))
    yvar = reform(data(6,*))
  endif

  nstreams = 14
  y0 = (findgen(nstreams)-float(nstreams/2)+0.5)*rsize/float(nstreams/2)
  x0 = y0*0.0+20.0

endelse

triangulate, x, y, tr

variable = reform(data(varnum,*))

new_v = trigrid(x,y,variable,tr,[min_size,min_size],			$
	[20.-3.0*rsize,-rsize,20,rsize])
new_x = trigrid(x,y,x,tr,[min_size,min_size],				$
	[20.-3.0*rsize,-rsize,20,rsize])
new_y = trigrid(x,y,y,tr,[min_size,min_size],				$
	[20.-3.0*rsize,-rsize,20,rsize])

if streamlines then begin

  xstream = trigrid(x,y,xvar,tr,[min_size,min_size],			$
	[20.-3.0*rsize,-rsize,20,rsize])
  ystream = trigrid(x,y,yvar,tr,[min_size,min_size],			$
	[20.-3.0*rsize,-rsize,20,rsize])
  smag = sqrt(xstream^2 + ystream^2)

endif

if filenum2 gt -1 then begin

  restore,file2

  x2 = reform(data(0,*))
  y2 = reform(data(1,*))
  z2 = reform(data(2,*))

  if (max(y2) eq 0) then y2 = z2

  triangulate, x2, y2, tr2

  variable2 = reform(data(varnum,*))

  new_v1 = new_v
  new_v2 = trigrid(x2,y2,variable2,tr2,[min_size,min_size],		$
	[20.-3.0*rsize,-rsize,20,rsize])

  new_v = new_v2 - new_v1

endif

maxi = max(new_v)
mini = min(new_v)

maxi = float(ask('maximum value for levels',string(maxi)))
mini = float(ask('minimum value for levels',string(mini)))

levels = findgen(21)/20.0*(maxi-mini) + mini

loc = where(new_v gt maxi,count)
if count gt 0 then new_v(loc) = maxi
loc = where(new_v lt mini,count)
if count gt 0 then new_v(loc) = mini

plotdumb
polyfill, [0,0,1,1,0],[0,1,1,0,0], color = ncolors

contour, new_v, new_x, new_y,  xrange = [20,20-3.0*rsize], 		$
	yrange = [-rsize,rsize], 					$
	levels = levels, c_colors = clevels,  /cell_fill, pos = pos,	$
	color = 0, /noerase, xtitle = xtitle, ytitle = ytitle

levels = levels(indgen(5)*4)

if not streamlines then 						$
  contour, new_v, new_x, new_y,  xrange = [20,20-3.0*rsize], 		$
	yrange = [-rsize,rsize], 					$
	levels = levels, pos = pos, /follow,				$
	color = 0, /noerase

if (not plotall) then 							$
	xyouts, pos(2), pos(3)+0.01, 'test', align = 1.0, /norm, color=0

theta = findgen(37)/36.0*2.0*!pi
polyfill, 3.0*cos(theta), 3.0*sin(theta), color = 0

ctpos = [pos(2)+0.02,pos(1),pos(2)+0.04,pos(3)]

title = ''
plotct, ncolors, ctpos, [mini,maxi], title, /right, 		$
	color = 0

if plotall then begin

  maxi = max([new_v1,new_v2])
  mini = min([new_v1,new_v2])

  maxi = float(ask('maximum value for levels',string(maxi)))
  mini = float(ask('minimum value for levels',string(mini)))

  for i=0,1 do begin

    pn = i*2

    get_position, ppp, space, sizes, pn, pos
    pos(2) = pos(2) + 0.5*(pos(2)-pos(0))
    diff = 0.5 - mean(pos([0,2]))
    pos([0,2]) = pos([0,2]) + diff

    if i eq 0 then new_v = new_v1 else new_v = new_v2

    loc = where(new_v gt maxi,count)
    if count gt 0 then new_v(loc) = maxi
    loc = where(new_v lt mini,count)
    if count gt 0 then new_v(loc) = mini

    levels = findgen(21)/20.0*(maxi-mini) + mini

    contour, new_v, new_x, new_y,  xrange = [20,20-3.0*rsize], 		$
	yrange = [-rsize,rsize], 					$
	levels = levels, c_colors = clevels,  /cell_fill, pos = pos,	$
	color = 0, /noerase, xtickname = strarr(20)+' ', xstyle = 1,	$
	ystyle = 1, ytitle = ytitle

    levels = levels(indgen(5)*4)

    contour, new_v, new_x, new_y,  xrange = [20,20-3.0*rsize], 		$
	yrange = [-rsize,rsize], 					$
	levels = levels, pos = pos, /follow,				$
	color = 0, /noerase, xstyle = 5,				$
	ystyle = 5

    if i eq 0 then 							$
	xyouts, pos(2), pos(3)+0.01, 'test', align = 1.0, /norm, color=0

    polyfill, 3.0*cos(theta), 3.0*sin(theta), color = 0

    ctpos = [pos(2)+0.02,pos(1),pos(2)+0.04,pos(3)]

    title = ''
    plotct, ncolors, ctpos, [mini,maxi], title, /right, 		$
	color = 0

  endfor

endif

if streamlines then begin

  dr = 0.25

  nx = n_elements(new_x(*,0))
  ny = n_elements(new_y(0,*))
  minx = min(new_x)
  maxx = max(new_x)
  xrange = maxx - minx
  miny = min(new_y)
  maxy = max(new_y)
  yrange = maxy - miny

  plot, [0,nx],[0,ny], /nodata, xstyle = 5, ystyle = 5, pos = pos,	$
	/noerase, xrange = [nx,0]

  for i=0,nstreams-1 do begin

    x = nx*(x0(i)-minx)/xrange
    y = ny*(y0(i)-miny)/yrange

    done = 0
    direct = 1.0
    loop = 0
    
    while not done do begin

      if (x lt 0.0) or (x gt nx) then done = 1
      if (y lt 0.0) or (y gt ny) then done = 1

      real_x = x*xrange/nx + minx
      real_y = y*yrange/ny + miny
      real_r = (real_x^2 + real_y^2)
      if (real_r lt 10) then begin
        if loop eq 1 then done = 1 else begin
          loop = 1
          direct = -1.0*direct
          x = nx*(x0(i)-minx)/xrange
          y = ny*(direct*y0(i)-miny)/yrange
        endelse
      endif

      if not done then begin

        sx = interpolate(xstream,x,y)
        sy = interpolate(ystream,x,y)
        sm = interpolate(smag,x,y)

        x2 = x + direct*dr*sx/sm
        y2 = y + direct*dr*sy/sm

        oplot, [x,x2], [y,y2], color = 0

        x = x2
        y = y2

      endif

    endwhile

  endfor

endif

closedevice

end

