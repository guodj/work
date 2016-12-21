
file = 'newbIOy_n001250.dat.save'
file = ask('file name (can use *)',file)

filelist = findfile(file)
nfiles = n_elements(filelist)

psfile = ''
psfile = ask('ps filename (null string for window)',psfile)

min_size = '0.25'
min_size = float(ask('cell size for rebinning (in Re)',min_size))

center_x = '-10.0'
center_x = float(ask('center X (in Re)',center_x))

center_y = '0.0'
center_y = float(ask('center Y (in Re)',center_y))

rsize = '20.0'
rsize = float(ask('maximum distance to plot along Y/Z axis (in Re)',rsize))

aspect = '1.5'
aspect = float(ask('aspect',aspect))

rBody = '3.01'
rBody = float(ask('Radius of Body',rBody))

print, rbody

delta_y = rsize
delta_x = aspect*rsize

ex = center_x + delta_x
sx = center_x - delta_x
ey = center_y + delta_y
sy = center_y - delta_y

grid = 'y'
grid = strpos(ask('whether you would like a grid or not (y/n)',grid),'y')+1

if nfiles gt 4 then ppp = '4' else ppp = tostr(nfiles)
ppp = fix(ask('plots per page',ppp))

for n=0,nfiles-1 do begin

  file = filelist(n)

  print, "Reading file ", file

  openr,1,file

  done = 0

  line = ''

  while not done do begin

    readf,1,line

    if strpos(line,'NUMERICAL') gt -1 then begin
      nvars = 0
      nlines = 0
      readf,1, nvars
      readf,1, nlines
    endif

    if strpos(line,'VARIABLE') gt -1 then begin
      if n_elements(nvars) eq 0 then begin
        print, 'File is in the wrong order, NUMERICAL VALUES must come'
        print, 'before VARIABLE LIST'
        stop
      endif else begin
        vars = strarr(nvars+1)
	vars(nvars) = 'Beta'
        for i=0,nvars-1 do begin
          readf,1,format="(I5,a)",j,line
          vars(i) = line
        endfor

        if (n eq 0) then begin
          for i=0,nvars-1 do print, string(i+1,format='(I2)'),'. ',vars(i)
	  print, string(i+1,format='(I2)'),'. Beta'
          var_to_plot = fix(ask('variable to plot','4'))-1
        endif

      endelse
    endif

    if strpos(line,'BEGIN') gt -1 then done = 1

  endwhile

  data = fltarr(nvars, nlines)
  tmp = fltarr(nvars)

  format = '('+tostr(nvars)+'E16.8)'

  for i=0,nlines-1 do begin
    readf,1,tmp, format=format
    data(*,i) = tmp
  endfor

  close,1

  x = reform(data(0,*))
  y = reform(data(1,*))
  z = reform(data(2,*))

  if (var_to_plot lt nvars) then begin
    input_var = alog10(reform(data(var_to_plot,*)))
  endif else begin
    tmp_bx = reform(data(11,*)) * 1.0e-9
    tmp_by = reform(data(12,*)) * 1.0e-9
    tmp_bz = reform(data(13,*)) * 1.0e-9
    tmp_b2  = tmp_bx^2 + tmp_by^2 + tmp_bz^2 + 1.0e-9^2.0 
    tmp_p  = reform(data(4,*)) * 1.0e-9
    input_var = !pi*4.0e-7*tmp_p/tmp_b2
  endelse

  ytitle = 'Y (Re)'
  if (max(y) eq 0) then begin
    yfac = 1.0
    y = z
    ytitle = 'Z (Re)'
    input_bx = reform(data(11,*))
    input_by = reform(data(13,*))
  endif else begin
    yfac = -1.0
    input_bx = reform(data(5,*))
    input_by = reform(data(6,*))
  endelse

  print, ytitle

  print, "=> Triangulating"

  triangulate, x, y, tr

  if (max(x) lt ex) then begin
    diff = ex - max(x)
    ex = max(x)
    sx = ex - 2.0*delta_x
  endif

  trigridded_input =trigrid(x,y,input_var,tr,[min_size,min_size],	$
	[sx,sy,ex,ey])

  trigridded_bx =trigrid(x,y,input_bx,tr,[min_size,min_size],	$
	[sx,sy,ex,ey])

  trigridded_by =trigrid(x,y,input_by,tr,[min_size,min_size],	$
	[sx,sy,ex,ey])

  if n eq 0 then begin
    nx = n_elements(trigridded_input(*,0))
    ny = n_elements(trigridded_input(0,*))
    new_v = fltarr(nfiles, nx, ny)
    new_bx = fltarr(nfiles, nx, ny)
    new_by = fltarr(nfiles, nx, ny)
  endif

  new_v(n,*,*) = trigridded_input
  new_bx(n,*,*) = trigridded_bx
  new_by(n,*,*) = trigridded_by

endfor

new_bt = sqrt(new_bx^2+new_by^2)
new_bx = new_bx/new_bt
new_by = new_by/new_bt

print, "=> Trigridding X"
new_x = trigrid(x,y,x,tr,[min_size,min_size],				$
	[sx,sy,ex,ey])
print, "=> Trigridding Y"
new_y = trigrid(x,y,y,tr,[min_size,min_size],				$
	[sx,sy,ex,ey])

; variables are:
; 0 - x
; 1 - y
; 2 - z
; 3 - density
; 4 - pressure
; 5 - Ux
; 6 - Uy
; 7 - Uz
; 8 - B1x
; 9 - B1y
;10 - B1z
;11 - Bx
;12 - By
;13 - Bz
;14 - jx
;15 - jy
;16 - jz
;17 - ex
;18 - ey
;19 - ez
;20 - t
;21 - 
;22 -
;23 - 

if strlen(psfile) eq 0 then window, 0 else begin
  if ppp gt 1 then landport = 'p' else landport = 'l'
  setdevice, psfile, landport, 4, 0.95
endelse

space = 0.02
pos_space, ppp, space, sizes, ny=ppp

if min(new_v) lt 0 then begin
  readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
endif else begin
  readct,ncolors, getenv("IDL_EXTRAS")+"jim6.ct"
endelse
clevels = (ncolors)*findgen(30)/29.0 + 1

maxi = max(new_v)
mini = min(new_v)

print, 'minimum value : ',mini
print, 'maximum value : ',maxi
mini = float(ask('new minimum',string(mini)))
maxi = float(ask('new maximum',string(maxi)))

mini_save = mini
maxi_save = maxi

n_trace = 20

if (strpos(ytitle, 'Z') gt -1) then begin

  pos_old_x = findgen(n_trace)*(nx-1)/float(n_trace-1) - 20
  pos_old_y = pos_old_x*0 + ny/2

  loc = where(new_x(*,0) eq 0)

  others_x = findgen(31)*12.0/30*4 - 28.0 + loc(0)
  others_y = others_x*0.0 + ny/2 + 20
  pos_old_x = [pos_old_x,others_x]
  pos_old_y = [pos_old_y,others_y]

  others_x = findgen(31)*12.0/30*4 - 28.0 + loc(0)
  others_y = others_x*0.0 + ny/2 - 20
  pos_old_x = [pos_old_x,others_x]
  pos_old_y = [pos_old_y,others_y]

endif else begin

  pos_old_x = fltarr(n_trace+1) + 7*nx/8
  pos_old_y = findgen(n_trace+1)*(ny)/n_trace

  others_x = fltarr(n_trace+1) + 1*nx/8
  others_y = findgen(n_trace+1)*(ny)/n_trace
  pos_old_x = [pos_old_x,others_x]
  pos_old_y = [pos_old_y,others_y]

  locx = where(new_x(*,0) eq 0.0)
  locy1 = where(new_y(0,*) eq -10)
  locy2 = where(new_y(0,*) eq 10)

  others_x = fltarr(n_trace) + locx(0)
  others_y = findgen(n_trace)*(locy2(0) - locy1(0))/(n_trace) + locy1(0)

  pos_old_x = [pos_old_x,others_x]
  pos_old_y = [pos_old_y,others_y]

endelse

n_trace = n_elements(pos_old_x)

pos_new_x = pos_old_x
pos_new_y = pos_old_y

pos_old_x_neg = pos_old_x
pos_old_y_neg = pos_old_y
pos_new_x_neg = pos_old_x
pos_new_y_neg = pos_old_y

new_plot_pos_x =  $
	interpolate(new_x,pos_old_x,pos_old_y)

new_plot_pos_y =  $
	interpolate(new_y,pos_old_x,pos_old_y)
old_plot_pos_x = new_plot_pos_x
old_plot_pos_y = new_plot_pos_y

new_plot_pos_x_neg =  $
	interpolate(new_x,pos_old_x_neg,pos_old_y_neg)
new_plot_pos_y_neg =  $
	interpolate(new_y,pos_old_x_neg,pos_old_y_neg)
old_plot_pos_x_neg = new_plot_pos_x_neg
old_plot_pos_y_neg = new_plot_pos_y_neg

    active = where(new_plot_pos_x le max(new_x) and   $
                   new_plot_pos_x ge min(new_x) and   $
                   new_plot_pos_y le max(new_y) and   $
                   new_plot_pos_y ge min(new_y),active_count)

    active_neg = where(new_plot_pos_x_neg le max(new_x) and   $
                       new_plot_pos_x_neg ge min(new_x) and   $
                       new_plot_pos_y_neg le max(new_y) and   $
                       new_plot_pos_y_neg ge min(new_y),active_count_neg)



for n = 0, nfiles-1 do begin

  mini = mini_save
  maxi = maxi_save

  pn = n mod ppp

  if pn eq ppp-1 then begin
    xtn = strarr(20) 
    xt = 'X (Re)'
  endif else begin
    xtn = strarr(20)+' '
    xt = ' '
  endelse

  get_position, ppp, space, sizes, pn, pos
  pos(2) = pos(0) + aspect*2*(pos(2)-pos(0))/2.0
  diff = 0.5 - mean(pos([0,2]))
  pos([0,2]) = pos([0,2]) + diff
  if (pos(2) gt 0.97) then begin
    xr = pos(2) - pos(0)
    if (xr gt 0.97) then begin
      ppp = ppp*2
      pos_space, ppp, space, sizes, nx = 2
      get_position, ppp, space, sizes, pn, pos
      pos(2) = pos(0) + aspect*2*(pos(2)-pos(0))/2.0
      diff = 0.5 - mean(pos([0,2]))
      pos([0,2]) = pos([0,2]) + diff
    endif else begin
      pos(2) = 0.97
      pos(0) = pos(2) - xr
    endelse
  endif

  if pn eq 0 then begin
    plotdumb
    polyfill, [0,0,1,1,0],[0,1,1,0,0], color = ncolors
  endif

  levels = findgen(30)/29.0*(maxi-mini) + mini

  data = reform(new_v(n,*,*))

  loc = where(data gt levels(27),count)
  if count gt 0 then data(loc) = levels(27)

  loc = where(data lt levels(1),count)
  if count gt 0 then data(loc) = levels(1)

  contour, data, new_x, new_y,  xrange = [ex,sx], $
	yrange = yfac*[sy,ey], 					$
	levels = levels, c_colors = clevels,  /fill, pos = pos,	$
	color = 0, /noerase, xtitle = xt, ytitle = ytitle,	$
	ystyle = 1, xstyle = 1, xtickname = xtn

  ctlevels = mm(levels)

;  if min(new_v) lt 0 then begin
;    maxi = max(abs(levels))
;    levels = [-maxi*0.75,-maxi*0.5,-maxi*0.25,0.0,maxi*0.25,maxi*0.5,maxi*0.75]
;  endif else levels = levels(indgen(5)*4)

;  contour, new_v(n,*,*), new_x, new_y,  xrange = [ex,sx], 	$
;	yrange = [sy,ey], 					$
;	levels = levels, pos = pos, /follow,				$
;	color = 0, /noerase, xstyle = 5, ystyle = 5

  xyouts, pos(2), pos(3)+0.01, vars(var_to_plot), align = 1.0, /norm, color=0
;  xyouts, pos(0), pos(3)+0.01, title, align = 0.0, /norm, color=0

  i = 0

  while ((active_count gt 0) or (active_count_neg gt 0)) and i lt 750 do begin

    if (active_count gt 0) then begin

      pos_new_x(active) = pos_old_x(active) + $
	0.5*interpolate(reform(new_bx(n,*,*)),$
	pos_old_x(active),pos_old_y(active))
      pos_new_y(active) = pos_old_y(active) + $
	0.5*interpolate(reform(new_by(n,*,*)),$
	pos_old_x(active),pos_old_y(active))

      old_plot_pos_x(active) = $
	interpolate(new_x,pos_old_x(active),pos_old_y(active))
      old_plot_pos_y(active) = $
	interpolate(new_y,pos_old_x(active),pos_old_y(active))
      new_plot_pos_x(active) = $
	interpolate(new_x,pos_new_x(active),pos_new_y(active))
      new_plot_pos_y(active) = $
	interpolate(new_y,pos_new_x(active),pos_new_y(active))

      for jj=0,active_count-1 do begin
        j = active(jj)
        oplot, [old_plot_pos_x(j), new_plot_pos_x(j)], $
		[old_plot_pos_y(j), new_plot_pos_y(j)]
      endfor

    endif

    if (active_count_neg gt 0) then begin

      pos_new_x_neg(active_neg) = pos_old_x_neg(active_neg) - $
	0.5*interpolate(reform(new_bx(n,*,*)), $
	 pos_old_x_neg(active_neg),pos_old_y_neg(active_neg))

      pos_new_y_neg(active_neg) = pos_old_y_neg(active_neg) - $
	0.5*interpolate(reform(new_by(n,*,*)), $
	 pos_old_x_neg(active_neg),pos_old_y_neg(active_neg))

      old_plot_pos_x_neg(active_neg) =  $
	interpolate(new_x,pos_old_x_neg(active_neg),pos_old_y_neg(active_neg))
      old_plot_pos_y_neg(active_neg) =  $
	interpolate(new_y,pos_old_x_neg(active_neg),pos_old_y_neg(active_neg))
      new_plot_pos_x_neg(active_neg) =  $
	interpolate(new_x,pos_new_x_neg(active_neg),pos_new_y_neg(active_neg))
      new_plot_pos_y_neg(active_neg) =  $
	interpolate(new_y,pos_new_x_neg(active_neg),pos_new_y_neg(active_neg))

    endif

    for jj=0,active_count_neg-1 do begin
      j = active_neg(jj)
      oplot, [old_plot_pos_x_neg(j), new_plot_pos_x_neg(j)], $
		[old_plot_pos_y_neg(j), new_plot_pos_y_neg(j)]
    endfor

    r = sqrt(new_plot_pos_x^2 + new_plot_pos_y^2)
    loc = where(r lt rBody, count)
    if (count gt 0) then begin
      new_plot_pos_x(loc) = max(new_x)*2.0
    endif
      
    r = sqrt(new_plot_pos_x_neg^2 + new_plot_pos_y_neg^2)
    loc = where(r lt rBody, count)
    if (count gt 0) then begin
      new_plot_pos_x_neg(loc) = max(new_x)*2.0
    endif
      
    active = where(new_plot_pos_x lt max(new_x) and   $
                   new_plot_pos_x gt min(new_x) and   $
                   new_plot_pos_y lt max(new_y) and   $
                   new_plot_pos_y gt min(new_y),active_count)

    active_neg = where(new_plot_pos_x_neg lt max(new_x) and   $
                       new_plot_pos_x_neg gt min(new_x) and   $
                       new_plot_pos_y_neg lt max(new_y) and   $
                       new_plot_pos_y_neg gt min(new_y),active_count_neg)

    pos_old_x = pos_new_x
    pos_old_y = pos_new_y

    pos_old_x_neg = pos_new_x_neg
    pos_old_y_neg = pos_new_y_neg

    i = i + 1

  endwhile

  theta = findgen(37)/36.0*2.0*!pi
  xcircle = rBody*cos(theta)
  ycircle = rBody*sin(theta)

  loc = where(xcircle lt max(new_x) and   $
              xcircle gt min(new_x) and   $
              ycircle lt max(new_y) and   $
              ycircle gt min(new_y),count)

  if count gt 18 then polyfill, rBody*cos(theta), rBody*sin(theta), color = 0

  if (grid) then oplot, x,y,psym=3, color = 0

  ctpos = [pos(2)+0.01,pos(1), pos(2)+0.03,pos(3)]
  plotct, ncolors, ctpos, ctlevels, vars(var_to_plot), /right, 	$
	color = 0

endfor

closedevice

end

