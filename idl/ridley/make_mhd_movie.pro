
pro ehandler, ev
  widget_control, /destroy, ev.top
end


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

file_list = findfile('*.save')
nfiles = n_elements(file_list)

print, '1. d B^2 / dt'
print, '2. alog10(density)'
print, '3. pressure'
type = fix(ask('thing to plot','1'))

if (type eq 1) then begin
  readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
  levels = 1600.0*(findgen(21)-10.0)/10.0
endif

if (type eq 2) then begin
  readct,ncolors, getenv("IDL_EXTRAS")+"white_red.ct"
  levels = 6*findgen(21)/20.0-3.0
endif

if (type eq 3) then begin
  readct,ncolors, getenv("IDL_EXTRAS")+"white_red.ct"
  levels = 5*findgen(21)/20.0-4.0
endif

;device, decomposed=0

gifpix = 400

set_plot, 'Z'
device, set_res=[1.5*gifpix,gifpix], Z=0
images = intarr(1.5*gifpix,gifpix, nfiles)

clevels = reverse((ncolors-2)*findgen(21)/20.0 + 1)

ppp = 2
space = 0.02
pos_space, ppp, space, sizes
get_position, ppp, space, sizes, 0, pos
pos(2) = pos(2) + 0.5*(pos(2)-pos(0))
diff = 0.5 - mean(pos([0,2]))
pos([0,2]) = pos([0,2]) + diff

for i=0,nfiles - 1 do begin

  t = i*30
  hour = t/60
  minu = t mod 60

  time = '00'+chopr('0'+tostr(hour),2)+':'+chopr('0'+tostr(minu),2)+' UT'

  print, 'reading file number ',i+1,' out of ',nfiles,' files.'

  file = file_list(i)

  restore,file

  if i eq 0 then begin

    b1x = reform(data(9,*))
    b1y = reform(data(10,*))
    b1z = reform(data(11,*))

    b1 = sqrt(b1x^2 + b1y^2 + b1z^2)

    x = reform(data(0,*))
    z = reform(data(2,*))

    triangulate, x, z, tr

  endif else begin

    b1x_p = b1x
    b1y_p = b1y
    b1z_p = b1z
    b1_p = b1
    x_p = x
    z_p = z
    tr_p = tr

    restore,file

    b1x = reform(data(9,*))
    b1y = reform(data(10,*))
    b1z = reform(data(11,*))
    b1 = sqrt(b1x^2 + b1y^2 + b1z^2)

    density = alog10(reform(data(3,*)))
    pressure = alog10(reform(data(4,*)))

    x = reform(data(0,*))
    z = reform(data(2,*))
    triangulate, x, z, tr

    bx = reform(data(11,*)) - b1x
    by = reform(data(12,*)) - b1y
    bz = reform(data(13,*)) - b1z
    b = sqrt(bx^2 + by^2 + bz^2)

    if type eq 1 then begin
      new_b = trigrid(x,z,b,tr,[0.5,0.5],[-40,-40,20,40])
      new_b1 = trigrid(x,z,b1,tr,[0.5,0.5],[-40,-40,20,40])
      old_b1 = trigrid(x_p,z_p,b1_p,tr_p,[0.5,0.5],[-40,-40,20,40])
    endif 

    if type eq 2 then new_r =trigrid(x,z,density,tr,[0.5,0.5],[-40,-40,20,40])
    if type eq 3 then new_p =trigrid(x,z,pressure,tr,[0.5,0.5],[-40,-40,20,40])

    new_x = trigrid(x,z,x,tr,[0.5,0.5],[-40,-40,20,40])
    new_z = trigrid(x,z,z,tr,[0.5,0.5],[-40,-40,20,40])

    if type eq 1 then begin

      dbtwo = new_b*(new_b1 - old_b1) + (new_b1^2 - old_b1^2)

      loc = where(dbtwo gt max(levels)*0.95, count)
      if count gt 0 then dbtwo(loc) = max(levels)*0.95
      loc = where(dbtwo lt min(levels), count)
      if count gt 0 then dbtwo(loc) = min(levels)

      variable = dbtwo

    endif

    if type eq 2 then variable = new_r
    if type eq 3 then variable = new_p

    plotdumb
    polyfill, [0,0,1,1,0],[0,1,1,0,0], color = 255

    contour, variable, new_x, new_z,  xrange = [20,-40], yrange = [-20,20], $
	levels = levels, c_colors = clevels,  /cell_fill, pos = pos,	$
	color = 0, /noerase

    contour, variable, new_x, new_z,  xrange = [20,-40], yrange = [-20,20], $
	levels = levels, pos = pos, /follow,	$
	color = 0, /noerase

    xyouts, pos(2), pos(3)+0.01, time, align = 1.0, /norm, color=0

    theta = findgen(37)/36.0*2.0*!pi
    polyfill, 3.0*cos(theta), 3.0*sin(theta), color = 0

    bytemap = tvrd()
    images(*,*,i) = bytemap

  endelse

endfor

set_plot, 'X'

for i=1,nfiles-1 do write_gif, 'test.gif', byte(reform(images[*,*,I])), /multiple
write_gif,'test.gif',/close

base = WIDGET_BASE(TITLE = 'Animation Widget')
animate = CW_ANIMATE(base, 1.5*gifpix, gifpix, nfiles)
for i=1,nfiles-1 do CW_ANIMATE_LOAD, animate, FRAME=I, IMAGE=images[*,*,I]

CW_ANIMATE_GETP, animate, pixmap_vect

widget_control, base, /realize 

CW_ANIMATE_RUN, animate

xmanager, 'cw_animate demo', base, event_handler = 'ehandler'

end

