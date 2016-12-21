
restore, file='Results_Weimer/9904/wei_9904.save'
rms_wei = reform(rms_save(0,*))

restore, file='Results_TD/mhd_9904_td.save'
rms_td = reform(rms_save(0,*))

restore, file='Results_9904/mhd_9904.save'

stime = min(time_save)
etime = max(time_save)

time_axis, stime, etime, s_time_range, e_time_range, 	$
	xtickname, xtitle, xtickvalue, xminor, xtickn

yrange = [0.0,max(rms_save)]

setdevice, 'rms_9904.ps','p'

ppp = 4
space = 0.01
pos_space, ppp, space, sizes, nx = 1

loc = where(north_save eq 1 and f13_save eq 1)

xtickname2 = strarr(60)+' '

get_position, ppp, space, sizes, 0, pos, /rect
plot, time_save(loc)-stime, rms_save(0,loc), yrange = yrange, $
	xrange = [s_time_range, e_time_range], xtickname = xtickname2, $
        xtickv=xtickvalue, xticks= xtickn, xminor = xminor,           $
        max_val = 10.0, title = 'RMS', pos = pos, ystyle = 1, psym = -5

oplot, time_save(loc)-stime, rms_save(1,loc), linestyle = 1, psym = -2
oplot, time_save(loc)-stime, rms_wei(loc), linestyle = 2, psym = -4
oplot, time_save(loc)-stime, rms_td(loc), linestyle = 3, psym = -1

dx = (e_time_range - s_time_range) * 0.02
dy = 0.05

y = dy

x = 1.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 1, psym = -2
xyouts, x+dx*2.5, dy, 'F13 North'

x = 11.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 0, psym = -5
xyouts, x+dx*2.5, dy, 'MHD Steady'

x = 21.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 3, psym = -1
xyouts, x+dx*2.5, dy, 'MHD T-A'

x = 31.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 2, psym = -4
xyouts, x+dx*2.5, dy, 'Weimer'

loc = where(north_save eq 1 and f13_save eq 0, count)

if (count gt 1) then begin

  pn = 1
  get_position, ppp, space, sizes, pn, pos, /rect
  plot, time_save(loc)-stime, rms_save(0,loc), yrange = yrange, $
	xrange = [s_time_range, e_time_range], xtickname = xtickname2, $
        xtickv=xtickvalue, xticks= xtickn, xminor = xminor,           $
        max_val = 10.0, title = ' ', pos = pos, /noerase, ystyle = 1, psym = -5

  oplot, time_save(loc)-stime, rms_save(1,loc), linestyle = 1, psym = -2
  oplot, time_save(loc)-stime, rms_wei(loc), linestyle = 2, psym = -4
  oplot, time_save(loc)-stime, rms_td(loc), linestyle = 3, psym = -1

x = 1.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 1, psym = -2
xyouts, x+dx*2.5, dy, 'F14 North'

x = 11.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 0, psym = -5
xyouts, x+dx*2.5, dy, 'MHD Steady'

x = 21.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 3, psym = -1
xyouts, x+dx*2.5, dy, 'MHD T-A'

x = 31.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 2, psym = -4
xyouts, x+dx*2.5, dy, 'Weimer'

endif else pn = 0

loc = where(north_save eq 0 and f13_save eq 1)

get_position, ppp, space, sizes, pn+1, pos, /rect
plot, time_save(loc)-stime, rms_save(0,loc), yrange = yrange, $
	xrange = [s_time_range, e_time_range], xtickname = xtickname2, $
        xtickv=xtickvalue, xticks= xtickn, xminor = xminor,           $
        max_val = 10.0, title = ' ', pos = pos, /noerase, ystyle = 1, psym = -5

oplot, time_save(loc)-stime, rms_save(1,loc), linestyle = 1, psym = -2
oplot, time_save(loc)-stime, rms_wei(loc), linestyle = 2, psym = -4
oplot, time_save(loc)-stime, rms_td(loc), linestyle = 3, psym = -1

x = 1.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 1, psym = -2
xyouts, x+dx*2.5, dy, 'F13 South'

x = 11.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 0, psym = -5
xyouts, x+dx*2.5, dy, 'MHD Steady'

x = 21.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 3, psym = -1
xyouts, x+dx*2.5, dy, 'MHD T-A'

x = 31.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 2, psym = -4
xyouts, x+dx*2.5, dy, 'Weimer'

loc = where(north_save eq 0 and f13_save eq 0)

get_position, ppp, space, sizes, pn+2, pos, /rect
plot, time_save(loc)-stime, rms_save(0,loc), yrange = yrange, $
	xrange = [s_time_range, e_time_range], xtickname = xtickname, $
        xtickv=xtickvalue, xticks= xtickn, xminor = xminor,           $
        max_val = 10.0, title = ' ', pos = pos, /noerase, xtitle = xtitle, $
	ystyle = 1, psym = -5

oplot, time_save(loc)-stime, rms_save(1,loc), linestyle = 1, psym = -2
oplot, time_save(loc)-stime, rms_wei(loc), linestyle = 2, psym = -4
oplot, time_save(loc)-stime, rms_td(loc), linestyle = 3, psym = -1

x = 1.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 1, psym = -2
xyouts, x+dx*2.5, dy, 'F14 South'

x = 11.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 0, psym = -5
xyouts, x+dx*2.5, dy, 'MHD Steady'

x = 21.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 3, psym = -1
xyouts, x+dx*2.5, dy, 'MHD T-A'

x = 31.0*dx
oplot, [x, x+2*dx],[y,y], linestyle = 2, psym = -4
xyouts, x+dx*2.5, dy, 'Weimer'

closedevice

end
