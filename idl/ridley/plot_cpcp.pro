
restore, file='Results_9904/pot.save'
;restore, file='Results_9812/pot.save'
;restore, file='Results_9811/pot.save'

nlines = 1440*2
data = fltarr(15,nlines)
time = dblarr(nlines)

openr,1,'IMF990416.dat'
;openr,1,'IMF981210.dat'
;openr,1,'IMF981105.dat'

line = ''
for i=0,8 do readf,1,line

readf,1,data

close,1

itime = intarr(6)

for i=0,nlines-1 do begin

  itime(*) = data(0:5,i)
  c_a_to_r, itime, rtime

  time(i) = rtime

endfor

nvals = n_elements(time_save)
sw = fltarr(8,nvals)

for i=0,nvals-1 do begin

  dt = abs(time - (time_save(i) - 3600.0))
  loc = where(dt eq min(dt))
  sw(*,i) = data(7:14, loc(0))

  c_r_to_a, itime, time_save(i) - 3600.0

  file = 'run.'+chopr('0'+tostr(itime(0)),2)+$
	chopr('0'+tostr(itime(1)),2)+$
	chopr('0'+tostr(itime(2)),2)+'.'+$
	chopr('0'+tostr(itime(3)),2)+$
	chopr('0'+tostr(itime(4)),2)+$
	chopr('0'+tostr(itime(5)),2)+$
	'/restartOUT/restart.H'

  openr,1,file

  line = ' '
  while strpos(line,'SOLARWIND') lt 0 do readf,1,line

  readf,1,sw(6,i)
  readf,1,sw(7,i)
  readf,1,sw(3,i)
  readf,1,sw(4,i)
  readf,1,sw(5,i)
  readf,1,sw(0,i)
  readf,1,sw(1,i)
  readf,1,sw(2,i)

  close,1

endfor

bz = sw(2,*)
by = sw(1,*)
bx = sw(0,*)
vx = sw(3,*)
btot = sqrt(bz^2+by^2)
den = sw(7,*)

clock = acos(bz/btot)

esw = -vx * btot*(sin(clock/2.0))

setdevice, 'cpcp.ps','p',4,0.9

ppp = 3
space = 0.09

!p.charsize = 1.2

pos_space, ppp, space, sizes, nx = 1

loc1 = where(bz lt -1)
loc2 = where(bz ge -1)

get_position, ppp, space, sizes, 0, pos, /rect
plot, bz(loc1), cpcp_save(0,loc1), psym = 2, pos = pos, $
	ytitle = 'Potential (kV)', xtitle = 'Bz (nT)', $
	xrange = [-25,15], xstyle = 1, title = 'MHD Only'

oplot, bz(loc2), cpcp_save(0,loc2), psym = 4

get_position, ppp, space, sizes, 1, pos, /rect
plot, by(loc1), cpcp_save(0,loc1), psym = 2, pos = pos, $
	/noerase, ytitle = 'Potential (kV)', xtitle = 'By (nT)', $
	xrange = [-15,25], xstyle = 1

oplot, by(loc2), cpcp_save(0,loc2), psym = 4

get_position, ppp, space, sizes, 2, pos, /rect
plot, btot(loc1), cpcp_save(0,loc1), psym = 2, pos = pos, $
	/noerase, ytitle = 'Potential (kV)', xtitle = 'Btot (nT)', $
	xrange = [0,25], xstyle = 1

oplot, btot(loc2), cpcp_save(0,loc2), psym = 4

end
