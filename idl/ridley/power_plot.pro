
if n_elements(filein) eq 0 then begin
   filelist = findfile('power*')
   filein = filelist(0)
endif

filein = ask('file to plot',filein)

if n_elements(psfile) eq 0 then psfile = filein+'.ps'
psfile = ask('ps file',psfile)

power_read, filein, time, power, North, nPts, notes

stime = min(time)
etime = max(time)

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

l = where(North eq 1)

maxi = max(power)*1.1

setdevice, psfile, 'p', 5

plot, time(l)-stime, power(l), $
      ytitle = 'Hemispheric Power (GW)', yrange = [0,maxi], $
      xstyle = 1, $
      xtickname=xtickname, xtitle=xtitle, xtickv=xtickv, $
      xminor=xminor, xticks=xtickn, thick = 4, $
      pos = [0.1,0.3,0.9,0.7]

l = where(North eq 0)

oplot, time(l)-stime, power(l), linestyle = 2, thick = 4

dt = (etime-stime)*0.05
dy = maxi*0.1

oplot, [dt,dt+dt], maxi-[dy,dy], thick = 4
oplot, [dt,dt+dt], maxi-2*[dy,dy], thick = 4, linestyle = 2

xyouts, dt+dt*1.5, maxi-dy, 'North'
xyouts, dt+dt*1.5, maxi-2*dy, 'South'

c_r_to_a, itime, stime
itime([3,4,5]) = 0
c_a_to_r, itime, basetime

nDays = fix((etime-basetime)/(24.0*3600.0))+1
for i=1,nDays-1 do $
   oplot, basetime-stime+double(i)*3600.0*24.0*[1.0,1.0], $
          [0,maxi*1.2], linestyle = 1

closedevice

end
