
f107_read, time, f107

itime = [2011,08,01,00,00,00]
c_a_to_r, itime, stime

itime = [2011,09,31,23,59,59]
c_a_to_r, itime, etime

l = where(time ge stime and time le etime)

timer = time(l)
f107r = f107(l)

f107a = f107r
nPts = n_elements(f107a)
for i=0,nPts-1 do begin
   j = l(i)
   js = j-40
   je = j+40
   f107a(i) = mean(f107(js:je))
endfor

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

setdevice, 'f107.ps', 'p', 5

plot, timer-stime, f107r, $
      xstyle = 1, $
      xtickname = xtickname,			$
      xtitle = xtitle,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn, thick = 4, $
      pos = [0.1, 0.3, 0.99, 0.8, 0.99], $
      ytitle = 'F10.7'

oplot, timer-stime, f107a, linestyle = 2, thick = 4

closedevice

end

