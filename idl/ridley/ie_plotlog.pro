
files = findfile('IE_*.log')

ie_readlog, files, time, data, Vars

stime = min(time)
etime = max(time)

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

setdevice, 'cpcp.ps','p', 5

plotdumb

if (n_elements(maxi) eq 0) then maxi = max(data(9:10,*))
maxi = float(ask('maximum value for the CPCP plot',string(maxi)))

yr = [0,maxi]
plot, time-stime, data(9,*), xstyle = 1, xrange = [btr,etr], $
  xtickname = xtickname, xminor=xminor, xticks=xtickn, $
  xtickv = xtickv, ystyle = 1, pos = [0.1,0.3, 0.9,0.7], $
  xtitle = xtitle, ytitle = 'CPCP (kV)', yrange = yr, /noerase, $
      thick = 4

oplot, time-stime, data(10,*), linestyle = 1, thick = 4, color = 100

dy = yr(1)*0.1
dx = (max(time)-min(time))*0.1
oplot, [dx, 2*dx], yr(1)-[dy,dy], thick = 4
xyouts, 2.05*dx, yr(1)-dy, 'North'
oplot, [dx, 2*dx], yr(1)-1.5*[dy,dy], thick = 4, linestyle = 1, color = 100
xyouts, 2.05*dx, yr(1)-1.5*dy, 'South', color = 100

if (etr gt 24.0*3600.0) then begin

  c_r_to_a, itime, stime
  itime(3:5) = 0
  c_a_to_r, itime, basetime

  ndays = round((max(time) - basetime)/24.0/3600.0)
  for i=1,ndays do begin
     t = basetime + i*24.0*3600.0 - stime
     oplot, [t,t], [0,maxi], linestyle = 2
  endfor

endif

spawn, 'pwd', pwd

xyouts, 0.0,0.0, /norm, pwd(0), charsize = 0.75

closedevice

end
