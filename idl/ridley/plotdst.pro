
logfilename = '/remotehome/gtoth/October/LowresFixedRCM60s/log_all.log'

getlogpro, logfilename, nlogfile, nVars, wlognames, wlog, wlog0, wlog1, wlog2
nPts = n_elements(wlog(*,0))

itime = [2003, 10, 29, 03, 00, 00]
c_a_to_r, itime, basetime
mhdtime = dblarr(nPts) + basetime
mhdtime = mhdtime + wlog(*,1)
mhddst = reform(wlog(*,16))
mhdctc = reform(wlog(*,17))

dstfile = '/r/ridley/MHD/20031029/Data/dst200310.html'

read_dst, dstfile, dst, dsttime

stime = min(mhdtime)
etime = max(mhdtime)

loc = where(dsttime ge stime and dsttime le etime)
mini = min([mhddst,mhdctc,dst(loc)])
maxi = max([mhddst,mhdctc,dst(loc)])

time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn


setdevice, 'dst_60.ps', 'p', 5, 0.95

ppp = 2
space = 0.05
pos_space, ppp, space, sizes, ny = ppp

makect, 'mid'

get_position, ppp, space, sizes, 0, pos, /rect
pos(0) = pos(0) + 0.05

plot, dsttime-stime, dst, yrange = [mini,maxi], $
  xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
  xminor = xminor, xticks = xtickn, xstyle = 1, $
  ytitle = 'Dst (nT)', pos = pos, thick = 3, charsize = 1.2, $
  psym = -4, xrange = [btr, etr]

oplot, mhdtime-stime, mhddst, linestyle = 1, thick = 3, color = 250

oplot, mhdtime-stime, mhdctc, color = 10, thick = 3, linestyle = 2

closedevice

end
