
;---------------------------------------------------------------------------
; Have to change these values:

logfilename = 'log_n005001.log'

logfilename = ask('logfilename',logfilename)

iYear  = fix(ask('start year of run','2003'))
iMonth = fix(ask('start Month of run','10'))
iDay   = fix(ask('start Day of run','29'))
iHour  = fix(ask('start Hour of run','03'))

itime = [iYear, iMonth, iDay, iHour, 00, 00]

dstfile = '/r/ridley/MHD/20031029/Data/dst200310.html'
dstfile = ask('dstfile',dstfile)

;---------------------------------------------------------------------------

getlogpro, logfilename, nlogfile, nVars, wlognames, wlog, wlog0, wlog1, wlog2
nPts = n_elements(wlog(*,0))

c_a_to_r, itime, basetime
mhdtime = dblarr(nPts) + basetime
mhdtime = mhdtime + wlog(*,1)
mhddst = reform(wlog(*,15))
mhdctc = reform(wlog(*,17))

read_dst, dstfile, dst, dsttime

stime = min(mhdtime)
etime = max(mhdtime)

endtime = tostrf((etime-stime)/3600.0)
endtime = float(ask('end time to plot (hours)',endtime))
etime   = stime + endtime * 3600.0

loc = where(dsttime ge stime and dsttime le etime)
mini = min([mhddst,mhdctc,dst(loc)])
maxi = max([mhddst,mhdctc,dst(loc)])

time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn


setdevice, 'dst.ps', 'p', 5, 0.95

makect, "mid"

ppp = 2
space = 0.05
pos_space, ppp, space, sizes, ny = ppp

get_position, ppp, space, sizes, 0, pos, /rect
pos(0) = pos(0) + 0.05

plot, mhdtime-stime, mhddst, yrange = [mini,maxi], $
  xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
  xminor = xminor, xticks = xtickn, xstyle = 1, $
  ytitle = 'Dst (nT)', pos = pos, thick = 3, charsize = 1.2, $
  xrange = [btr, etr]

oplot, mhdtime-stime, mhdctc, linestyle = 1, color = 250

oplot, dsttime-stime, dst, linestyle = 2, psym = -4, color = 10

closedevice

end
