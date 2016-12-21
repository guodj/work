
common ffinfo, header

maxnumfil = 2
exten = 'hed cdf'
initpath = '/home/hao/ridleya/d.c/d.flatfiles/'
;get_files, datafile, numoffil, maxnumfil, filtitle, initpath, exten
numoffil = 1
datafile = strarr(6)
datafile(0) = initpath+'test'

setup_var_menu, numoffil, datafile, timeint, titles,            $
	units, varcount, ntotrows, rowlen
header = {nf : numoffil, df : datafile, ti : timeint,           $
	na : titles, un : units, nr : ntotrows, rl : rowlen,          $
	unit : intarr(numoffil)}

c_s_to_a, stime, timeint(0,0)
c_s_to_a, etime, timeint(0,1)

print, stime
print, etime

c_a_to_r, stime, rstime
c_a_to_r, etime, retime

print, rstime
print, retime

unit = 0
header.unit(unit) = 5
openr,header.unit(unit),datafile(unit)+'.dat'

print, 'first time : ',loc_ffunix_row(unit,rstime)
print, 'last time : ',loc_ffunix_row(unit,retime)
print, 'midtime time : ',loc_ffunix_row(unit,(rstime+retime)/2.0)
print, 'midtime time + 10 min : ',loc_ffunix_row(unit,(rstime+retime)/2.0+600)

close,header.unit(unit)

numofvar = 1
grafvar = intarr(2,4)
grafvar(0,0) = 0
grafvar(0,1) = 1
grafvar(1,0) = 0
grafvar(1,1) = 2
plottype = 1
times = intarr(6,2,6)
times(0,0,*) = stime
times(0,1,*) = etime
times(1,0,*) = stime
times(1,1,*) = etime

get_ffunix_data, single, t, dt, basetime, numofvar, numofpts, 	$
	plottype, times, datafile, grafvar, same

help, single

end