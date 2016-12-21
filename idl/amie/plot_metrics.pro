
if (n_elements(amie_file) eq 0) then begin

  initial_guess = findfile('-t b*_metrics')
  initial_guess = initial_guess(0)
  if strlen(initial_guess) eq 0 then initial_guess='b970101_metrics'

endif else initial_guess = amie_file

amie_file = ask('AMIE metrics file name',initial_guess)
psfile = ask('ps file',amie_file+'.ps')

setdevice,psfile,'p',4,0.9

threshold = float(ask('threshold for cc w/AVE cut off','0.5'))

openr, 2, amie_file

readf,2,nFiles
readf,2,nMetrics

Metrics = fltarr(nMetrics,2,nFiles)
MetricsTimes = dblarr(nFiles)
MetricTitle = strarr(nMetrics)
line = ''
for i=0,nMetrics-1 do begin
  readf,2,line
  MetricTitle(i) = line
endfor
readf,2,metrics
readf,2,MetricsTimes

close,2

loc = where(MetricsTimes gt 0.0)

stime = min(MetricsTimes(loc))
etime = max(MetricsTimes(loc))

c_r_to_a, itime, stime
itime([3,4,5]) = 0
c_a_to_r, itime, stime

c_r_to_a, eitime, etime
if (eitime(2) eq itime(2)) then begin
  itime([3,4,5]) = [23,59,59]
  c_a_to_r, itime, etime
endif

time_axis, stime, etime, s_time_range, e_time_range,        	$
	xtickname, xtitle, xtickvalue, xminor, xtickn

plotdumb

ppp = nMetrics
space = 0.005
pos_space, ppp, space, sizes, nx = 1

for iMetric = 0,nMetrics-1 do begin

  get_position, ppp, space, sizes, iMetric, pos, /rect
  pos(0) = pos(0) + 0.05

  tmp = reform(metrics(iMetric,*,*))
  loc = where(tmp gt -100.0)
  yrange = mm(tmp(loc))

  if (iMetric eq 0) then yrange = [-0.5,1.0]

  if iMetric eq nMetrics-1 then begin
    xtn = xtickname
    xt = xtitle
  endif else begin
    xtn = strarr(20)+' '
    xt = ' '
  endelse

  loc = where(metrics(3,0,*) ge threshold and $
              metrics(iMetric,0,*) gt -999.0 and $
              metrics(5,0,*) ge 5.0)

  if (loc(0) gt -1) then begin
    plot, MetricsTimes(loc)-stime,metrics(iMetric,0,loc), pos = pos,	 $
	xstyle = 1, /noerase,		$
	xtickname = xtn, xtickv=xtickvalue, 			$
	xticks = xtickn, xminor = xminor, xtitle = xt,		$
	xrange = [s_time_range, e_time_range], psym = 4, $
        ytitle = MetricTitle(iMetric), yrange = yrange

    oplot, [s_time_range, e_time_range],$
	fltarr(2)+mean(metrics(iMetric,0,loc)), linestyle = 0
  endif

  if iMetric eq 0 then begin
    loc = where(metrics(3,0,*) ge threshold and $
                metrics(iMetric,0,*) gt -999.0 and $
		metrics(iMetric,0,*) lt -0.5 and $
		metrics(5,0,*) ge 5.0)
    if (loc(0) gt -1) then $
      oplot, MetricsTimes(loc)-stime,metrics(iMetric,0,loc)*0.0-0.5, psym = 5
  endif

  ;------------------------------------------------

  loc = where(metrics(3,1,*) ge threshold and $
              metrics(iMetric,1,*) gt -999.0 and $
              metrics(5,1,*) ge 5.0)

  if (loc(0) gt -1) then begin
    oplot, MetricsTimes(loc)-stime,metrics(iMetric,1,loc), psym = 2
    oplot, [s_time_range, e_time_range],[0,0], linestyle = 1
    oplot, [s_time_range, e_time_range],$
	  fltarr(2)+mean(metrics(iMetric,1,loc)), linestyle = 2
  endif

  if iMetric eq 0 then begin
    loc = where(metrics(3,1,*) ge threshold and $
                metrics(iMetric,1,*) gt -999.0 and $
		metrics(iMetric,1,*) lt -0.5 and $
		metrics(5,1,*) ge 5.0)
    if (loc(0) gt -1) then $
      oplot, MetricsTimes(loc)-stime,metrics(iMetric,1,loc)*0.0-0.5, psym = 5
  endif

endfor

xyouts, 0.0, -0.05, amie_file, /norm, charsize = 0.75

plotdumb

ppp = nMetrics
space = 0.005
pos_space, ppp, space, sizes, nx = 1

for iMetric = 0,nMetrics-1 do begin

  get_position, ppp, space, sizes, iMetric, pos, /rect
  pos(0) = pos(0) + 0.05

  tmp = reform(metrics(iMetric,*,*))
  loc = where(tmp gt -100.0)
  yrange = mm(tmp(loc))

  if (iMetric eq 0) then yrange = [-0.5,1.0]

  if iMetric eq nMetrics-1 then begin
    xtn = xtickname
    xt = xtitle
  endif else begin
    xtn = strarr(20)+' '
    xt = ' '
  endelse

  loc = where(metrics(3,0,*) ge threshold and $
              metrics(iMetric,0,*) gt -999.0 and $
              metrics(5,0,*) lt 5.0)

  plot, MetricsTimes(loc)-stime,metrics(iMetric,0,loc), pos = pos,	 $
	xstyle = 1, /noerase,		$
	xtickname = xtn, xtickv=xtickvalue, 			$
	xticks = xtickn, xminor = xminor, xtitle = xt,		$
	xrange = [s_time_range, e_time_range], psym = 4, $
        ytitle = MetricTitle(iMetric), yrange = yrange

  oplot, [s_time_range, e_time_range],$
	fltarr(2)+mean(metrics(iMetric,0,loc)), linestyle = 0

  if iMetric eq 0 then begin
    loc = where(metrics(3,0,*) ge threshold and $
                metrics(iMetric,0,*) gt -999.0 and $
		metrics(iMetric,0,*) lt -0.5 and $
		metrics(5,0,*) lt 5.0)
    if (loc(0) gt -1) then $
      oplot, MetricsTimes(loc)-stime,metrics(iMetric,0,loc)*0.0-0.5, psym = 5
  endif

  ;------------------------------------------------

  loc = where(metrics(3,1,*) ge threshold and $
              metrics(iMetric,1,*) gt -999.0 and $
              metrics(5,1,*) lt 5.0)

  oplot, MetricsTimes(loc)-stime,metrics(iMetric,1,loc), psym = 2
  oplot, [s_time_range, e_time_range],[0,0], linestyle = 1
  oplot, [s_time_range, e_time_range],$
	fltarr(2)+mean(metrics(iMetric,1,loc)), linestyle = 2

  if iMetric eq 0 then begin
    loc = where(metrics(3,1,*) ge threshold and $
                metrics(iMetric,1,*) gt -999.0 and $
		metrics(iMetric,1,*) lt -0.5 and $
		metrics(5,1,*) lt 5.0)
    if (loc(0) gt -1) then $
      oplot, MetricsTimes(loc)-stime,metrics(iMetric,1,loc)*0.0-0.5, psym = 5
  endif

endfor

xyouts, 0.0, -0.05, amie_file, /norm, charsize = 0.75


end
