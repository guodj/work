
if (n_elements(file) eq 0) then begin

  list = findfile("*.save")
  file = list(0)
  if (strlen(file) lt 1) then file = "unknown.save"

endif

file = ask("save file to plot",file)

restore, file

icpcp = n_elements(indices(0,0,*,0))-1
cpcp_n = reform(indices(0,*,icpcp,1) - indices(0,*,icpcp,0))
cpcp_s = reform(indices(1,*,icpcp,1) - indices(1,*,icpcp,0))

;time = time-time(0)

stime = min(time)
etime = max(time)

if (stime lt 3600.0) then begin
    if (max(time) gt 7200.0) then begin
        time = time/3600.0
        xtitle = 'Simulation Time (Hours)'
    endif else begin
        if (max(time) gt 60.0) then begin
            time = time/60.0
            xtitle = 'Simulation Time (Minutes)'
        endif else begin

            xtitle = 'Simulation Time (Seconds)'

        endelse
    endelse
endif else begin

    c_r_to_a, itime, stime
    itime(4) = 0
    itime(5) = 0
    c_a_to_r, itime, stime

    c_r_to_a, itime, etime
    itime(3) = itime(3) + 1
    itime(4) = 0
    itime(5) = 0
    c_a_to_r, itime, etime

    time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn
    
endelse

time = time-stime

if (stime-etime lt 1.0) then begin

    nTimes = n_elements(time)
    time = findgen(nTimes)
    xtitle = 'Iterations'
    time = time-min(time)
    stime = 0.0

endif


ytitle = 'Potential (kV)'

setdevice, 'cpcp_north.ps', 'p', 5, 0.95

yrange = mm([cpcp_n, cpcp_s])
;yrange = [0.0,150.0]


ppp = 2
space = 0.05
pos_space, ppp, space, sizes, ny = ppp

get_position, ppp, space, sizes, 0, pos, /rect
pos(0) = pos(0) + 0.05

xrange = mm(time)
;xrange = [0.0,120.0]

if (stime lt 3600.0) then begin
    plot, time, cpcp_n, xstyle = 1, xtitle = xtitle, ytitle = ytitle, $
      charsize = 1.2, thick = 3, pos = pos, yrange = yrange, $
      title = 'Northern Hemisphere Cross Polar Cap Potential', $
      xrange = xrange
endif else begin


print, xtickn

    plot, time, cpcp_n, ytitle = ytitle, $
      charsize = 1.2, thick = 3, pos = pos, yrange = yrange, $
      title = 'Northern Hemisphere Cross Polar Cap Potential', $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1

endelse


;loc = where(cpcp_n lt 20.0,count)
;
;if count gt 0 then begin
;
;    m = median(cpcp_n(loc))
;;    l = where(cpcp_n eq m)
;;    oplot, [time(l),time(l)], [0.0,1000]
;    d = (20.0 + m)/2.0
;
;    l = loc(count-1)
;    i = l
;    while (i gt 0 and cpcp_n(i) gt d) do i = i - 1
;    if (i gt 1) then l = i
;    oplot, [time(l),time(l)], [0.0,1000]
;    t1 = time(l)
;    t = 't ='+string(time(l), format='(f5.1)')+'m'
;    xyouts, time(l)-1, 110.0, t, alignment = 1.0
;
;
;    while (i lt n_elements(cpcp_n)-1 and cpcp_n(i+1) gt cpcp_n(i)) do i = i + 1
;    if (i lt n_elements(cpcp_n)-1) then l = i
;    oplot, [time(l),time(l)], [0.0,1000]
;    t2 = time(l)
;    t = 't ='+string(time(l), format='(f5.1)')+'m'
;    xyouts, time(l)+1, 110.0, t
;
;    oplot, [t1,t2], [10.0,10.0], linestyle = 2
;    td = 'dt='+string(t2-t1, format='(f5.1)')+'m'
;    xyouts, (t1+t2)/2.0, 12.0, td, alignment = 0.5, charsize = 0.8
;
;
;    m = cpcp_n(l)
;    ms = 'Max = '+string(m, format='(f6.1)')+' kV'
;
;    oplot, [0.0,120.0],[m,m], linestyle = 1
;    xyouts, 100.0, m+2, ms, alignment = 0.5
;
;endif

spawn, 'pwd',pwd

xyouts, pos(0), pos(1)-0.1, pwd, /norm

closedevice 

setdevice, 'cpcp_south.ps', 'p', 5, 0.95

if (stime lt 3600.0) then begin

    plot, time, cpcp_s, xstyle = 1, xtitle = xtitle, ytitle = ytitle, $
      charsize = 1.2, thick = 3, pos = pos, yrange = yrange, $
      title = 'Southern Hemisphere Cross Polar Cap Potential'

endif else begin

    plot, time, cpcp_s, ytitle = ytitle, $
      charsize = 1.2, thick = 3, pos = pos, yrange = yrange, $
      title = 'Southern Hemisphere Cross Polar Cap Potential', $
      xtickname = xtickname, xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, xstyle = 1

endelse

closedevice 

end

