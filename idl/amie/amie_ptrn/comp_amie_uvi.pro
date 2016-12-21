
amie_file = ask('AMIE file','amie.0600.244.hal')
UVI_file = ask('UVI file','d.data/UVI_124cond.gen')
psfile = ask('ps file',amie_file+'.uvi_comp.ps')
mspall = float(ask('MSPALL value in amiein file','3.0'))

read_UVI, UVI_file, ut, hall, ped, lats, mlts, /indicate

nlats = n_elements(lats)
nmlts = n_elements(mlts)

lat2d = fltarr(nmlts,nlats)
lon2d = fltarr(nmlts,nlats)

for i=0,nlats-1 do lon2d(*,i) = mlts*!pi/12.0 - !pi/2.0
for i=0,nmlts-1 do lat2d(i,*) = lats

x = (90.0-lat2d)*cos(lon2d)
y = (90.0-lat2d)*sin(lon2d)

openr,1,amie_file
t = byte(0)
readu,1,t
close,1
if (t ge 32) then begin
  read_barbara, amie_file, amie_ltpos, amie_lnpos, 			$
                amie_data, amie_time, amie_date, 			$
	        amie_lats, amie_lons, amie_type, spd, by, bz
endif else begin
  print, '1. Hall'
  print, '2. Pedersen'
  type = fix(ask('Conductance to plot','1'))
  if type eq 1 then field = 4 else field = 2
  read_amie_binary, amie_file, amie_data, amie_lats, amie_lons,		$
                    amie_time, amie_type, imf,				$
                    ae, dst, hp, cpcp, date = amie_date,		$
                    ltpos = amie_ltpos, lnpos = amie_lnpos, 		$
                    /plotapot, field = field
endelse

if (strpos(mklower(amie_type),"hall") gt -1) then UVI_data = hall 	$
else UVI_data = ped

print, "Done reading AMIE file"

n_amie_times = n_elements(amie_data(*,0,0))

amie_lnpos = reform(amie_lnpos(*,0))*24.0/360.0
amie_ltpos = reform(amie_ltpos(0,*))

change_amie_resolution, amie_data, amie_ltpos, amie_lnpos, lats, mlts

nl = n_elements(amie_ltpos)
nm = n_elements(amie_lnpos)

alat2d = fltarr(nm,nl)
alon2d = fltarr(nm,nl)

for i=0,nl-1 do alon2d(*,i) = amie_lnpos*!pi/12.0 - !pi/2.0
for i=0,nm-1 do alat2d(i,*) = amie_ltpos

ax = (90.0-alat2d)*cos(alon2d)
ay = (90.0-alat2d)*sin(alon2d)

maxi = float(fix(max(UVI_data)/10.0))*10.0
levels = findgen(21)*maxi/10.0 - maxi

ppp = 15
space = 0.005
mr = fix(max(90.0-lats)/10.0)*10.0 + 10.0

setdevice, psfile,'p'

readct, ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
clevels = (ncolors-1)*findgen(21)/20.0 + 1

pos_space, ppp, space, sizes, nx = 3

pn = -3

amie_pattern = fltarr(nmlts,nlats)
UVI_pattern  = fltarr(nmlts,nlats)
n_pattern    = intarr(nmlts,nlats)
time_stats   = fltarr(4,n_amie_times)-9999.0
time_save    = dblarr(n_amie_times)

for n = 0, n_amie_times-1 do begin

  stime = strmid(amie_date(n),4,2)+'-'+strmid(amie_date(n),0,3)+'-'+$
          strmid(amie_date(n),10,2)
  stime = stime + ' ' + strmid(amie_time(n),0,5)
  c_s_to_a,itime,stime
  c_a_to_r,itime,rtime

  time_save(n) = rtime
  time_stats(0,n) = max(amie_data(n,*,*))

; The way AMIE works is that it will take the first time within the mspall window.
; This needs to be reflected in the IDL plotting.

  d = ut - (rtime-mspall)
  loc = where(d gt 0, count)

  if count gt 0 then begin
    nUVI = loc(0)
    nminutes = d(nUVI)/60.0
  endif else nminutes = mspall*2.0

  if (nminutes le mspall) then begin

    pn = (pn + 3) mod ppp
    if pn eq 0 then begin
      plotdumb
      xyouts, 0.0, 1.04, amie_type,/norm
      xyouts, 0.0, 1.01, amie_date(n),/norm
    endif

    get_position, ppp, space, sizes, pn, pos

    contour, reform(UVI_data(nUVI,*,*)),x,y,/follow, xstyle = 5, ystyle = 5,$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	min_value=0.001, pos = pos, /noerase, /cell_fill, c_color = clevels

    xyouts, pos(0)-0.01,pos(3)-0.01, amie_time(n),charsize=0.9, /norm

    if (pn+3 ge ppp) or (n eq n_amie_times-1) then 			$
      plotmlt, mr, /no06, /no12 else					$
      if pn-3 lt 0 then begin
        plotmlt, mr, /no06, /no00 
	xyouts, mean(pos([0,2])), 1.01, "UVI", /norm, align=0.5
      endif else begin
        plotmlt, mr, /no06, /no00, /no12
      endelse

    xyouts, pos(0)-space*5,mean(pos([1,3])),				$
	'!9D!XT:'+tostr(fix(nminutes))+" min.",	$
	align = 0.5, orient = 90, /norm

    get_position, ppp, space, sizes, pn+1, pos
    contour, reform(amie_data(n,*,*)),x,y,/follow, xstyle = 5, ystyle = 5,$
         xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, /noerase,	$
	pos = pos, /cell_fill, c_color = clevels

    if pn+3 ge ppp then plotmlt, mr, /no06, /no12, /no18 else		$
      if pn-3 lt 0 then begin
        plotmlt, mr, /no06, /no00, /no18 
	xyouts, mean(pos([0,2])), 1.01, amie_file, /norm, align=0.5
      endif else begin
        plotmlt, mr, /no06, /no00, /no12, /no18
      endelse

    diff = reform(amie_data(n,*,*)-UVI_data(nUVI,*,*))

    loc = where(reform(UVI_data(nUVI,*,*)) eq 0, count)
    if count gt 0 then diff(loc) = 0.0

    UVI_tmp = reform(UVI_data(nUVI,*,*))
    amie_tmp = reform(amie_data(n,*,*))
    loc = where(UVI_tmp gt 0, count)
    if count gt 0 then begin
      amie_pattern(loc) = amie_pattern(loc) + amie_tmp(loc)
      UVI_pattern(loc)  = UVI_pattern(loc) + UVI_tmp(loc)
      n_pattern(loc)    = n_pattern(loc) + 1
      time_stats(1,n) = max(UVI_data(nUVI,*,*))
      time_stats(2,n) = min(diff(loc))
      time_stats(3,n) = max(diff(loc))
    endif

    get_position, ppp, space, sizes, pn+2, pos
    contour, diff,x,y,/follow, xstyle = 5, ystyle = 5,$
         xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase, /cell_fill, c_color = clevels

    if pn+3 ge ppp then begin
      plotmlt, mr, /no12, /no18,/no06
      ctpos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
      plotct, ncolors, ctpos, [-maxi,maxi], "(mhos)", /right
    endif else if pn-3 lt 0 then begin
      plotmlt, mr, /no00, /no18 
      xyouts, mean(pos([0,2])), 1.01, "Diff.", /norm, align=0.5
    endif else begin
      plotmlt, mr, /no00, /no12, /no18
    endelse

  endif

endfor

closedevice

loc = where(n_pattern gt 0,count)
if count gt 0 then begin
  amie_pattern(loc) = amie_pattern(loc)/float(n_pattern(loc))
  UVI_pattern(loc) = UVI_pattern(loc)/float(n_pattern(loc))
endif

setdevice, strmid(psfile,0,strpos(psfile,'.ps'))+'.diff.ps','p'

plotdumb

xyouts, 0.5, 1.01, amie_type,/norm, align=0.5

pn = 3
get_position, ppp, space, sizes, pn, pos_temp
pos([0,1]) = pos_temp([0,1])
pn = 5
get_position, ppp, space, sizes, pn, pos_temp
pos([2,3]) = pos_temp([2,3])
pos([1,3]) = pos([1,3]) - 0.06

stime = min(time_save)
etime = max(time_save)

time_axis, stime, etime, s_time_range, e_time_range,        		$
	xtickname, xtitle, xtickvalue, xminor, xtickn

plot, time_save-stime,time_stats(0,*), pos = pos, 			$
	yrange = [-maxi,maxi], xstyle = 1, ystyle = 1, 			$
	ytitle = 'mhos', /noerase, min_value = -maxi,			$
	xtickname = xtickname, xtickv=xtickvalue, 			$
	xticks = xtickn, xminor = xminor, xtitle = xtitle,		$
	xrange = [s_time_range, e_time_range]

loc = where(time_stats(2,*) gt -maxi,count)
if count gt 0 then begin
  xyouts, pos(0), pos(1) - 0.08, 'Average Maximum Deviation : '+		$
	string(mean(time_stats(2,loc))), /norm
  xyouts, pos(0), pos(1) - 0.10, 'Average Minimum Deviation : '+		$
	string(mean(time_stats(3,loc))), /norm
endif

oplot, time_save-stime, time_stats(1,*), linestyle = 2, min_value = -maxi
oplot, time_save-stime, time_stats(2,*), linestyle = 1, min_value = -maxi
oplot, time_save-stime, time_stats(3,*), linestyle = 1, min_value = -maxi

maxi = maxi/1.5
levels = findgen(21)*maxi/10.0 - maxi

pn = 0
get_position, ppp, space, sizes, pn, pos
pos([1,3]) = pos([1,3]) - 0.02
contour, UVI_pattern,x,y,/follow, xstyle = 5, ystyle = 5,$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	min_value=0.001, pos = pos, /noerase, /cell_fill, c_color = clevels
plotmlt, mr, /no06

pn = 1
get_position, ppp, space, sizes, pn, pos
pos([1,3]) = pos([1,3]) - 0.02
contour, amie_pattern,x,y,/follow, xstyle = 5, ystyle = 5,$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	min_value=0.001, pos = pos, /noerase, /cell_fill, c_color = clevels
plotmlt, mr, /no06, /no18

pn = 2
get_position, ppp, space, sizes, pn, pos
pos([1,3]) = pos([1,3]) - 0.02
contour, amie_pattern-UVI_pattern,x,y,/follow, 		$
	xstyle = 5, ystyle = 5,						$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase, /cell_fill, c_color = clevels
plotmlt, mr, /no18, /no06

ctpos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
plotct, ncolors, ctpos, [-maxi,maxi], "(mhos)", /right

closedevice




end