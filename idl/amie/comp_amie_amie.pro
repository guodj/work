
amie_file_1 = ask('1st AMIE file','b.x1.m.uvi.r1.121')
amie_file_2 = ask('2nd AMIE file','b.x1.m.uvi.r1.244')
psfile = ask('ps file',amie_file_1+'.comp.ps')

openr,1,amie_file_1
t = byte(0)
readu,1,t
close,1
if (t ge 32) then begin
  read_barbara, amie_file_1, amie_ltpos_1, amie_lnpos_1, 		$
                amie_data_1, amie_time_1, amie_date_1, 			$
	        amie_lats_1, amie_lons_1, amie_type_1, spd, by, bz
endif else begin
  print, '1. Potential'
  print, '2. Hall'
  print, '3. Pedersen'
  print, '4. Simple Joule Heating'
  type = fix(ask('field to compare','1'))
  if type eq 1 then begin
    field = 0
    cttitle = 'Potential (kV)'
  endif
  if type eq 2 then begin
    field = 4
    cttitle = 'Conductance (mhos)'
  endif
  if type eq 3 then begin
    field = 2
    cttitle = 'Conductance (mhos)'
  endif
  if type eq 4 then begin
    field = 8
    cttitle = 'Joule Heating (W/m!E2!N)'
  endif
  read_amie_binary, amie_file_1, amie_data_1, amie_lats_1, amie_lons_1,	$
                    amie_time_1, amie_type_1, imf,			$
                    ae, dst, hp, cpcp, date = amie_date_1,		$
                    ltpos = amie_ltpos_1, lnpos = amie_lnpos_1, 	$
                    /plotapot, field = field
endelse

;stop

openr,1,amie_file_2
t = byte(0)
readu,1,t
close,1
if (t ge 32) then begin
  read_barbara, amie_file_2, amie_ltpos_2, amie_lnpos_2, 		$
                amie_data_2, amie_time_2, amie_date_2, 			$
	        amie_lats_2, amie_lons_2, amie_type_2, spd, by, bz
endif else begin

  if strpos(amie_type_1,"Potential") gt -1 then begin
    field = 0
    if strpos(amie_file_1, amie_file_2) gt -1 then begin
      print, "Assuming that comparison between stat and derived is desired"
      field = 14
    endif
  endif
  if strpos(amie_type_1,"Hall") gt -1 then field = 4
  if strpos(amie_type_1,"Pedersen") gt -1 then field = 2
  if strpos(amie_type_1,"Joule") gt -1 then field = 8
  read_amie_binary, amie_file_2, amie_data_2, amie_lats_2, amie_lons_2,	$
                    amie_time_2, amie_type_2, imf,			$
                    ae, dst, hp, cpcp, date = amie_date_2,		$
                    ltpos = amie_ltpos_2, lnpos = amie_lnpos_2, 	$
                    /plotapot, field = field
endelse

;stop

n_amie_times_1 = n_elements(amie_data_1(*,0,0))
n_amie_times_2 = n_elements(amie_data_2(*,0,0))

amie_lnpos_1 = reform(amie_lnpos_1(*,0))*24.0/360.0
amie_ltpos_1 = reform(amie_ltpos_1(0,*))

amie_lnpos_2 = reform(amie_lnpos_2(*,0))*24.0/360.0
amie_ltpos_2 = reform(amie_ltpos_2(0,*))

min_lat = max([min(amie_ltpos_1),min(amie_ltpos_2)])
dlat = min([amie_ltpos_1(0)-amie_ltpos_1(1),amie_ltpos_2(0)-amie_ltpos_2(1)])
nlats = fix((90.0-min_lat)/dlat)+1
lats = 90.0-findgen(nlats)*dlat
nmlts = max([n_elements(amie_lnpos_1),n_elements(amie_lnpos_2)])
mlts = findgen(nmlts)*24.0/(nmlts-1)

lat2d = fltarr(nmlts,nlats)
lon2d = fltarr(nmlts,nlats)

for i=0,nlats-1 do lon2d(*,i) = mlts*!pi/12.0 - !pi/2.0
for i=0,nmlts-1 do lat2d(i,*) = lats

x = (90.0-lat2d)*cos(lon2d)
y = (90.0-lat2d)*sin(lon2d)

if (n_elements(amie_lats_2) ne n_elements(amie_lats_1)) or 		$
   (n_elements(amie_lons_2) ne n_elements(amie_lons_1)) then begin

  change_amie_resolution, amie_data_1, amie_ltpos_1, amie_lnpos_1, lats, mlts
  change_amie_resolution, amie_data_2, amie_ltpos_2, amie_lnpos_2, lats, mlts

endif

nl = n_elements(amie_ltpos_1)
nm = n_elements(amie_lnpos_1)

alat2d = fltarr(nm,nl)
alon2d = fltarr(nm,nl)

for i=0,nl-1 do alon2d(*,i) = amie_lnpos_1*!pi/12.0 - !pi/2.0
for i=0,nm-1 do alat2d(i,*) = amie_ltpos_1

ax = (90.0-alat2d)*cos(alon2d)
ay = (90.0-alat2d)*sin(alon2d)

maxi_1 = float(fix(max(abs(amie_data_1))/10.0))*10.0
maxi_2 = float(fix(max(abs(amie_data_2))/10.0))*10.0
maxi = max([maxi_1,maxi_2])
levels = findgen(21)*maxi/10.0 - maxi

ppp = 15
space = 0.005
mr = fix(max(90.0-lats)/10.0)*10.0 + 10.0

setdevice, psfile,'p'

readct, ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
clevels = (ncolors-1)*findgen(21)/20.0 + 1

pos_space, ppp, space, sizes, nx = 3

pn = -3

amie_pattern_1 = fltarr(nmlts,nlats)
amie_pattern_2 = fltarr(nmlts,nlats)
n_pattern      = intarr(nmlts,nlats)
time_stats     = fltarr(4,n_amie_times_1)-9999.0
time_save      = dblarr(n_amie_times_1)

ut = dblarr(n_amie_times_2)

for i=0, n_amie_times_2-1 do begin

  stime = strmid(amie_date_2(i),4,2)+'-'+strmid(amie_date_2(i),0,3)+'-'+$
          strmid(amie_date_2(i),10,2)
  stime = stime + ' ' + strmid(amie_time_2(i),0,5)
  c_s_to_a,itime,stime
  c_a_to_r,itime,rtime
  ut(i) = rtime

endfor

for n = 0, n_amie_times_1-1 do begin

  stime = strmid(amie_date_1(n),4,2)+'-'+strmid(amie_date_1(n),0,3)+'-'+$
          strmid(amie_date_1(n),10,2)
  stime = stime + ' ' + strmid(amie_time_1(n),0,5)
  c_s_to_a,itime,stime
  c_a_to_r,itime,rtime

  time_save(n) = rtime
  time_stats(0,n) = max(amie_data_1(n,*,*))-min(amie_data_1(n,*,*))

  d = abs(ut - rtime)
  loc = where(d eq min(d), count)

  namie2 = loc(0)
  nminutes = d(namie2)/60.0

  pn = (pn + 3) mod ppp
  if pn eq 0 then begin
    plotdumb
    xyouts, 0.0, 1.04, amie_type_1,/norm
    xyouts, 0.0, 1.01, amie_date_1(n),/norm
  endif

  get_position, ppp, space, sizes, pn, pos
  ;print, namie2
  ;help, x
  ;help, y
  contour, reform(amie_data_2(namie2,*,*)),x,y,/follow, 		$
	xstyle = 5, ystyle = 5,$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase, /cell_fill, c_color = clevels

  xyouts, pos(0)-0.01,pos(3)-0.01, amie_time_1(n),charsize=0.9, /norm

  if (pn+3 ge ppp) or (n eq n_amie_times_1-1) then 			$
    plotmlt, mr, /no06, /no12 else					$
    if pn-3 lt 0 then begin
      plotmlt, mr, /no06, /no00 
      xyouts, mean(pos([0,2])), 1.01, amie_file_2, /norm, align=0.5
    endif else begin
      plotmlt, mr, /no06, /no00, /no12
    endelse

  xyouts, pos(0)-space*5,mean(pos([1,3])),				$
	'!9D!XT:'+tostr(fix(nminutes))+" min.",	$
	align = 0.5, orient = 90, /norm

  get_position, ppp, space, sizes, pn+1, pos
  contour, reform(amie_data_1(n,*,*)),x,y,/follow, xstyle = 5, ystyle = 5,$
         xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, /noerase,	$
	pos = pos, /cell_fill, c_color = clevels

  if pn+3 ge ppp then plotmlt, mr, /no06, /no12, /no18 else		$
    if pn-3 lt 0 then begin
      plotmlt, mr, /no06, /no00, /no18 
      xyouts, mean(pos([0,2])), 1.01, amie_file_1, /norm, align=0.5
    endif else begin
      plotmlt, mr, /no06, /no00, /no12, /no18
    endelse

  diff = reform(amie_data_1(n,*,*)-amie_data_2(namie2,*,*))

  amie_pattern_1(*,*) = amie_pattern_1(*,*) + amie_data_1(n,*,*)
  amie_pattern_2(*,*) = amie_pattern_2(*,*) + amie_data_2(namie2,*,*)
  n_pattern           = n_pattern + 1
  time_stats(1,n) = max(amie_data_2(namie2,*,*))-min(amie_data_2(namie2,*,*))
  time_stats(2,n) = abs(time_stats(1,n) - time_stats(0,n))
  time_stats(3,n) = max(abs(diff))

  get_position, ppp, space, sizes, pn+2, pos
  contour, diff,x,y,/follow, xstyle = 5, ystyle = 5,$
         xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase, /cell_fill, c_color = clevels

  if pn+3 ge ppp then begin
    plotmlt, mr, /no12, /no18,/no06
    ctpos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
    plotct, ncolors, ctpos, [-maxi,maxi], cttitle, /right
  endif else if pn-3 lt 0 then begin
    plotmlt, mr, /no00, /no18 
    xyouts, mean(pos([0,2])), 1.01, "Diff.", /norm, align=0.5
  endif else begin
    plotmlt, mr, /no00, /no12, /no18
  endelse

endfor

closedevice

loc = where(n_pattern gt 0,count)
if count gt 0 then begin
  amie_pattern_1(loc) = amie_pattern_1(loc)/float(n_pattern(loc))
  amie_pattern_2(loc) = amie_pattern_2(loc)/float(n_pattern(loc))
endif

setdevice, strmid(psfile,0,strpos(psfile,'.ps'))+'.diff.ps','p'

dy = float(!d.y_ch_size)/float(!d.y_size)

plotdumb

xyouts, 1.0, 1.01, strcompress(amie_type_1),/norm, align=1.0

maxi = float(fix(max([max(abs(amie_pattern_1)),				$
                      max(abs(amie_pattern_2))])/5.0)+1)*5.0
levels = findgen(21)*maxi/10.0 - maxi

pn = 0
get_position, ppp, space, sizes, pn, pos
pos([1,3]) = pos([1,3]) - 0.02
x1 = mean(pos([0,2]))
contour, amie_pattern_2,x,y,/follow, xstyle = 5, ystyle = 5,$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase, /cell_fill, c_color = clevels
contour, amie_pattern_2,x,y,/follow, xstyle = 5, ystyle = 5,$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase
plotmlt, mr, /no06
xyouts, mean(pos([0,2])), pos(3)+1.5*dy, amie_file_2, /norm, alignment = 0.5

pn = 1
get_position, ppp, space, sizes, pn, pos
pos([1,3]) = pos([1,3]) - 0.02
x2 = mean(pos([0,2]))
y2 = mean(pos([1,3]))
size_y = pos(3) - pos(1)
size_x = pos(2) - pos(0)
contour, amie_pattern_1,x,y,/follow, xstyle = 5, ystyle = 5,$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase, /cell_fill, c_color = clevels
contour, amie_pattern_1,x,y,/follow, xstyle = 5, ystyle = 5,$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase
plotmlt, mr, /no06, /no18
xyouts, mean(pos([0,2])), pos(3)+1.5*dy, amie_file_1, /norm, alignment = 0.5

ctpos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
plotct, ncolors, ctpos, [-maxi,maxi], cttitle, /right

between  = (x2-x1)/2.0
x_center = x2 - between
y_center = y2 - sqrt(3)*((size_y+space)/2.0)

pos(0) = x_center - size_x/2
pos(2) = pos(0) + size_x
pos(1) = y_center - size_y/2
pos(3) = pos(1) + size_y

contour, amie_pattern_1-amie_pattern_2,x,y,/follow, 		$
	xstyle = 5, ystyle = 5,						$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase, /cell_fill, c_color = clevels
contour, amie_pattern_1-amie_pattern_2,x,y,/follow, 		$
	xstyle = 5, ystyle = 5,						$
	xrange = [-mr,mr],yrange=[-mr,mr], levels = levels, 		$
	pos = pos, /noerase
plotmlt, mr,/no00

pos(3) = pos(1) - space
pos(1) = pos(3) - size_y*0.3
pos(0) = pos(0) - size_x/4.0
pos(2) = ctpos(2) + space*10.0

stime = min(time_save)
etime = max(time_save)

time_axis, stime, etime, s_time_range, e_time_range,        		$
	xtickname, xtitle, xtickvalue, xminor, xtickn

plot, time_save-stime,imf(*,2), pos = pos,	 			$
	xstyle = 1, ytitle = 'Bz (nT)', /noerase,			$
	xtickname = strarr(10)+' ', xtickv=xtickvalue, 			$
	xticks = xtickn, xminor = xminor, xtitle = ' ',		$
	xrange = [s_time_range, e_time_range], 				$
	ytickname = [' ','',' ','',' ','',' ','',' ','',' ']

oplot, mm(time_save-stime),[0.0,0.0], linestyle = 1

pos(3) = pos(1) - space
pos(1) = pos(3) - size_y*0.6

maxi = float(fix(max([time_stats(0,*),time_stats(1,*)])/10.0)+1)*10.0

plot, time_save-stime,time_stats(0,*), pos = pos, 			$
	yrange = [0,maxi], xstyle = 1, ystyle = 1, 			$
	ytitle = cttitle, /noerase, min_value = -maxi,			$
	xtickname = xtickname, xtickv=xtickvalue, 			$
	xticks = xtickn, xminor = xminor, xtitle = xtitle,		$
	xrange = [s_time_range, e_time_range]

loc = where(time_stats(2,*) gt -maxi,count)
if count gt 0 then begin
  xyouts, pos(0), pos(1) - 0.08, 'Average Deviation : '+	$
	string(mean(time_stats(2,loc))), /norm
  xyouts, pos(0), pos(1) - 0.10, 'Average Maximum Deviation : '+	$
	string(mean(time_stats(3,loc))), /norm
endif

oplot, time_save-stime, time_stats(1,*), linestyle = 2, min_value = -maxi
oplot, time_save-stime, time_stats(2,*), linestyle = 1, min_value = -maxi

closedevice

end
