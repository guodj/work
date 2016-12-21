
;--------------------------------------------------------------
; Get Inputs from the user
;--------------------------------------------------------------

initial_guess = findfile('-t b*')
initial_guess = initial_guess(0)
if strlen(initial_guess) eq 0 then initial_guess='b970101'

amie_file = ask('AMIE binary file name',initial_guess)
psfile = ask('ps file',amie_file+'.ps')

ppp = 6

read_amie_binary, amie_file, data, lats, mlts, time, fields, 		$
                  imf, ae, dst, hp, cpcp

;--------------------------------------------------------------
; figure out grid:
;--------------------------------------------------------------

nlats = n_elements(lats)
nmlts = n_elements(mlts)
ntimes = n_elements(time)

;--------------------------------------------------------------
; Need to figure out what to plot, so list to fields to the user
;--------------------------------------------------------------

nfields = n_elements(fields)
;for i=0,nfields-1 do print, tostr(i+1)+'. '+fields(i)

;--------------------------------------------------------------
; Get field to be contoured
;--------------------------------------------------------------

;type_1 = fix(ask('field to plot','1'))-1
;if (type_1 lt 0) or (type_1 gt nfields-1) then type_1 = 0

type_1 = 6

;--------------------------------------------------------------
; Get start time and end time, with defaults as the file times
;--------------------------------------------------------------

c_r_to_a, itime_start, min(time)
c_r_to_a, itime_end, max(time)

c_a_to_s, itime_start, stime_start
c_a_to_s, itime_end, stime_end

start_time = ask('start time of plotting',strmid(stime_start,0,15))

if (strlen(start_time) lt 9) then $
  start_time = strmid(stime_start,0,9)+' '+start_time

c_s_to_a, itime_start, start_time
c_a_to_r, itime_start, rtime_start

;--------------------------------------------------------------
; I got sick of typing in the ending date, so if the date is
; the same, assume the user just wants to enter the time
;--------------------------------------------------------------

sdate = strmid(start_time,0,9)
if strpos(stime_end,sdate) gt -1 then 					$
  end_time_default = strmid(stime_end,10,5)				$
else end_time_default = strmid(stime_end,0,15)

end_time   = ask('end time of plotting',end_time_default)

;--------------------------------------------------------------
; If the user entered a short string, assume it is just a time
; and add the date on the front
;--------------------------------------------------------------

if (strlen(end_time) lt 9) then begin
  if (strpos(end_time,'+') lt 0) then					$
    end_time = strmid(start_time,0,9)+' '+end_time			$
  else begin
    rtime_end = rtime_start+double(strmid(end_time,1,strlen(end_time)-1))*60.0
    c_r_to_a, itime_end, rtime_end
    c_a_to_s, itime_end, end_time
  endelse
endif

;--------------------------------------------------------------
; Now figure out where in the file these things are, with the
; default to give the user everything
;--------------------------------------------------------------

c_s_to_a, itime_end, end_time
c_a_to_r, itime_end, rtime_end

n_start = where(time ge rtime_start)
if n_start(0) ge 0 then n_start = n_start(0) else n_start = 0

n_end = where(time ge rtime_end)
if n_end(0) ge 0 then n_end = n_end(0) else n_end = n_elements(time)-1

;--------------------------------------------------------------
; Put the contour data into data_1 array and get field name
;--------------------------------------------------------------

data_1 = reform(data(*,type_1,*,*))
field_1 = strcompress(fields(type_1))

;--------------------------------------------------------------
; figure out contour levels:
;--------------------------------------------------------------

maxi_array = fltarr(n_end-n_start+1)
mini_array = fltarr(n_end-n_start+1)

for i = n_start, n_end do begin
  mini_array(i-n_start) = min(data_1(i,*,*))
  maxi_array(i-n_start) = max(data_1(i,*,*))
endfor


mini_1  = min(data_1(n_start:n_end,*,*))
maxi_1  = max(data_1(n_start:n_end,*,*))

count = 0

onepercent = n_elements(data_1(n_start:n_end,*,*)) * 0.0001

while count lt onepercent do begin

  mini_1 = mini_1 * 0.95
  maxi_1 = maxi_1 * 0.95

  loc = where(data_1(n_start:n_end,*,*) le mini_1 or $
              data_1(n_start:n_end,*,*) ge maxi_1, count)

endwhile

if (mini_1 ge 0.0) then begin
  mini_1  = 0.0
  range_1 = maxi_1
endif else range_1 = max([abs(mini_1),maxi_1])*2.0
dc      = 10.0^fix(alog10(range_1/20.0))
factor  = 0.0
while (range_1 gt dc*20.0*factor) do factor=factor+0.05
dc = factor*dc
;dc = float(ask('contour level for '+field_1,strcompress(string(dc))))
if (mini_1) eq 0.0 then 						$
  levels_1 = findgen(21)*dc						$
else levels_1 = (findgen(21) - 10.0)*dc

;--------------------------------------------------------------
; Set up device
;--------------------------------------------------------------

setdevice, psfile, 'p', 4, 0.95

plotdumb

;--------------------------------------------------------------
; Read color table. Blue to white to red for -/+ data,
; red to white for + only data.
;--------------------------------------------------------------

readct, ncolors, getenv("IDL_EXTRAS")+"/blue_white_red.ct"	;	$
clevels = (ncolors-1)*findgen(21)/20.0 + 1

;--------------------------------------------------------------
; Set up plot sizes for the circles
;--------------------------------------------------------------

space = 0.01
pos_space, ppp, space, sizes, ny = ppp

;--------------------------------------------------------------
; Determine the character size in normalized coordinates
;--------------------------------------------------------------

dy  = float(!d.y_ch_size)/float(!d.y_size)

;--------------------------------------------------------------
; If we have a whole bunch of plots, scale the size down.
;--------------------------------------------------------------

if (ppp gt 12) then ch_size = 0.65 else ch_size = 1.0

ntimes = n_end - n_start + 1

image = fltarr(ntimes, nmlts)
latimage = fltarr(ntimes, nmlts)
area = fltarr(ntimes)
cpcp = fltarr(ntimes)
efield = fltarr(ntimes)
efe = reform(data(*,7,*,*))
efn = reform(data(*,8,*,*))
eft = sqrt(efn^2 + efe^2)

dl = (nmlts-1)/4
ef_dayside = reform(eft(*,dl+dl/2:dl*3-dl/2,4:nlats-1))

dl = 24.0 / (nmlts-1) * 15.0

for n = n_start, n_end do begin
  nn = n-n_start
  area(nn) = 0.0
  efield(nn) = max(eft(n,*,*)) ; max(data(n,0,*,*)) - min(data(n,0,*,*))
  cpcp(nn) = max(data(n,0,*,*)) - min(data(n,0,*,*))
  for j=0,nmlts-1 do begin
    maxi = max(abs(data_1(n,j,*)))
    mini = min(data_1(n,j,*))
    if (abs(mini) eq maxi) then image(nn,j) = mini else image(nn,j) = maxi
    loc = where(data_1(n,j,*) eq image(nn,j))
    latimage(nn,j) = lats(loc(0))

    ilat = loc(0)-1
    done = 0
    while (not done) do begin
      if (ilat le 0) then begin
        done = 1 
      endif else begin
        if (abs(data_1(n,j,ilat)) gt 0.25*abs(data_1(n,j,loc(0)))) then begin
          ilat = ilat - 1 
        endif else done = 1
      endelse
    endwhile

    if (ilat ge 0) then l = lats(ilat) else l = 90.0

    area(nn) = area(nn) + (90.0 - l)*111.0 * dl * 111.0 * sin((90.0 - l)*!dtor)
  endfor
endfor

image = shift(image,0,-(nmlts-1)/2)
image = alog10(image)

latimage = shift(latimage,0,-(nmlts-1)/2)

mini_1 = min(image)
maxi_1 = max(image)
range_1 = maxi_1-mini_1
image = (image - mini_1)/range_1 * float(ncolors-1)

mm_latimage = mm(latimage)

latimage = (latimage - min(latimage)) / $
           (max(latimage)-min(latimage)) * float(ncolors-1)

time_axis, rtime_start, rtime_end, btr, etr, $
           xtickname, xtitle, xtickv, xminor, xtickn

get_position, ppp, space, sizes, 0, pos, /rect
pos(2) = pos(2) - 0.05
pos(0) = pos(0) + 0.05
xsize = pos(2)-pos(0)-0.001
ysize = pos(3)-pos(1)-0.001
tv,image, pos(0)+0.0005,pos(1)+0.0005,xsize=xsize,ysize=ysize,/norm

plot, [btr,etr], [0.0,1.0], /nodata, $
      xstyle = 1, ystyle = 1, pos = pos, /noerase, $
      xtickname = strarr(30)+' ',	$
      xtitle = xtitle,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn, xrange = [btr,etr],  $
      ytickv = [0.0, 0.25, 0.5, 0.75, 1.0], $
      ytickname = ['12','18','24','06','12'], yticks = 4, yminor = 6, $
      ytitle = 'MLT'

ctpos = [pos(2)+0.01, pos(1), pos(2)+0.04, pos(3)]
ctmm = [mini_1, maxi_1]
cttitle = 'log('+field_1+')'
plotct, ncolors, ctpos, ctmm, cttitle, /right

get_position, ppp, space, sizes, 1, pos, /rect
pos(2) = pos(2) - 0.05
pos(0) = pos(0) + 0.05
xsize = pos(2)-pos(0)-0.001
ysize = pos(3)-pos(1)-0.001
tv,latimage, pos(0)+0.0005,pos(1)+0.0005,xsize=xsize,ysize=ysize,/norm

plot, [btr,etr], [0.0,1.0], /nodata, $
      xstyle = 1, ystyle = 1, pos = pos, /noerase, $
      xtickname = strarr(30)+' ', $
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn, xrange = [btr,etr],  $
      ytickv = [0.0, 0.25, 0.5, 0.75, 1.0], $
      ytickname = ['12','18','24','06','12'], yticks = 4, yminor = 6, $
      ytitle = 'MLT'

ctpos = [pos(2)+0.01, pos(1), pos(2)+0.04, pos(3)]
ctmm = mm_latimage
cttitle = 'Latitude'
plotct, ncolors, ctpos, ctmm, cttitle, /right

get_position, ppp, space, sizes, 2, pos, /rect
pos(2) = pos(2) - 0.05
pos(0) = pos(0) + 0.05
xsize = pos(2)-pos(0)-0.001
ysize = pos(3)-pos(1)-0.001

plot, time(n_start:n_end)-rtime_start, area/1.0e6, $
      xstyle = 1, ystyle = 1, pos = pos, /noerase, $
      xtickname = strarr(30)+' ',	$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn, xrange = [btr,etr]

xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, 'Area (m^2/1.0e12)', $
        orientation = 270, alignment = 0.5, /norm

get_position, ppp, space, sizes, 3, pos, /rect
pos(2) = pos(2) - 0.05
pos(0) = pos(0) + 0.05
xsize = pos(2)-pos(0)-0.001
ysize = pos(3)-pos(1)-0.001

plot, time(n_start:n_end)-rtime_start, ae(n_start:n_end,2), $
      xstyle = 1, ystyle = 1, pos = pos, /noerase, $
      xtickname = strarr(30)+' ',	$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn, xrange = [btr,etr]

xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, 'AL (nT)', $
        orientation = 270, alignment = 0.5, /norm

get_position, ppp, space, sizes, 4, pos, /rect
pos(2) = pos(2) - 0.05
pos(0) = pos(0) + 0.05
xsize = pos(2)-pos(0)-0.001
ysize = pos(3)-pos(1)-0.001

plot, time(n_start:n_end)-rtime_start, cpcp, $
      xstyle = 1, ystyle = 1, pos = pos, /noerase, $
      xtickname = strarr(30)+' ',	$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn, xrange = [btr,etr]

xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, 'Potential (kV)', $
        orientation = 270, alignment = 0.5, /norm

get_position, ppp, space, sizes, 5, pos, /rect
pos(2) = pos(2) - 0.05
pos(0) = pos(0) + 0.05
xsize = pos(2)-pos(0)-0.001
ysize = pos(3)-pos(1)-0.001

plot, time(n_start:n_end)-rtime_start, imf(n_start:n_end,2), $
      xstyle = 1, ystyle = 1, pos = pos, /noerase, $
      xtickname = xtickname,			$
      xtitle = xtitle,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn, xrange = [btr,etr], yrange = mm(imf(n_start:n_end,1:2))

oplot, time(n_start:n_end)-rtime_start, imf(n_start:n_end,1), linestyle = 2

oplot, [btr,etr], [0.0,0.0], linestyle = 1

xyouts, pos(2)+0.01, (pos(1)+pos(3))/2.0, 'Bz & By (nT)', $
        orientation = 270, alignment = 0.5, /norm

closedevice

end




