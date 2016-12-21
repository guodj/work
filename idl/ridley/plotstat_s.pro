space = 0.06
ppp = 4
pos_space, ppp, space, sizes

master=''
print,'Enter master station file :'
read, master

close,1
openr,1,master
;openr,1,'/d/data/d.mag/d.flatfiles/d.nocurve/amiein.910518a'

line = ''
readf,1,line
nstat = fix(strmid(line,0,3))

readf,1,line

lat = fltarr(nstat)
lon = fltarr(nstat)
list = intarr(50)
slist = strarr(50)
n = 0

for i=0,nstat-1 do begin

  readf,1,line
  lat(i) = float(strmid(line,36,6))
  lon(i) = float(strmid(line,43,6))*!pi/180.0
  if (lat(i) gt 71.0) and (lat(i) lt 76.0) then begin
    list(n) = i
    slist(n) = strmid(line,4,3)
    n = n+1
  endif

endfor

close,1

list = list(0:n-1)

print,"Enter time (hh:mm) :"
read,line
h = float(strmid(line,0,2))*60.0
m = (float(strmid(line,3,2))+h - 4.0*60.0 - 47.0)/4.0 * !pi/180.0
lon = lon + m - !pi/2.0

r = 40.0

get_position, ppp, space, sizes, 0, pos
xsize = pos(2) - pos(0)
ysize = pos(3) - pos(1)
xs = 0.5 - xsize/2.0
ys = 0.5 - ysize/2.0
pos = [xs,ys,xs+xsize,ys+ysize]

plot, [-r,r],[-r,r], xstyle=5, ystyle=5, /nodata, pos = pos

t = findgen(361)*!pi/180.0
oplot, r*cos(t),r*sin(t)
for i=0.0,r-10.0,10.0 do oplot, i*cos(t),i*sin(t), linestyle=2
oplot, [-r,r],[0.0,0.0]
oplot, [0.0,0.0],[-r,r]

xyouts, 0.0, -r-r/15, '00', alignment = 0.5
xyouts, 0.0, r+r/100, '12', alignment = 0.5
xyouts, -r-r/100, -r/40, '18', alignment = 1.0
xyouts, +r+r/100, -r/40, '06', alignment = 0.0
xyouts, +r-r/20, +r+r/100, line+' UT', alignment = 1.0
xyouts, 0.5, 0.8, 'Approximate locations of Magnetometers', 	$
	alignment=0.5, /norm

lat = 90 + lat
loc = where(lat le r, count)

if count gt 0 then begin
  oplot, [lat(loc)*cos(lon(loc))], [lat(loc)*sin(lon(loc))], psym = 4
endif

oplot, [0.0,r*cos(min(lon(loc)))],[0.0,r*sin(min(lon(loc)))]
oplot, [0.0,r*cos(max(lon(loc)))],[0.0,r*sin(max(lon(loc)))]

lon(list) = lon(list)-!pi/2.0
loc = where(lon(list) lt -2.0*!pi, count)
if count gt 0 then lon(list(loc)) = lon(list(loc)) + 2.0*!pi
loc = where(lon(list) gt 2.0*!pi, count)
if count gt 0 then lon(list(loc)) = lon(list(loc)) - 2.0*!pi
loc = where(lon(list) lt -1.0*!pi, count)
if count gt 0 then lon(list(loc)) = lon(list(loc)) + 2.0*!pi

for i=0,n_elements(list)-1 do begin

  loc = where(abs(lon(list)) eq min(abs(lon(list))))
  ii = list(loc(0))
;  xyouts, lat(ii)*cos(lon(ii)+!pi/2.0), lat(ii)*sin(lon(ii)+!pi/2.0), 	$
;	tostr(i+1)
  lon(ii) = lon(ii)*12.0/!pi + 12.0
  h = chopr('0'+tostr(fix(lon(ii))),2)
  m = chopr('0'+tostr(fix(60.0*(lon(ii) - fix(lon(ii))))),2)
;  xyouts, 0.75, 0.80-i*0.04, 			$
;	chopr('0'+tostr(i+1),2)+'. '+slist(loc(0))+' '+h+m+' MLT', /norm
  lon(ii) = 500.0

endfor

end