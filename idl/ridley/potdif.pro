dmin = fltarr(200)
dmax = fltarr(200)

print,'enter file name : '
filein = ''
read,filein

print, 'Enter start time of the file (hhmm) : '
stime = ''
read,stime

hh = float(strmid(stime,0,2))
mm = float(strmid(stime,2,4))

openr,1,filein

dumb1 = 0.0
dumb2 = 0.0
n = 0

while not eof(1) do begin

  readf,1,dumb1, dumb2
  dmin(n) = dumb1
  dmax(n) = dumb2
  n = n + 1

endwhile

close,1

plot, dmax(0:n-1)-dmin(0:n-1), xstyle=1

print, 'Enter first data point to use in fit : '
read, m
print, 'Enter last data point to use in fit : '
read, n

print, 'Postscript file name (return for screen) ? '
psfile = ''
read, psfile
if (strlen(psfile) gt 0) then setdevice, psfile,'l',4,0.6

print,'time resolution (1 or 5) :'
dt = 0
read,dt

title = ''
print, 'Enter a title :'
read, title

data = dmax(m:n) - dmin(m:n)
time = findgen(n_elements(data))*float(dt)

mm = mm + m*dt
if mm ge 60 then begin
  mm = mm - 60
  hh = hh + 1
endif
xtitle = 'Minutes past '+chopr('0'+tostr(hh),2)+chopr('0'+tostr(mm),2)+' UT'
plot, time, data, xtitle = xtitle, xstyle = 1, title = title, 		$
	charsize=1.2, ytitle = 'kV'

n = (n-m)*dt
errors = 100000.0
ntimes = 10000

for ii=0,ntimes do begin

  if ii lt 3000 then begin
    x1 = float(n)*randomu(seed,1)
    x2 = (float(n)-x1)*randomu(seed,1) + x1
  endif else begin
    fac = 1.0 - float(ii)/float(ntimes)
    xp = fac*float(n)/10.0
    x1 = x1s + xp*(randomu(seed,1)-0.5)
    x2 = x2s + xp*(randomu(seed,1)-0.5)
  endelse

  x1 = x1(0)
  x2 = x2(0)
  if x2-x1 lt dt then begin
    if x1 lt 1 then x2 = x2 + float(dt)
    if x1 gt n-1 then x1 = x1 - float(dt)
  endif

  aved1 = 0.0
  a = 0
  loc = where(time lt fix(x1),count)
  if count gt 0 then b = loc(count-1) else b = 0
  for i=a,b do 				$
    aved1 = aved1 + data(i)/float(b-a+1)

  aved2 = 0.0
  loc = where(time gt fix(x2),count)
  if count gt 0 then a = loc(0) else a = n_elements(data)-1
  b = n_elements(data)-1
  for i=a,b do			$
    aved2 = aved2 + data(i)/float(b-a+1)

  slope = (aved2-aved1)/(x2-x1)
  yi2 = aved2 - slope*x2
  yi1 = aved1 - slope*x1
  yi = (yi1+yi2)/2.0

  error = 0.0
  for j=0,n_elements(data)-1 do begin
    t = float(time(j))
    if t lt x1 then error = error + (aved1 - data(j))^2.0
    if t gt x2 then error = error + (aved2 - data(j))^2.0
    if (t ge x1 and t le x2) then 		$
	error = error + (slope*t+yi - data(j))^2.0
  endfor

  error = error/float(n_elements(data))^2.0

  if error lt errors then begin
    slopes = slope
    yis = yi
    x1s = x1
    x2s = x2
    aved2s = aved2
    aved1s = aved1
    errors = error
    print, ii, slopes, yis, x1s, x2s, errors
  endif

endfor

oplot, [0,x1s],[aved1s,aved1s], linestyle = 2
oplot, [x1s,x2s],[slopes*x1s+yis,slopes*x2s+yis], linestyle = 1
oplot, [x2s,n],[aved2s,aved2s], linestyle = 2
oplot, [x1s,x1s], [0.0,300.0]
oplot, [x2s,x2s], [0.0,300.0]

m1 = fix(x1s+0.5)
m1 = m1 + mm
h1 = hh
if m1 ge 60 then begin
  h1 = h1 + 1
  m1 = m1 - 60
endif

m2 = fix(x2s+0.5)
m2 = m2 + mm
h2 = hh
if m2 ge 60 then begin
  h2 = h2 + 1
  m2 = m2 - 60
endif

a1 = tostr(fix(aved1s))
a2 = tostr(fix(aved2s))

t1 = chopr('0'+tostr(h1),2)+chopr('0'+tostr(m1),2)+' UT'
t1s = '('+a1+' kV)'
t2 = chopr('0'+tostr(h2),2)+chopr('0'+tostr(m2),2)+' UT'
t2s = '('+a2+' kV)'

yp = 0.3*(max(data)-aved1s) + aved1s
xyouts, x1s-0.2,yp,t1, charsize = 1.1, alignment = 1.0
yp = 0.22*(max(data)-aved1s) + aved1s
xyouts, x1s-0.2,yp,t1s, charsize = 1.1, alignment = 1.0

if (x2s le 0.75*max(time)) then begin

  yp = 0.7*aved2s
  xyouts, x2s+0.2,yp,t2, charsize = 1.1, alignment = 0.0
  yp = 0.62*aved2s
  xyouts, x2s+0.2,yp,t2s, charsize = 1.1, alignment = 0.0

  yp = 0.1*max(data)
  xyouts, x2s-0.2,yp,'slope : ',alignment=1.0,charsize=1.1
  xyouts, x2s+0.2,yp,strmid(tostrf(slopes),0,4)+' kV/min',	$
	charsize = 1.1

endif else begin

  yp = 0.7*aved2s
  xyouts, x2s-0.2,yp,t2, charsize = 1.1, alignment = 1.0
  yp = 0.62*aved2s
  xyouts, x2s-0.2,yp,t2s, charsize = 1.1, alignment = 1.0

  yp = 0.1*max(data)
  xyouts, x2s-0.2,yp,'slope : '+strmid(tostrf(slopes),0,4)+' kV/min',	$
	alignment=1.0,charsize = 1.1

endelse

xyouts, 1.0,0.9*max(data),'err: '+strmid(tostrf(errors),0,4)+' kV',	$
	charsize=0.9

if !d.name eq 'PS' then begin
  device, /close
  set_plot,'X'
endif

end
