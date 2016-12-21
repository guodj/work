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

print, 'Postscript (y or n) ? '
que = ''
read, que
que = mklower(que)

data = dmax(m:n) - dmin(m:n)

mm = mm + m
if mm ge 60 then begin
  mm = mm - 60
  hh = hh + 1
endif
xtitle = 'Minutes past '+chopr('0'+tostr(hh),2)+chopr('0'+tostr(mm),2)+' UT'

plot, data, xtitle = xtitle, xstyle = 1, ytitle = 'kV',		$
	charsize = 1.2

n = n-m
errors = 100000.0

x = findgen(n+1)
y = fltarr(n+1)

ntimes = 10000

for ii=0,ntimes do begin

  if ii lt ntimes/3 then begin

    x1 = float(n/2)*randomu(seed,1)
    a = max(data)*randomu(seed,1)/2.0 + max(data)/2.0
    t = (float(n)-x1)*randomu(seed,1)

  endif else begin

    fac = 1.0-float(ii)/float(ntimes)
    x1p = fac*x1s
    ap = fac*as
    tp = fac*ts
    x1 = x1s + x1p*(randomu(seed,1)-0.5)
    a = as + ap*(randomu(seed,1)-0.5)
    t = ts + tp*(randomu(seed,1)-0.5)

  endelse

  x1 = x1(0)
  if x1 eq 0 then x1 = 1

  aved1 = 0.0
  for i=0,fix(x1) do aved1 = aved1 + data(i)/float(fix(x1)+1)

  a = a(0)
  t = t(0)
  b = (1.0-aved1/a)/exp((n-(x(x1)-x(x1)))/t)

  y(0:x1-1) = aved1
  y(x1:n) = a*(1.0-b*exp((n-(x(x1:n)-x(x1)))/t))

  error = 0.0
  for j=0,n do error = error + (y(j) - data(j))^2.0
  error = error/float(n)^2.0

  if error lt errors then begin
    x1s = x1
    as = a
    bs = b
    ts = t
    ys = y
    errors = error
    print, ii, as, bs, ts, errors
    plot, x, data, xstyle=1
    oplot, x, ys
    wait,0.2
  endif

endfor
print, 'done'

if (strmid(que,0,1) eq 'y') then setdevice, filein+'.ps','l',4,0.6

plot, data, xtitle = xtitle, xstyle = 1, ytitle = 'kV',		$
	charsize = 1.2
oplot, x, ys, linestyle = 1
oplot, [x(x1s),x(x1s)],[0.0,300.0]
m1 = fix(x(x1s))
s1 = fix((x(x1s)-fix(x(x1s)))*60.0)
m1 = m1 + mm
h1 = hh
if m1 ge 60 then begin
  h1 = h1 + 1
  m1 = m1 - 60
endif
t1 = chopr('0'+tostr(h1),2)+chopr('0'+tostr(m1),2)+':'+chopr('0'+tostr(s1),2)
xyouts, x1s+0.2,3.0,t1
xyouts, 1.0,max(data),'A : '+tostrf(as)
xyouts, 1.0,0.95*max(data),'B : '+tostrf(bs)
xyouts, 1.0,0.90*max(data),'T : '+tostrf(ts)
xyouts, 1.0,0.85*max(data),'error : '+tostrf(errors)

if !d.name eq 'PS' then begin
  device, /close
  set_plot,'X'
endif

end
