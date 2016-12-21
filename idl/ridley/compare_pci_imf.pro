
if (n_elements(iyear) eq 0) then iyear='2001' else iyear=tostr(iyear)
if (n_elements(imonth) eq 0) then imonth='1' else imonth=tostr(imonth)
if (n_elements(iday) eq 0) then iday='1' else iday=tostr(iday)

iyear  = fix(ask('year to compare',iyear))
imonth = fix(ask('month to compare',imonth))
iday   = fix(ask('day to compare',iday))

itime = [iyear,imonth,iday,0,0,0]
c_a_to_s, itime, stime

if (iyear lt 65)  then iyear = iyear + 2000
if (iyear lt 100) then iyear = iyear + 1900
syear = tostr(iyear)
smonth = chopr('0'+tostr(imonth),2)
sday   = chopr('0'+tostr(iday),2)

pcifile = '/f/amie/Data/Indices/ftp.dmi.dk/data/anonymousdata/pub/wdcc1/indices/pcn/p15n'+syear+'.wdc'

read_pci_wdc, pcifile, pci, time


cdffile = '/f/amie/Data/Indices/IMF/cdaweb.gsfc.nasa.gov/'+ $
          'pub/istp/ace/mfi_h0/'+syear+'/ac_h0_mfi_'+ $
          syear+smonth+sday+'_v04.cdf'

id = cdf_open(cdffile)
re = cdf_inquire(id)

npts1 = re.maxrec+1

cdf_varget, id, 0, imfepoch1, rec_count = npts1
cdf_varget, id, 5, imf1, rec_count = npts1

cdf_close, id

; --------------------------------------------------------------

itime2 = [iyear,imonth,iday+1,0,0,0]
c_a_to_r, itime2, rtime2
c_r_to_a, itime2, rtime2
syear2 = tostr(itime2[0])
smonth2 = chopr('0'+tostr(itime2[1]),2)
sday2   = chopr('0'+tostr(itime2[2]),2)

cdffile2 = '/f/amie/Data/Indices/IMF/cdaweb.gsfc.nasa.gov/'+ $
          'pub/istp/ace/mfi_h0/'+syear+'/ac_h0_mfi_'+ $
          syear2+smonth2+sday2+'_v04.cdf'

id = cdf_open(cdffile2)
re = cdf_inquire(id)

npts2 = re.maxrec+1

cdf_varget, id, 0, imfepoch2, rec_count = npts2
cdf_varget, id, 5, imf2, rec_count = npts2

cdf_close, id

npts = npts1+npts2
imfepoch = dblarr(3,npts)
imf      = fltarr(3,npts)

for i=0,npts1-1 do begin
  imfepoch(*,i) = imfepoch1(*,i)
  imf(*,i)      = imf1(*,i)
endfor

for i=0,npts2-1 do begin
  ii = i + npts1
  imfepoch(*,ii) = imfepoch2(*,i)
  imf(*,ii)      = imf2(*,i)
endfor

; --------------------------------------------------------------

cdffile = '/f/amie/Data/Indices/IMF/cdaweb.gsfc.nasa.gov/'+ $
          'pub/istp/ace/swe_h0/'+syear+'/ac_h0_swe_'+ $
          syear+smonth+sday+'_v05.cdf'

id = cdf_open(cdffile)
re = cdf_inquire(id)

nptsswe1 = re.maxrec+1

cdf_varget, id,  0, sweepoch1, rec_count = nptsswe1
cdf_varget, id,  5, np1, rec_count = nptsswe1
cdf_varget, id, 13, vgsm1, rec_count = nptsswe1
cdf_varget, id, 17, scpgsm1, rec_count = nptsswe1

cdf_close, id

cdffile2 = '/f/amie/Data/Indices/IMF/cdaweb.gsfc.nasa.gov/'+ $
          'pub/istp/ace/swe_h0/'+syear+'/ac_h0_swe_'+ $
          syear2+smonth2+sday2+'_v05.cdf'

id = cdf_open(cdffile2)
re = cdf_inquire(id)

nptsswe2 = re.maxrec+1

cdf_varget, id,  0, sweepoch2, rec_count = nptsswe2
cdf_varget, id,  5, np2, rec_count = nptsswe2
cdf_varget, id, 13, vgsm2, rec_count = nptsswe2
cdf_varget, id, 17, scpgsm2, rec_count = nptsswe2

cdf_close, id

nptsswe = nptsswe1+nptsswe2
sweepoch = dblarr(3,nptsswe)
np       = fltarr(3,nptsswe)
vgsm     = fltarr(3,nptsswe)
scpgsm   = fltarr(3,nptsswe)

for i=0,nptsswe1-1 do begin
  sweepoch(*,i) = sweepoch1(*,i)
  np(*,i)       = np1(*,i)
  vgsm(*,i)     = vgsm1(*,i)
  scpgsm(*,i)   = scpgsm1(*,i)
endfor

for i=0,nptsswe2-1 do begin
  ii = i + nptsswe1
  sweepoch(*,ii) = sweepoch2(*,i)
  np(*,ii)       = np2(*,i)
  vgsm(*,ii)     = vgsm2(*,i)
  scpgsm(*,ii)   = scpgsm2(*,i)
endfor

; --------------------------------------------------------------

imftime = dblarr(npts)
swetime = dblarr(nptsswe)

vx = vgsm(0,*)
xx = scpgsm(0,*)
delay = fltarr(nptsswe)
loc = where(vx gt -10000.0,count)
if count gt 0 then begin
  delay(loc) = xx(loc)/(-vx(loc))
endif else delay(*) = 3600.0

loc = where(delay eq 0.0,count)

if count gt 0 then begin
  for i=0,count-1 do begin
    ii = loc(i)
    if (ii eq 0) then begin
      loc2 = where(delay ne 0.0)
      delay(ii) = delay(loc2(0))
    endif else delay(ii) = delay(ii-1)
  endfor
endif

t1 = reform(sweepoch(0,*))
t2 = reform(imfepoch(0,*))

imfdelay = fltarr(npts)

for i=0,npts-1 do begin

  d = abs(t2(i)-t1)
  loc = where(d lt 600.0, count)
  if count gt 1 then imfdelay(i) = mean(delay(loc)) else begin
    loc = where(d eq min(d), count)
    imfdelay(i) = delay(loc(0))
  endelse
endfor



for i=0,npts-1 do begin

  cdf_epoch, imfepoch(0,i), year, month, day, hour, minute, second, /break

  itime = [year,month,day,hour,minute,second]
  c_a_to_r, itime, rtime
  imftime(i) = rtime + imfdelay(i)

endfor

for i=0,nptsswe-1 do begin

  cdf_epoch, sweepoch(0,i), year, month, day, hour, minute, second, /break

  itime = [year,month,day,hour,minute,second]
  c_a_to_r, itime, rtime
  swetime(i) = rtime + delay(i)

endfor

mini = min(imftime)
maxi = max(imftime)
loc = where(time ge mini and time le maxi, count)

if count gt 0 then begin

  ; this means we have overlapping data sets!

  pcinew = pci(loc)
  timenew = time(loc)

  dt = (timenew(1) - timenew(0))/2.0
  imfnew = fltarr(3, count)
  vswnew = fltarr(3, count)
  scpnew = fltarr(3, count)
  npnew  = fltarr(count)

  for i=0,count-1 do begin

    d = abs(timenew(i) - imftime)
    l = where(d lt dt and imf(0,*) gt -1000.0,c)
    if (c gt 0) then begin
      for j=0,2 do begin
        imfnew(j,i) = mean(imf(j,l))
      endfor
    endif else imfnew(*,i) = -1000.0

    d = abs(timenew(i) - swetime)
    l = where(d lt dt and vgsm(0,*) gt -10000.0,c)
    if (c gt 0) then begin
      for j=0,2 do begin
        vswnew(j,i) = mean(vgsm(j,l))
        scpnew(j,i) = mean(scpgsm(j,l))
      endfor
      npnew(i) = mean(np(0,l))
    endif else begin
      vswnew(*,i) = -10000.0
      npnew(i) = -10000.0
      scpnew(*,i) = -10000.0
    endelse

  endfor

endif

loc = where(vswnew(0,*) gt -10000.0 and imfnew(0,*) gt -1000.0, count)

if (count eq 0) then begin
  vswnew(0,*) = -1100.0
  npnew(*) = 2.0
  loc = where(vswnew(0,*) gt -10000.0 and imfnew(0,*) gt -1000.0, count)
endif

bt = sqrt(imfnew(1,loc)^2 + imfnew(2,loc)^2)
b  = sqrt(imfnew(0,loc)^2 + imfnew(1,loc)^2 + imfnew(2,loc)^2)

theta = acos(imfnew(2,loc) / bt)

ekl = vswnew(0,loc) * bt * sin(theta/2.0)/1000.0
pci = pcinew(loc)

vxsw = vswnew(0,loc)
nsw = npnew(loc)

mu0 = 4.0*!pi*1.0e-7
k   = 1.3807e-23
mp  = 1.6726e-27

pdyn = nsw * 100.0^3 * mp * (vxsw*1000.0)^2
pb   = 2.0*(b*1.0e-9)^2/mu0

esw = abs(ekl)
psw = (pdyn+pb) * 1.0e9

;psw = pdyn * 1.0e9

d = 1.0
f = 1.0

s0 = 12.0
zeta = 21.7/4

siscoe = 57.6 * esw * (psw)^(1.0/3.0) * d * f / $
	  ((psw)^(1.0/2.0)*d + 0.0125*zeta*s0*esw*f)

va = b*1.0e-9 / sqrt(mu0*nsw*100.0^3*mp)
ma = abs(vxsw*1000.0)/va

setdevice, 'naga'+syear+smonth+sday+'_siscoe2.ps','p',4

ppp = 3
space = 0.075
pos_space, ppp, space, sizes, ny = ppp

plotdumb

get_position, ppp, space, sizes, 0, pos, /rect
pos(0) = pos(0) + 0.2
pos(2) = pos(2) - 0.1

maxi = max([max(siscoe),max(19.35*pci+8.78),450])

xm = fix(max(esw)/5.0+1)*5.0

plot, esw, 19.35*pci+8.78, psym = 1, yrange = [0,maxi], $
      pos = pos, /noerase, ytitle = 'CPCP (kV)', xtitle = 'E!DM!N (mV/m)', $
      title = strmid(stime,0,9), ystyle = 1, xrange = [0,xm], xstyle=1

loc_ma = where(ma lt 2.0,count)
if (count gt 0) then oplot, esw(loc_ma), abs(19.35*pci(loc_ma)+8.78), psym = 2
oplot, esw, siscoe, psym = 4

loc = where(19.35*pci+8.78 gt 0.0, c)

esw2 = esw(loc)
cpcp2 = 19.35*pci(loc)+8.78

vip = 20.0
sp = 120.0
efold = 5.0

;efold = 100.0

cpcp3 = sp-(sp-vip)*exp(-esw2/efold)

error = mean(abs(cpcp2 - cpcp3))/float(c)
es = error
efs = efold
sps = sp
vis = vip

for i=0,10000 do begin

  sp    = sps + (randomu(s,1) - 0.5)*10.0
  efold = efs + (randomu(s,1) - 0.5)*1.0

;efold = 100.0

  vip   = vis + (randomu(s,1) - 0.5)*2.0
  sp = sp(0)
  efold = efold(0)
  vip = vip(0)
  cpcp3 = sp-(sp-vip)*exp(-esw2/efold)
  error = mean(abs(cpcp2 - cpcp3))/float(c)

  if (error lt es) then begin
    sps = sp
    efs = efold
    es  = error
    vis = vip
  endif

endfor

print, sps, efs, vis, es

en = findgen(200)/5.0
cpcp3 = sps-(sps-vis)*exp(-en/efs)

oplot, en, cpcp3, thick = 3

ssps = string(sps, format='(f5.1)')
sefs = string(efs, format='(f4.1)')
svis = string(vis, format='(f4.1)')
ses  = string(es, format='(f4.2)')

loc = where(19.35*pci+8.78 gt 0.0 and esw gt 2.0, c)

esw2 = esw(loc)
cpcp2 = 19.35*pci(loc)+8.78

cpcp3 = sps-(sps-vis)*exp(-en/efs)
error = mean(abs(cpcp2 - cpcp3))/float(c)

print, error
ses2 = string(error, format='(f4.2)')

out1 = 'CPCP = '+ssps+'-('+ssps+'-'+svis+')e'
out2 = '(-E!Dm!N/'+sefs+')'
out3 = ';E='+ses+' ['+ses2+']'

x1 = 3.0/40.0 * xm
x2 = 5.0/40.0 * xm

oplot, [x1, x2], [410.0, 410.0]

x1 =  5.5/40.0 * xm
x2 = 23.7/40.0 * xm
x3 = 30.0/40.0 * xm

xyouts, x1, 400.0, out1
xyouts, x2, 415.0, out2, charsize = 0.9
xyouts, x3, 400.0, out3

vip = 20.0
sp = 120.0
efold = 20.0

cpcp3 = sp-(sp-vip)*exp(-esw2/efold)

error = mean(abs(cpcp2 - cpcp3))/float(c)
es = error
efs = efold
sps = sp
vis = vip

for i=0,10000 do begin

  sp    = sps + (randomu(s,1) - 0.5)*10.0
  efold = efs + (randomu(s,1) - 0.5)*1.0
  vip   = vis + (randomu(s,1) - 0.5)*2.0
  sp = sp(0)
  efold = efold(0)
  vip = vip(0)
  cpcp3 = sp-(sp-vip)*exp(-esw2/efold)
  error = mean(abs(cpcp2 - cpcp3))/float(c)

  if (error lt es) then begin
    sps = sp
    efs = efold
    es  = error
    vis = vip
  endif

endfor

print, sps, efs, vis, es

en = findgen(200)/5.0
cpcp3 = sps-(sps-vis)*exp(-en/efs)

oplot, en, cpcp3, thick = 3, linestyle = 1

oplot, [2.0,2.0],[0.0,340.0], linestyle = 2

ssps = string(sps, format='(f5.1)')
sefs = string(efs, format='(f4.1)')
svis = string(vis, format='(f4.1)')
ses  = string(es, format='(f4.2)')

loc = where(19.35*pci+8.78 gt 0.0, c)

esw2 = esw(loc)
cpcp2 = 19.35*pci(loc)+8.78

cpcp3 = sps-(sps-vis)*exp(-en/efs)
error = mean(abs(cpcp2 - cpcp3))/float(c)

print, error
ses2 = string(error, format='(f4.2)')

out1 = 'CPCP = '+ssps+'-('+ssps+'-'+svis+')e'
out2 = '(-E!Dm!N/'+sefs+')'
out3 = ';E='+ses2+' ['+ses+']'

x1 = 3.0/40.0 * xm
x2 = 5.0/40.0 * xm

oplot, [x1, x2], [370.0, 370.0], linestyle = 1

x1 =  5.5/40.0 * xm
x2 = 23.7/40.0 * xm
x3 = 30.0/40.0 * xm

xyouts, x1, 360.0, out1
xyouts, x2, 375.0, out2, charsize = 0.9
xyouts, x3, 360.0, out3

;get_position, ppp, space, sizes, 1, pos, /rect
;pos(0) = pos(0) + 0.2
;pos(2) = pos(2) - 0.1
;
;loc = where(esw ge 7.0 and esw lt 10.0)
;plot, bt(loc), pci(loc), psym = 1, yrange = [1,9], xrange = [25,45], $
;      ystyle = 1, xstyle = 1, $
;      pos = pos, /noerase, xtitle = 'B!DT!N (nT)', ytitle = 'PCI (mV/m)'
;
;loc = where(esw ge 10.0 and esw lt 15.0)
;oplot, bt(loc), pci(loc), psym = 2
;
;loc = where(esw ge 15.0 and esw lt 20.0)
;oplot, bt(loc), pci(loc), psym = 4
;
;loc = where(esw ge 20.0)
;oplot, bt(loc), pci(loc), psym = 5
;
;psw = pdyn * 1.0e9
;
;get_position, ppp, space, sizes, 2, pos, /rect
;pos(0) = pos(0) + 0.2
;pos(2) = pos(2) - 0.1
;
;loc = where(esw ge 7.0 and esw lt 10.0)
;plot, psw(loc), pci(loc), psym = 1, yrange = [1,9], xrange = [0,6], $
;      ystyle = 1, xstyle = 1, $
;      pos = pos, /noerase, xtitle = 'P!Dsw!N (nP)', ytitle = 'PCI (mV/m)'
;
;loc = where(esw ge 10.0 and esw lt 15.0)
;oplot, psw(loc), pci(loc), psym = 2
;
;loc = where(esw ge 15.0 and esw lt 20.0)
;oplot, psw(loc), pci(loc), psym = 4
;
;loc = where(esw ge 20.0)
;oplot, psw(loc), pci(loc), psym = 5

closedevice

end
