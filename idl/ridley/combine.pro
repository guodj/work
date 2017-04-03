
filelist = findfile('*.map.pot')
display, filelist
if (n_elements(iFile) eq 0) then iFile = 0
iFile = fix(ask('superdarn file to convert',tostr(iFile)))
file = filelist(iFile)
read_superdarn, file, time, potential, lats, mlts

c_r_to_a, itime, time(0)
fileout = 'sd'+tostr(itime(0),4)+tostr(itime(1),2)+tostr(itime(2),2)+'.bin'
print, fileout

dlat = lats(0,1) - lats(0,0)
dlon = (mlts(1,0) - mlts(0,0))*180.0/12.0

area = 120000.0^2 * dlat * dlon * cos(lats*!dtor)

lats = reform(lats(0,*))
mlts = reform(mlts(*,0))

nLatsSD = n_elements(lats)
nMlts = n_elements(mlts)

MinLat = min(lats)-10.0
nLatsAdd = 10.0/dLat
lats = [findgen(nLatsAdd)+MinLat,lats]
nLats = n_elements(lats)


imffile = 'imf.dat'
imf_read, imffile, imftime, mag, vel, den, temp, nPts, notes

calc_ec, imftime, mag, vel, den, Ec4Hour

load_ovation, ovation

nTimesSD = n_elements(time)

nTimesSD = 20

by = fltarr(nTimesSD+2)
bz = fltarr(nTimesSD+2)
vx = fltarr(nTimesSD+2)
ec = fltarr(nTimesSD+2)
dipoletilt = fltarr(nTimesSD+2)

Vars = ['Potential (V)', $
        'Energy Flux (ergs/cm2)','Mean Energy (eV)','Number Flux (/cm2)']

nVars = n_elements(Vars)
data = fltarr(nTimesSD+2, nVars, nmlts, nlats)

cpcp = fltarr(nTimesSD+2)
hp   = fltarr(nTimesSD+2,2)
imf = fltarr(nTimesSD+2,4)

newtime = dblarr(nTimesSD+2)

for i=0,nTimesSD+1 do begin

   ii = min([max([i-1,0]),nTimesSD-1])

   t = time(ii)
   c_r_to_a, itime, t
   if (i eq 0) then begin
      itime(3:5) = 0
      c_a_to_r, itime, t
   endif
   if (i eq nTimesSD+1) then begin
      itime(3:5) = 0
      c_a_to_r, itime, t
      t = t + 24.0*3600.0
      c_r_to_a, itime, t
   endif

   newtime(i) = t

   ;dipoletilt(i) = tilt(itime(0),itime(1),itime(2),itime(3),itime(4),itime(5))

   d = abs(imftime - t)
   l = where(d eq min(d))
   ec(i) = ec4Hour(l(0))
   by(i) = mag(1,l(0))
   bz(i) = mag(2,l(0))
   vx(i) = vel(0,l(0))

   print, itime
;   print, by(i), bz(i)

   run_ovation, ovation, ec(i), t, north, south

   interpolate_ovation, north.eflux, north.lats, north.mlts, $
                        lats, mlts, eflux

   interpolate_ovation, north.avee, north.lats, north.mlts, $
                        lats, mlts, avee

   interpolate_ovation, north.nflux, north.lats, north.mlts, $
                        lats, mlts, nflux

   data(i,0,*,nLatsAdd:nLats-1) = potential(ii,*,*)
   for iMlt = 0, nMlts-1 do begin
      for iLat = 0, nLatsAdd-1 do begin
         fac = float(iLat)/float(nLatsAdd)
         data(i,0,iMlt,iLat) = data(i,0,iMlt,nLatsAdd) * fac
      endfor
   endfor

   l = where(eflux lt 0.1,c)
   if (c gt 0) then eflux(l) = 0.1

   l = where(avee lt 0.5,c)
   if (c gt 0) then avee(l) = 0.5

   l = where(nflux lt 10.0,c)
   if (c gt 0) then nflux(l) = 10.0

   data(i,1,*,*) = eflux
   data(i,2,*,*) = avee
   data(i,3,*,*) = nflux

   data(i,1,*,0:nLatsAdd-1) = 0.0
   data(i,2,*,0:nLatsAdd-1) = 0.0
   data(i,3,*,0:nLatsAdd-1) = 0.0

   hp(i,0) = total(area * eflux/1000.0)/1.0e9

   cpcp(i) = max(potential(ii,*,*))-min(potential(ii,*,*))
   imf(i,0) = vx(i)
   imf(i,1) = 0
   imf(i,2) = by(i)
   imf(i,3) = bz(i)

   print, "power, cpcp : ",hp(i,0), cpcp(i)

endfor

Version = 0.1
colats = 90.0-lats

amie_write_binary, fileout, Vars, lats, mlts, newtime, data, $
                   imf = imf, ae = ae, dst = dst, hpi = hp, cpcp = cpcp, $
                   Version = Version


end
