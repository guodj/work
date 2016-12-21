
if n_elements(data) eq 0 then getv = 1 else begin
  que = ''
  read,'Would you like to re-read in the data (y/n)?',que
  if strmid(que,0,1) eq 'y' then getv = 1 else getv = 0
endelse

if getv then begin

  filein = ask('data file from TIEGCM run','tiegcm.3d.save')

  que = ask('y if this file an ascii file','n')

  if (strmid(que,0,1) eq 'y') then begin
    ntmax = 350
    rdascii, filein, data, vars, utdes, axdes, xax, yax, ntmax
  endif else begin
    if strpos(filein,'save') eq -1 then 				$
      rdbin, filein, data, vars, utdes, axdes, xax, yax			$
    else restore, filein
  endelse

endif

n = n_elements(vars)

for i=0,n-1 do begin
  vars(i) = stripend(vars(i))
  axdes(i) = stripend(axdes(i))
  utdes(i) = stripend(utdes(i))
endfor

done = 0
i = 1
while not done do if vars(i) eq vars(0) then i = i + 1 else done = 1
nz = i

for i=0,n-1,nz do print, i,'. ',vars(i)

psfile = ask('ps file name','test.ps')

minlat = float(ask('minimum latitude to plot','40.0'))
maxran = 90.0-minlat

que=ask('y to smooth the images (y/n)','n')
if strmid(que,0,1) eq 'y' then smooth = 1 else smooth = 0

print, 'You can also smooth the image by expanding the array size.'
xfac = fix(ask('expansion factor (1 - normal size, 3 - 3x normal ...)','1'))

nx = n_elements(data(0,*,0))
ny = n_elements(data(0,0,*))

ed = fltarr(nx,ny)
et = fltarr(nx,ny)
nm = fltarr(nx,ny)
td = fltarr(nx,ny)
vx = fltarr(nx,ny)
vy = fltarr(nx,ny)
ux = fltarr(nx,ny)
uy = fltarr(nx,ny)

b = fltarr(nx,ny)
bz = fltarr(nx,ny)
bx = fltarr(nx,ny)

re = 6371.0
mn2 = 12.0*2.0
mo  = 16.0
mo2 = mo*2

print, 'There are ',tostr(nz),' height slices.'
print, 'Enter height slice to examine (-1 for height integrated) : '
hs = 0
read, hs

if hs eq -1 then begin

  hss  = 0
  hse  = nz-1
  hint = 1

  h1 = float(strmid(utdes(0),22,6))
  dh = 1000.0*(float(strmid(utdes(1),22,6))-h1)
  print, 'dh = ',dh,' meters.'

endif else begin

  hss  = hs
  hse  = hs
  hint = 0

endelse

q    = fltarr(nx,ny)
sigp = fltarr(nx,ny)
sigh = fltarr(nx,ny)

ut = float(strmid(utdes(0),9,5))
hs = fix(ut)
ms = fix(60.0*(ut - hs))
uts = chopr('0'+tostr(hs),2)+chopr('0'+tostr(ms),2)+' UT'

for hs=hss,hse do begin

  height = float(strmid(utdes(hs),22,6))
  print, 'Processing height : ',tostr(height),' km'

  ied  = 0*nz + hs
  iet  = 7*nz + hs
  inzw = 1*nz + hs
  inmw = 2*nz + hs
  iizw = 3*nz + hs
  iimw = 4*nz + hs
  in2d = 8*nz + hs
  io2d = 9*nz + hs
  iod  =10*nz + hs

  ed(*,*) = data(ied,*,*)
  et(*,*) = data(iet,*,*)
  vx(*,*) = data(iizw,*,*)
  vy(*,*) = data(iimw,*,*)
  ux(*,*) = data(inzw,*,*)
  uy(*,*) = data(inmw,*,*)

  td(*,*) = data(in2d,*,*) + data(io2d,*,*) + data(iod,*,*)

  nm(*,*) = (data(in2d,*,*)*mn2 + data(io2d,*,*)*mo2 + 	$
	     data(iod,*,*)*mo)/td

  for i=0,nx-1 do for j=0,ny-1 do begin
    z = height
    t = yax(j)*!pi/180.0
    re3 = ((z+re)/re)^(-3.0)
    b(i,j) = 0.3*re3*(1.0+3.0*(sin(t))^2.0)^0.2
    bz(i,j) = -0.6*re3*sin(t)
    bx(i,j) = 0.3*re3*cos(t)
  endfor

  ; --------------
  ; Electric Field
  ex = -vy*bz*1.0e-5  ; gauss to T
  ey =  vx*bz*1.0e-5  ; gauss to T
  ez = -vy*bx*1.0e-5  ; gauss to T
  ; --------------

  ; --------------
  ; UxB
  uxbx =  uy*bz*1.0e-5  ; gauss to T
  uxby = -vx*bz*1.0e-5  ; gauss to T
  uxbz =  vy*bx*1.0e-5  ; gauss to T
  ; --------------

  ; --------------
  ; E + UxB
  tex = ex + uxbx
  tey = ey + uxby
  tez = ez + uxbz
  ; --------------

  ; --------------
  ; squared
  e2 = tex^2.0 + tey^2.0 + tez^2.0
  ; --------------

  e = 1.609e-16
  sm = 9.1094e-31
  bm = nm*1.6726e-27

  ; ---------------------
  ; Collision Frequencies
  vin = 2.6e-9 * (td+ed)*(nm^(-0.5))
  ven = 5.4e-10 * td * et^0.5
  vei = (34.0 + 4.18*log((et^3.0)/ed))*ed*et^(-3.0/2.0)
  vi = vin
  ve = ven + vei
  ; ---------------------

  ; ---------------------
  ; Gyro Frequencies
  oi = 2.0*!pi*(1.52e3) * b/nm
  oe = -2.0*!pi*(2.80e6) * b
  ; ---------------------

  ; ---------------------
  ; Specific conductivity
  so = ed*(e^2.0)*((1.0/(sm*ve)) + (1.0/(bm*vi)))

  ; ---------------------
  ; Hall and Pederson Conductivities
  p1 = ve/(oe^2.0+ve^2.0)
  p2 = vi/(oi^2.0+vi^2.0)
  h1 = oi/(oi^2.0+vi^2.0)
  h2 = oe/(oe^2.0+ve^2.0)
  sigp = sigp + ed*e^2.0*(p1/sm + p2/bm)
  sigh = sigh + ed*e^2.0*(h1/bm - h2/sm)

  ; --------------
  ; Joule Heating
  q = q + (ed*e^2.0*(p1/sm + p2/bm))*e2
  ; --------------

endfor

if hint then begin
  q = q * dh
  sigh = sigh * dh
  sigp = sigp * dh
endif

if smooth then begin

  hf = 4.0
  lf = 1.0
  df = hf+4.0*lf

  i = 0
  for j=1,ny-2 do begin

    q(i,j) = (hf*q(i,j) + lf*q(nx-1,j) + 		$
			  lf*q(i+1,j) +			$
			  lf*q(i,j-1) +			$
			  lf*q(i,j+1))/df

    sigp(i,j) = (hf*sigp(i,j) + lf*sigp(nx-1,j) + 	$
			  lf*sigp(i+1,j) +		$
			  lf*sigp(i,j-1) +		$
			  lf*sigp(i,j+1))/df

  endfor

  for i=1,nx-2 do for j=1,ny-2 do begin

    q(i,j) = (hf*q(i,j) + lf*q(i-1,j) + 		$
			  lf*q(i+1,j) +			$
			  lf*q(i,j-1) +			$
			  lf*q(i,j+1))/df

    sigp(i,j) = (hf*sigp(i,j) + lf*sigp(i-1,j) + 	$
			  lf*sigp(i+1,j) +		$
			  lf*sigp(i,j-1) +		$
			  lf*sigp(i,j+1))/df

  endfor

  i = nx-1
  for j=1,ny-2 do begin

    q(i,j) = (hf*q(i,j) + lf*q(i-1,j) + 		$
			  lf*q(0,j) +			$
			  lf*q(i,j-1) +			$
			  lf*q(i,j+1))/df

    sigp(i,j) = (hf*sigp(i,j) + lf*sigp(i-1,j) + 	$
			  lf*sigp(0,j) +		$
			  lf*sigp(i,j-1) +		$
			  lf*sigp(i,j+1))/df

  endfor

endif

;expand_array, xax, xfac
;expand_array, yax, xfac*2
;expand_array, q, [xfac,xfac*2]
;expand_array, sigp, [xfac,xfac*2]

nq = float(n_elements(q))

if strlen(psfile) gt 0 then setdevice, psfile,'l',4

ct_dir = getenv('IDL_EXTRAS')
ctname = ct_dir+'jim6.ct'
readct, ncolors, ctname

plotdumb

ppp = 4
space = 0.05
pos_space, ppp, space, sizes
get_position, ppp, space, sizes, 0, pos1
get_position, ppp, space, sizes, 2, pos2

lon = xax + (ut-12.0)*360.0/24.0 + 90.0

if hint then nf = 2.0e-3 else begin

  nf = 10.0^float(round(log(max(q)))+1.0)
  done = 0
  while not done do begin
    loc = where(q gt nf,count)
    per = float(count)/nq
    if per gt 5.0/nq then nf = nf*1.25		$
    else if per lt 1.0/nq then nf = nf*0.75		$
    else done = 1
  endwhile

endelse

image = q*float(ncolors)/nf
mltpoly, image, yax, lon, maxran, pos1
plotmlt, maxran, /white

xyouts, pos1(0), pos1(3), uts, /norm

posc = [pos1(2)+0.12,pos1(1),pos1(2)+0.17,pos1(3)]
if hint then begin
  cttitle = 'Joule Heating mW/m2' 		
  xyouts, (pos1(0)+pos1(2))/2.0, pos1(3)+0.05, 		$
  'Height Integrated Values', /norm, alignment = 0.5
endif else begin
  cttitle = 'Joule Heating mW/m3'
  xyouts, (pos1(0)+pos1(2))/2.0, pos1(3)+0.05, 		$
    'Values at '+tostr(fix(height))+' km',		$
    /norm, alignment = 0.5
endelse
plotct, ncolors, posc, nf*1e3, cttitle

nf = 10.0^float(round(log(max(sigp)))+1.0)

done = 0
while not done do begin
  loc = where(sigp gt nf,count)
  per = float(count)/nq
  if per gt 5.0/nq then nf = nf*1.25		$
  else if per lt 1.0/nq then nf = nf*0.75		$
  else done = 1
endwhile

image = sigp*float(ncolors)/nf
mltpoly, image, yax, lon, maxran, pos2
plotmlt, maxran, /white

posc = [pos2(2)+0.12,pos2(1),pos2(2)+0.17,pos2(3)]
if hint then 					$
  cttitle = 'Pederson Conductivity (mhos)'	$
else						$
  cttitle = 'Pederson Conductivity (mhos/m)'
plotct, ncolors, posc, nf, cttitle

if !d.name eq 'PS' then begin
  device, /close
  set_plot,'X'
endif

end


