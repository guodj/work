
;--------------------------------------------------------------
; Get Inputs from the user
;--------------------------------------------------------------

amie_file = ask('AMIE binary file name','b910603.all')
psfile = ask('ps file',amie_file+'.ps')

read_amie_binary, amie_file, data, lats, mlts, time, fields, 		$
                  imf, ae, dst, hp, cpcp

nped = 2
nhal = 4
nfac = 7

fac = reform(data(*,nfac,*,*))
hal = reform(data(*,nhal,*,*))
ped = reform(data(*,nped,*,*))

ntimes = n_elements(fac(*,0,0))
nmlts = n_elements(fac(0,*,0))
nlats = n_elements(fac(0,0,*))

peak_hal = fltarr(ntimes,nmlts)
peak_ped = fltarr(ntimes,nmlts)
peak_fac = fltarr(ntimes,nmlts)

for i=0,ntimes-1 do for j=0,nmlts-1 do begin
  peak_hal(i,j) = max(hal(i,j,*))
  peak_ped(i,j) = max(ped(i,j,*))
  peak_fac(i,j) = max(abs(fac(i,j,*)))
endfor

end




