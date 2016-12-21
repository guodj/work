

VarsToGet = [0,1,2,4,5,6,15]
gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                    VarsToGet = VarsToGet

lons = reform(data(0,*,*,*))
lats = reform(data(1,*,*,*))
alts = reform(data(2,*,*,*))
n    = reform(data(3,*,*,*))+reform(data(4,*,*,*))+reform(data(5,*,*,*))
t    = reform(data(6,*,*,*))

p = n*ck_*t

nAlts = n_elements(alts(0,0,*))
p1d = fltarr(nAlts)
for iAlt = 0,nAlts-1 do begin
   p1d(iAlt) = mean(p(*,*,iAlt)*cos(lats(*,*,iAlt)))/mean(cos(lats(*,*,iAlt)))
endfor

end
