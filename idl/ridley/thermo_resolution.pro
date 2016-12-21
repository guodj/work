
filelist = findfile("-t *.bin")
filename = filelist(0)

read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
                        vars, data, rb, cb, bl_cnt, iTime, Version

lats = reform(data(1,0,*,0)) / !dtor
nLats = n_elements(lats)

dlat = fltarr(nLats)
dlat(1:nLats-2) = (lats(2:nLats-1)-lats(0:nLats-3))/2
dlat(0) = lats(1)-lats(0)
dLat(nLats-1) = Lats(nLats-1)-Lats(nLats-2)

plot, lats, dlat

end
