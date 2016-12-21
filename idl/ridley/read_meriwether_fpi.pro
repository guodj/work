
pro read_meriwether_fpi, file, time, vz, vze, tn, tne, in

  utoff = 9.0 * 3600.0

  close,1

  openr,1,file

  nPts = 0

  readf,1,nPts

  time = dblarr(nPts)
  vz = fltarr(nPts)
  vze = fltarr(nPts)
  tn = fltarr(nPts)
  tne = fltarr(nPts)
  in = fltarr(nPts)

  l = ''
  tmp = fltarr(9)

  for i=0,nPts-1 do begin
     readf,1,tmp,l
     time(i) = tmp(0)*3600.0 + utoff
     vz(i)   = tmp(1)
     vze(i)  = tmp(2)
     tn(i)   = tmp(3)
     tne(i)  = tmp(4)
     in(i)   = tmp(5)
  endfor

  close,1

end
