
pro ccmc_results, filename, time, lat, lon, data, shift, ndata

  spawn, 'wc '+filename, wc
  wc = wc(0)
  nPts = long(wc)-shift

  print, filename, ': ',nPts

  if (nPts le 0) then nPts = 1

  time = dblarr(nPts)
  lat  = fltarr(nPts)
  lon  = fltarr(nPts)
  alt  = fltarr(nPts)
  data = fltarr(nPts,ndata)

  if (nPts eq 1) then return

  close,1
  openr,1,filename

  line = ''
  for i=0,shift-1 do readf,1,line

  it = intarr(5)

  nDa = 3
  da = fltarr(3+ndata)

  if (strpos(filename,'hmf2') gt 0) then begin
     da = fltarr(2+ndata+2)
     nDa = 2
  endif
  if (strpos(filename,'nmf2') gt 0 or strpos(filename,'mF2') gt 0) then begin
     da = fltarr(2+ndata)
     nDa = 2
  endif

  for i=0,nPts-1 do begin
     readf,1,it,da
     itime = [it(0), 1, it(1), it(2), it(3), it(4)]
     c_a_to_r, itime, rtime
     time(i) = rtime
     lat(i)  = da(0)
     lon(i)  = da(1)
     if (nDa gt 2) then alt(i)  = da(2)
     for j=0,nData-1 do begin
        data(i,j) = da(nDa+j)
     endfor
  endfor

  close,1

end
