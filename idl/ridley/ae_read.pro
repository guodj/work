
pro ae_read, file, time, ae, starttime = starttime, endtime = endtime

  spawn, 'wc '+file, wc
  nPts = long(wc(0))-1

  ae   = fltarr(nPts, 8)
  time = dblarr(nPts)

  close,1
  openr,1,file

  itime = intarr(6)
  tmp   = fltarr(8)

  for i=0L, nPts-1 do begin

     readf,1,itime,tmp
     c_a_to_r, itime, rtime
     time(i) = rtime
     ae(i,*) = tmp

  endfor

  if (n_elements(starttime) gt 0) then begin
     l = where(time ge starttime,c)
     if (c gt 0) then begin
        time = time(l)
        ae = ae(l,*)
     endif
  endif

  if (n_elements(endtime) gt 0) then begin
     l = where(time le endtime,c)
     if (c gt 0) then begin
        time = time(l)
        ae = ae(l,*)
     endif
  endif

  close, 1

end
