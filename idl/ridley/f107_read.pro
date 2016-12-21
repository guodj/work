
pro f107_read, time, f107

  file = '/raid3/Data/f107.txt'
  spawn, 'wc '+file,wc
  wc = fix(wc)
  wc = wc(0)

  close,1
  openr,1,file

  line = ''
  while (strpos(line,'yyyy') lt 0) do readf,1,line

  time = dblarr(wc)
  f107 = fltarr(wc)
  i = 0
  while (not eof(1)) do begin

     readf,1,yyyy,mm,dd,hh,mi, f, format='(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.1)'
     itime = [yyyy,mm,dd,hh,mi,0]
     c_a_to_r, itime, rtime
     time(i) = rtime
     f107(i) = f
     i++

  endwhile

  close,1

  time = time(0:i-1)
  f107 = f107(0:i-1)

end

