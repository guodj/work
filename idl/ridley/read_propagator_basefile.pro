
pro read_propagator_basefile, file, vars, time, data

  vars = ['x',$
          'y',$
          'z',$
          'vx',$
          'vy',$
          'vz',$
          'lon',$
          'lat',$
          'alt',$
          'sma',$
          'inc',$
          'ecc',$
          'trueano',$
          'raan',$
          'argper',$
          'rightasc',$
          'lt']

  spawn, 'wc '+file, wc
  wc = long(wc(0))-5

  time = dblarr(wc)
  data = fltarr(17,wc)

  openr,1,file
  line = ''
  readf,1,line
  readf,1,line
  while strpos(line,'-----') eq -1 do readf,1,line

  nLines = 0L
  itime = intarr(6)
  while not eof(1) do begin

     readf,1,line
     x = strsplit(line,' ',/extract)
     itime(0:2) = fix(strsplit(x(0),'/',/extract))
     itime(3:5) = fix(strsplit(x(1),':',/extract))
     c_a_to_r, itime, rtime
     time(nLines) = rtime
     data(*,nLines) = float(x(2:18))
     nLines++

  endwhile

  close,1

  nLines=nLines-1L
  data = data(*,0:nLines-1L)
  time = time(0:nLines-1L)

end
