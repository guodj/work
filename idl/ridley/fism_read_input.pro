
pro fism_read_input, file, time, data

  close,1

  openr,1,file

  nBins = 59L
  nTimesMax = 1440L*32L

  nTimes = 0L

  line = ''
  while strpos(line,'#START') lt 0 do readf,1,line

  tmp = fltarr(nBins)
  iTime = intarr(6)
  data = fltarr(nBins,nTimesMax)
  time = dblarr(nTimesMax)

  while not eof(1) do begin

     readf,1,iTime,tmp
     c_a_to_r, iTime, rTime
     time(nTimes) = rTime
     data(*,nTimes) = tmp
     
     nTimes = nTimes + 1
     
  endwhile

  time = time(0:nTimes-1)
  data = data(*,0:nTimes-1)

  close,1

end
