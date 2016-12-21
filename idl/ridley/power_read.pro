
pro split_line_old, line, time1, power1, IsNorth

  ; sample: 'NOAA-17(N)1340011    13.2  5   1.367'

  doy = float(strmid(line,10,3))
  hr  = float(strmid(line,13,2))
  mi  = float(strmid(line,15,2))
  time1 = doy*86400.0 + hr*3600.0 + mi*60.0
  power1 = float(strmid(line,17,8))
  if (strpos(line,'(N)') gt 0) then IsNorth = 1 else IsNorth = 0

end

pro split_line_new, line, time1, power1, IsNorth

  ; sample: '2011-01-01 00:38:02 NOAA-17 (S)  3   4.80   0.90'

  itime = intarr(6)
  itime(0) = fix(strmid(line,0,4))
  for i=1,5 do itime(i) = fix(strmid(line,2+i*3,2))
  c_a_to_r, itime, time1
  power1 = float(strmid(line,34,7))
  if (strpos(line,'(N)') gt 0) then IsNorth = 1 else IsNorth = 0

end

pro power_read, file, time, power, North, nPts, notes

  spawn,'wc '+file, wc
  nLines = fix(wc)
  nLines = nLines(0)

  time  = dblarr(nLines)
  power = fltarr(nLines)
  North = intarr(nLines)

  close,1
  openr,1,file

  line = ''
  notes = ['']

  while (strpos(line, 'Normalizing factor') lt 0) do begin
     readf,1,line
     notes = [notes,line]
  endwhile

  readf,1,line
  readf,1,line

  itime = intarr(6)

  ; now we see which format the data is in....
  if (strlen(line) lt 5) then IsOldFormat = 1 else IsOldFormat = 0

  nPts = 0

  if IsOldFormat then begin
     itime(0) = fix(line)
     itime(1) = 1
     itime(2) = 1
     itime([3,4,5]) = 0
     c_a_to_r, itime, basetime
  endif else begin
     split_line_new, line, time1, power1, IsNorth
     power(nPts) = power1
     time(nPts) = time1
     North(nPts) = IsNorth
     nPts = nPts + 1
  endelse

  while not eof(1) do begin

     readf,1,line

     if IsOldFormat then begin
        split_line_old, line, time1, power1, IsNorth
        time1 = basetime + double(time1)        
     endif else split_line_new, line, time1, power1, IsNorth

     time(nPts) = time1
     power(nPts) = power1
     North(nPts) = IsNorth
     nPts = nPts + 1

  endwhile

  close,1

  time = time(0:npts-1)
  power = power(0:npts-1)
  north = north(0:npts-1)

end
