function read_izmem, filename

  openr,1, filename

  line = ''

  readf,1,line

  imfx = strmid(line,0,9)
  imfy = strmid(line,12,9)
  imfz = strmid(line,24,9)
  timedum = strmid(line,92,3)+strmid(line,96,2)
  if strmid(timedum,3,1) eq ' ' then 			$
    timedum = strmid(timedum,0,3)+'0'+strmid(timedum,4,1)
  datedum = strmid(line,78,10)

  datadum = fltarr(24,89-57+1)

  for i=89,57,-1 do begin

   readf,1,line

   for j=0,23 do datadum(j,89-i) = float(strmid(line,j*10+14,10))

  endfor

  close,1

  data = {pot:datadum, imf:[imfx,imfy,imfz],time:timedum, date:datedum}

  return, data

end
