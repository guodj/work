
pro parse_line, line, subdst, subtime

;  print, line

  basetime = subtime
  c_r_to_a, itime, subtime
  itime(2) = fix(line)

  subdst  = fltarr(24)
  subtime = dblarr(24)

;  itime = [year, month, day, 0, 0, 0]

  j = 3
  for i=0,23 do begin

      subdst(i) = float(strmid(line,j,4))
      itime(3) = i
      c_a_to_r, itime, rtime
      subtime(i) = rtime

      j = j + 4
      if (i eq  7) then j = j + 1
      if (i eq 15) then j = j + 1

  endfor

end

pro read_dst, file, dst, time

  openr,1,file

  line = ''
  while (strpos(line,'hourly') lt 0 and not eof(1)) do begin
     readf,1,line
     line = mklower(line)
;     print,line
  endwhile
     
  realtime = 0
  if (strpos(mklower(line),'real-time') gt -1) then realtime = 1
  if (strpos(mklower(line),'provisional') gt -1) then realtime = 1

  readf,1,line
;  if (not realtime) then readf,1,line
  dateline = line

;  if (realtime) then begin
     line = strcompress(dateline)
     l = strlen(line)
     year = fix(strmid(line,l-4,4))
;  endif else year = fix(strmid(dateline,10,4))

  if (strpos(dateline,'JAN') gt -1) then month = 1
  if (strpos(dateline,'FEB') gt -1) then month = 2
  if (strpos(dateline,'MAR') gt -1) then month = 3
  if (strpos(dateline,'APR') gt -1) then month = 4
  if (strpos(dateline,'MAY') gt -1) then month = 5
  if (strpos(dateline,'JUN') gt -1) then month = 6
  if (strpos(dateline,'JUL') gt -1) then month = 7
  if (strpos(dateline,'AUG') gt -1) then month = 8
  if (strpos(dateline,'SEP') gt -1) then month = 9
  if (strpos(dateline,'OCT') gt -1) then month = 10
  if (strpos(dateline,'NOV') gt -1) then month = 11
  if (strpos(dateline,'DEC') gt -1) then month = 12

  if (not realtime) then begin
     while (strpos(line,'DAY') ne 0 and $
            not eof(1)) do readf,1,line
;     print, line
  endif else begin
     readf,1,line
     readf,1,line
  endelse

  dst  = fltarr(31*24)
  time = dblarr(31*24)

  iDay = 0
;  str_to_match = 'CODE'
;  if (realtime) then 
  str_to_match = 'vvvv'
  while (strpos(line,str_to_match) lt 0 and $
         not eof(1)) do begin
      readf,1,line
      if (strpos(line,str_to_match) lt 0) then begin
          if (fix(line) gt 0) then begin
              itime = [year,month,1,0,0,0]
              c_a_to_r, itime, subtime
              parse_line, line, subdst, subtime
              dst(iDay*24:iDay*24+23)  = subdst
              time(iDay*24:iDay*24+23) = subtime
              iDay = iDay + 1
          endif
      endif
  endwhile

  dst  = dst(0:iDay*24-1)
  time = time(0:iDay*24-1)

  close,1

end

