
pro fpi_readfile, filename, time, data, ntimes, station, vars, codes

  ntimes = intarr(9)
  utoff = 9.0 * 3600.0

  file = findfile(filename)
  if (strlen(file(0)) gt 0) then begin

     openr,1,file
     tmp = fltarr(9)
     data = fltarr(9,8,1000)
     time = dblarr(9,1000)
     codes = strarr(9)

     OldCode = 'blerg'
     station = ''
     date = ''

     readf,1,station
     station = strcompress(mklower(station),/remove_all)
     readf,1,date
     date = strcompress(date)
     date = strmid(date,1,strlen(date))
     if (fix(date) lt 10) then date = '0'+date
     date = strmid(date,0,7)+strmid(date,9,2)
     print, date
     c_s_to_a, itime, date
     c_a_to_r, itime, basetime

     line = ''
     code = ''
     readf,1,line
     iPos = -1

     while (not eof(1)) do begin

        readf,1,tmp,code
        if (strpos(code,oldcode) eq -1) then iPos = iPos + 1

        codes(iPos) = code
        oldcode = code
        data(iPos,*,ntimes(iPos)) = tmp(1:8)
        time(iPos,ntimes(iPos)) = basetime + tmp(0)*3600.0 + utoff

        ntimes(iPos) = ntimes(iPos) + 1

     endwhile

     close,1

     vars = ['Velocity(m/s)','Ap?','Temp (K)','Unknown','Intensity(R)','Unknown','Azimuth(deg)','Elevation(deg)']

  endif

end
