
pro iono_read_file, filename, nvars, nlats, nlons, vars, time, data

  line = ""
  nVars = 0
  nLats = 0
  nLons = 0

  openr,1,filename

  done = 0

  while (not done) do begin

      readf,1, line

      if (strpos(mklower(line),"numerical") gt -1) then begin

          readf,1, nvars
          readf,1, nlats
          readf,1, nlons

          tmp = fltarr(nvars)
          vars = strarr(nvars)

      endif

      if (strpos(mklower(line),"variable") gt -1) then begin

          for i=0,nvars-1 do begin
              readf,1,line
              vars(i) = strmid(line,6,strlen(line)-6)
          endfor

      endif

      if (strpos(mklower(line),"time") gt -1 and $
          strpos(mklower(line),"simulation") lt 0) then begin

          itime = intarr(6)
          int_tmp = 0
          for i=0,5 do begin
              readf, 1, int_tmp
              itime(i) = int_tmp
          endfor

          c_a_to_r, itime, time

      endif

      if (strpos(mklower(line),"northern") gt -1) then begin

          data = fltarr(2,nvars,nlons,nlats)
          for j=0,nlons-1 do for i=0,nlats-1 do begin
              readf,1,tmp
              data(0,*,j,i) = tmp
          endfor

      endif

      if (strpos(mklower(line),"southern") gt -1) then begin

          for j=0,nlons-1 do for i=0,nlats-1 do begin
              readf,1,tmp
              data(1,*,j,i) = tmp
          endfor

      endif

      if eof(1) then done = 1
          
  endwhile

  close,1

end

