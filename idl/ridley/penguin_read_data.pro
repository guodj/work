
pro penguin_read_data, filein, vars, data, time

  file = findfile(filein)
  if (strlen(file(0)) eq 0) then begin
      print, "File ",filein," not found!!!"
      stop
  endif

  nFiles = n_elements(file)

  spawn, "wc "+filein,list
  nPts = 0L
  for i=0,nFiles-1 do nPts = nPts + long(list(i))-1

  print, nPts, ' expected'

  iPt = 0L
  close,1
  for i=0,nFiles-1 do begin

      line = ''
      print, "Reading File : ", file(i)
      openr,1,file(i)
      readf,1,line

      if (i eq 0) then begin

          l = strsplit(line,',')
          nVars = n_elements(l)
          Vars = strarr(nVars)
          for j=0,nVars-2 do begin
              Vars(j) = strmid(line,l(j),l(j+1)-l(j)-1)
          endfor
          Vars(j) = strmid(line,l(j),strlen(line)-l(j))

          tmp = fltarr(nVars)

          Time = dblarr(nPts)
          Data = fltarr(nVars,nPts)

      endif

      while not eof(1) do begin

          readf,1,tmp
          itime = fix(tmp(1:7))
          c_a_to_r, itime, rtime
          time(iPt) = rtime
          data(*,iPt) = tmp
          iPt = iPt + 1L

      endwhile

      close,1

  endfor

  print, iPt, ' read'
  time = time(0:iPt-1)
  data = data(*,0:iPt-1)

end
