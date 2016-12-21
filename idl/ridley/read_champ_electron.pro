
pro read_champ_electron, file, time, data

  f = findfile(file)
  if (strlen(f(0)) le 1) then begin

     time = dblarr(1)
     time(0) = -999.0

  endif else begin

     nF = n_elements(f)
     file = f(nF-1)

     spawn, 'wc '+file,wc
     nLines = fix(wc)
     nLines = nLines(0)

     close,1
     openr,1,file
     line = ''
     nHeader = 0

     while strpos(line,'[s]') lt 0 do begin
        readf,1,line
        nHeader++
     endwhile
     readf,1,line
     nHeader++
     nLines = nLines-nHeader
     data = fltarr(4,nLines)
     time = dblarr(nLines)

     for iLine = 0, nLines-1 do begin

        readf,1,line
        tmp = strsplit(line,' ',/extract)
        itime = fix(tmp(1:6))
        c_a_to_r, itime, rtime
        time(iLine) = rtime
        data(0,iLine) = float(tmp( 7))
        data(1,iLine) = float(tmp( 8))
        data(2,iLine) = float(tmp( 9))
        data(3,iLine) = float(tmp(10))*1.0e6

     endfor

  endelse

end

