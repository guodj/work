
pro read_champ_winds, file, time, data

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
     readf,1,line
     readf,1,line
     nHeader = 2
     nLines = nLines-nHeader
     data = fltarr(6,nLines)
     time = dblarr(nLines)

     itime = lonarr(6)
     for iLine = 0, nLines-1 do begin

        readf,1,line
        tmp = strsplit(line,' ',/extract)
        itime(0) = fix(tmp(0))
        itime(1) = 1
        itime(2) = fix(tmp(1))
        itime(3) = 0
        itime(4) = 0
        itime(5) = long(tmp(2))

        c_a_to_r, itime, rtime
        time(iLine) = rtime
        data(0,iLine) = float(tmp( 6))
        data(1,iLine) = float(tmp( 8))
        data(2,iLine) = float(tmp( 9))
        data(3,iLine) = float(tmp(10))
        data(4,iLine) = float(tmp(11))
        data(5,iLine) = float(tmp( 5))

     endfor

  endelse

end

