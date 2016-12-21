
pro get_gitm_times, filelist, times

  nFiles = n_elements(filelist)
  times = dblarr(nFiles)
  lTime = lonarr(7)

  for iFile = 0L, nFiles-1 do begin

     if ((iFile+1) mod 100 eq 0) then print, "Progress : ",iFile+1, nFiles

     filein = filelist(iFile)

     close, 1
     openr, 1, filein, /f77

     if (iFile eq 0) then begin

        version = 0.0D

        nLons = 0L
        nLats = 0L
        nAlts = 0L
        nVars = 0L
     
        IsFirstTime = 1

        DoAppendFile = 0

        readu, 1, version
        readu, 1, nLons, nLats, nAlts
        readu, 1, nVars

        line = bytarr(40)
        for iVars = 0, nVars-1 do begin
           readu, 1, line
        endfor

        point_lun, -1, position

     endif else begin
        
        point_lun, 1, position

     endelse

     readu, 1, lTime

     iTime = fix(lTime(0:5))
     c_a_to_r, itime, rtime
     Time = rTime + lTime(6)/1000.0

     times(iFile) = time

     close,1

  endfor


end
