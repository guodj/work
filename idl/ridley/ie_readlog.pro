
pro ie_readlog, files, time, data, Vars

  nFiles = n_elements(files)

  nLinesTotal = 0L
  nLinesSub = lonarr(nFiles)
  for i=0,nFiles-1 do begin
     spawn, 'wc '+files(i), wc
     wc = long(wc)
     nLinesTotal=nLinesTotal+wc(0)-2L
     nLinesSub(i) = wc(0)-2
  endfor

  time = dblarr(nLinesTotal)

  line = ''
  iLineT = 0L
  for iFile=0L,nFiles-1 do begin

     openr,1,files(iFile)

     readf,1,line
     readf,1,line
     Vars = strsplit(line, ' ', /extract)

     nVars = n_elements(Vars)

     if (iFile eq 0) then data = fltarr(nVars, nLinesTotal)

     if (Vars(1) eq 'yy') then ConvertTime = 1 else ConvertTime = 0
     tmp = fltarr(nVars)
     for iLine = 0L, nLinesSub(iFile)-1 do begin
        readf,1,tmp
        data(*, iLineT) = tmp
        if (ConvertTime) then begin
           itime = tmp(1:6)
           c_a_to_r, itime, rtime
           time(iLineT) = rtime
        endif
        iLineT = iLineT+1L
     endfor

     close,1

  endfor

end
