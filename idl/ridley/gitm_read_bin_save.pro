
pro gitm_read_bin, file, data, time, nVars, Vars, version, skip = skip

  if (n_elements(skip) eq 0) then skip = 0
  
  if (n_elements(file) eq 1) then filelist = findfile(file) $
  else filelist = file

  nFiles = n_elements(filelist)

  if (nFiles gt 1) then Time = dblarr(nFiles) else Time = 0.0D

  DoAppendFile = 0

  OldTime = -1.0d32

  for iFile = 0L, nFiles-1 do begin

     filein = filelist(iFile)

     valid = 0
     on_ioerror, skipfile

     close, 1
     openr, 1, filein, /f77

     version = 0.0D

     nLons = 0L
     nLats = 0L
     nAlts = 0L
     nVars = 0L
     
     nTimes = 0L

     while (not eof(1)) do begin

        readu, 1, version
        readu, 1, nLons, nLats, nAlts
        readu, 1, nVars

        Vars = strarr(nVars)
        line = bytarr(40)
        for iVars = 0, nVars-1 do begin
           readu, 1, line
           Vars(iVars) = strcompress(string(line),/remove)
        endfor

        lTime = lonarr(7)
        readu, 1, lTime

        iTime = fix(lTime(0:5))
        c_a_to_r, itime, rtime
        Time(iFile) = rTime + lTime(6)/1000.0

        if (nFiles eq 1) then begin
           Data = dblarr(nVars, nLons, nLats, nAlts)
        endif else begin
           if (iFile eq 0) then $
              Data = dblarr(nFiles, nVars, nLons, nLats, nAlts)
        endelse

        tmp = dblarr(nLons, nLats, nAlts)
        for i=0,nVars-1 do begin
           readu,1,tmp
           if (nFiles eq 1) then data(i,*,*,*) = tmp $
           else data(iFile,i,*,*,*) = tmp
        endfor

        if (not eof(1)) then DoAppendFile = 1

        if (DoAppendFile and nFiles eq 1) then begin

           if (nTimes eq 0) then begin
              spawn, 'ls -s '+filein,filesize
              filesize = long(filesize(0))*1024.0
              datasize = (nVars*nLons*nLats*nAlts*8L + $
                          8 + 3*4 + 4 + nVars*40 + 7*4) + $
                         2*4*(3+nVars+1+nVars)
              nTimesMax = filesize/datasize
              print, 'Assuming that file contains nTimes : ',nTimesMax
              DataAllNew = dblarr(nTimesMax, nVars, nLons, nLats, nAlts)
              TimeAllNew = dblarr(nTimesMax)
           endif
           if (nTimes lt nTimesMax) then begin
              DataAllNew(nTimes,*,*,*,*) = Data
              TimeAllNew(nTimes) = Time
              nTimes = nTimes+1L
           endif else begin
              print, "Estimate of nTimesMax is wrong!"
              print, "nTimesMax, nTimes : ", nTimesMax, nTimes
           endelse

        endif

     endwhile

     valid = 1

     skipfile: 
     if (not valid) then begin
        print, "Bad file found : ",filein
        if (nFiles gt 1) then time(ifile) = 0.0
     endif

     close, 1

  endfor

  if (nTimes gt 1) then begin
     data = DataAllNew
     time = TimeAllNew
  endif

  ; If we had some bad files, we could be running into trouble with
  ; missing data.  Clean it up here.

  if (nTimes gt 1 or nFiles gt 1) then begin
     nTimes = n_elements(time)
     l = where(time gt 0.0,c)
     if (c lt nTimes) then begin
        print, "Found some empty times.  Correcting."
        DataAllNew = dblarr(c, nVars, nLons, nLats, nAlts)
        TimeAllNew = dblarr(c)
        for i=0L,c-1 do begin
           DataAllNew(i,*,*,*,*) = data(l(i),*,*,*,*)
           TimeAllNew(i) = time(l(i))
        endfor
        data = DataAllNew
        time = TimeAllNew
     endif
  endif



end

