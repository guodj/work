
pro get_gitm_points, Times, Lons, Lats, Alts, VarsToGet, Data

  nVarsToGet = n_elements(VarsToGet)

  print, "Getting GitmTimes..."
  filelist = findfile('3DALL*bin')
  get_gitm_times, filelist, gitmtimes

  nTimes = n_elements(Times)
  iFiles = intarr(nTimes)

  print, "Finding which data goes with which GITM file..."
  for i=0L,nTimes-1 do begin

     l = where(gitmtimes lt Times(i), c)
     if (c gt 0) then iFiles(i) = c-1 else iFiles(i) = 0

  endfor

  iFilesPlus = iFiles+1

  iFilesToRead = [iFiles,iFilesPlus]

  nFilesEstimate = max(iFilesToRead)-min(iFilesToRead)+1
  iFilesUnique = intarr(nFilesEstimate)

  print, "Determining unique GITM files..."
  iFile = 0
  while (max(iFilesToRead) gt 0) do begin
     m = max(iFilesToRead)
     l = where(iFilesToRead eq m)
     iFilesUnique(iFile) = m
     iFile++
     iFilesToRead(l) = -1
  endwhile

  nFiles = iFile
  iFilesUnique = reverse(iFilesUnique(0:nFiles-1))

  print, "Reading GITM files and interpolating to points..."
  for iFile = 0, nFiles-2 do begin

     t1 = gitmtimes(iFilesUnique(iFile))
     t2 = gitmtimes(iFilesUnique(iFile)+1)
     l = where(times ge t1 and times lt t2, c)

     if (c gt 0) then begin

        if (iFile eq 0) then begin
           AddedVars = [0,1,2]
           file = filelist(iFilesUnique(iFile))
           print, "Reading File : ",file
           gitm_read_bin_1var, file, gitmcoords, gitmtime1, nVars, Vars, ver,$
                               VarsToGet = AddedVars
           gitmLons = reform(gitmcoords(0,*,0,0))/!dtor
           dLon = gitmLons(1)-gitmLons(0)
           gitmLats = reform(gitmcoords(1,0,*,0))/!dtor
           dLat = gitmLats(1)-gitmLats(0)
           gitmAlts = reform(gitmcoords(2,0,0,*))/1000.0

           gitm_read_bin_1var, file, gitmdata1, gitmtime1, nVars, Vars, ver,$
                               VarsToGet = VarsToGet

           iLon = (Lons-gitmLons(0))/dlon
           xLon = iLon - fix(iLon)
           iLon = fix(iLon)
           iLat = (Lats-gitmLats(0))/dlat
           xLat = iLat - fix(iLat)
           iLat = fix(iLat)

           nAlts = n_elements(gitmAlts)
           Data = fltarr(nTimes,nVarsToGet,nAlts)

        endif

        file = filelist(iFilesUnique(iFile)+1)
        print, "Reading File : ",file
        gitm_read_bin_1var, file, gitmdata2, gitmtime2, nVars, Vars, ver, $
                            VarsToGet = VarsToGet

        for iD = 0, c-1 do begin
           iT = l(iD)
           x = (times(iT) - gitmtime1)/(gitmtime2-gitmtime1)
           t = (1.0-x) * gitmtime1 + x * gitmtime2
           if (gitmtime2-gitmtime1 eq 0) then stop
           
           for iVar = 0,nVarsToGet-1 do begin
              d1 = $ 
                 (1.0-xLon(iT)) * (1.0-xLat(iT)) * gitmdata1(iVar,iLon(iT),  iLat(iT)  ,*) + $
                 (    xLon(iT)) * (1.0-xLat(iT)) * gitmdata1(iVar,iLon(iT)+1,iLat(iT)  ,*) + $
                 (1.0-xLon(iT)) * (    xLat(iT)) * gitmdata1(iVar,iLon(iT),  iLat(iT)+1,*) + $
                 (    xLon(iT)) * (    xLat(iT)) * gitmdata1(iVar,iLon(iT)+1,iLat(iT)+1,*)
              d2 = $ 
                 (1.0-xLon(iT)) * (1.0-xLat(iT)) * gitmdata2(iVar,iLon(iT),  iLat(iT)  ,*) + $
                 (    xLon(iT)) * (1.0-xLat(iT)) * gitmdata2(iVar,iLon(iT)+1,iLat(iT)  ,*) + $
                 (1.0-xLon(iT)) * (    xLat(iT)) * gitmdata2(iVar,iLon(iT),  iLat(iT)+1,*) + $
                 (    xLon(iT)) * (    xLat(iT)) * gitmdata2(iVar,iLon(iT)+1,iLat(iT)+1,*)
              Data(iT,iVar,*) = (1.0-x) * d1 + x * d2
           endfor

        endfor

     endif

     gitmtime1 = gitmtime2
     gitmdata1 = gitmdata2

  endfor

  if (max(Alts) gt 0.0) then begin
     DataNew = fltarr(nTimes,nVarsToGet)
     for iTime = 0, nTimes-1 do begin
        l=where(gitmAlts gt Alts(iTime),c)
        if (c eq 0) then begin
           DataNew(iTime,*) = Data(iTime,*,nAlts-1)
        endif else begin
           iAlt = l(0)-1
           if (iAlt le 1) then begin
              DataNew(iTime,*) = Data(iTime,*,iAlt)
           endif else begin
              x = (Alts(iTime)-gitmAlts(iAlt))/(gitmAlts(iAlt+1) - gitmAlts(iAlt))
              DataNew(iTime,*) = (1.0-x)*Data(iTime,*,iAlt) + $
                                 (    x)*Data(iTime,*,iAlt+1)
           endelse
        endelse

     endfor

     Data = DataNew

  endif

end
