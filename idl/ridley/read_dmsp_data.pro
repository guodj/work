
pro read_dmsp_data, DmspRootDir, DmspData, starttime, endtime

  nDays = round((endtime-starttime)/86400.0)

  nLinesMax = 0
  nFilesTotal = 0

  ; The first time through the loop, just figure out how many 
  ; files exist and the maximum number of lines in the files

  print, nDays

  for iDay = 0, nDays-1 do begin

     print, iDay

     t = starttime + iDay*86400.0
     c_r_to_a, itime, t
     c_a_to_ymd, itime, ymd

     dir = DmspRootDir+"/"+tostr(itime(0),4)+"/"+tostr(itime(1),2)+"/"
     print, dir+"f*_"+ymd+"_*.txt"
     filelist = findfile(dir+"f*_"+ymd+"_*.txt")

     nFiles = n_elements(filelist)
     if (nFiles gt 1) then begin

        for iFile = 0, nFiles-1 do begin
           file = filelist(iFile)
           spawn, "wc "+file, wc
           if (fix(wc(0)) gt nLinesMax) then nLinesMax = fix(wc(0))
           print, 'wc '+file+': ',fix(wc(0)), nLinesMax
           nFilesTotal++
        endfor

        print, "nLinesMax : ",nLinesMax

     endif

  endfor

  DmspSubData = {nPts        : 0L, $
                 file        : filelist(0), $
                 satellite   : "f13", $
                 time        : dblarr(nLinesMax), $
                 quality     : intarr(2,nLinesMax), $
                 GeoPos      : fltarr(3,nLinesMax), $
                 MagPos      : fltarr(2,nLinesMax), $
                 Velocity    : fltarr(3,nLinesMax), $
                 Sigma       : fltarr(3,nLinesMax), $
                 NDensity    : fltarr(nLinesMax), $
                 Frac        : fltarr(3,nLinesMax), $
                 Temperature : fltarr(2,nLinesMax)}

  DmspData = replicate(DmspSubData,nFilesTotal)

  iFileTotal = 0

  for iDay = 0, nDays-1 do begin

     t = starttime + iDay*86400.0
     c_r_to_a, itime, t
     c_a_to_ymd, itime, ymd

     dir = DmspRootDir+"/"+tostr(itime(0),4)+"/"+tostr(itime(1),2)+"/"
     filelist = findfile(dir+"f*_"+ymd+"_*.txt")

     nFiles = n_elements(filelist)
     if (nFiles gt 1) then begin

        for iFile = 0, nFiles-1 do begin
           file = filelist(iFile)

           print, 'reading file ',file
           read_dmsp_webfile, file, satellite, nPts, $
                              time, Quality, GeoPos, MagPos, Velocity, $
                              Sigma, NDensity, Frac, Temperature, Points
 
           DmspData(iFileTotal).nPts = nPts
           DmspData(iFileTotal).file = file
           DmspData(iFileTotal).satellite = satellite
           DmspData(iFileTotal).time(0:nPts-1) = time(0:nPts-1)
           DmspData(iFileTotal).quality(*,0:nPts-1) = quality(*,0:nPts-1)
           DmspData(iFileTotal).GeoPos(*,0:nPts-1) = GeoPos(*,0:nPts-1)
           DmspData(iFileTotal).MagPos(*,0:nPts-1) = MagPos(*,0:nPts-1)
           DmspData(iFileTotal).Velocity(*,0:nPts-1) = Velocity(*,0:nPts-1)
           DmspData(iFileTotal).Sigma(*,0:nPts-1) = Sigma(*,0:nPts-1)
           DmspData(iFileTotal).NDensity(0:nPts-1) = NDensity(0:nPts-1)
           DmspData(iFileTotal).Frac(*,0:nPts-1) = Frac(*,0:nPts-1)
           DmspData(iFileTotal).Temperature(*,0:nPts-1) = Temperature(*,0:nPts-1)
           iFileTotal++
        endfor

     endif

  endfor

end
