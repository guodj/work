
pro read_thermosphere_file, filelist, nvars, nalts, nlats, nlons, $
                            vars, data, nBLKlat, nBLKlon, nBLK, $
                            iTime, Version

  if (n_elements(nBLKlat) eq 0) then nBLKlat = 0
  if (n_elements(nBLKlon) eq 0) then nBLKlon = 0
  if (n_elements(nBLK)    eq 0) then nBLK    = 0

  filelist = findfile(filelist)

  Version = -1.0

  if (strpos(filelist(0), "save") gt 0) then begin
      restore, filelist(0)
      nBLK = 1
      if (n_elements(iTime) eq 0) then begin
          p = strpos(filelist(0),".save")-1
          while (strpos(strmid(filelist(0),p,1),'.') eq -1) do p = p-1
          iYear   = fix(strmid(filelist(0),p-13,2))
          iMonth  = fix(strmid(filelist(0),p-11,2))
          iDay    = fix(strmid(filelist(0),p-9,2))
          iHour   = fix(strmid(filelist(0),p-6,2))
          iMinute = fix(strmid(filelist(0),p-4,2))
          iSecond = fix(strmid(filelist(0),p-2,2))
          iTime = [iYear, iMonth, iDay, iHour, iMinute, iSecond]
      endif
      return
  endif else begin
      if (strpos(filelist(0), "bin") gt 0) then begin
          gitm_read_bin, filelist, data, time, nVars, Vars, version
          s = size(data)
          if (s(0) eq 5) then begin
              nLons = s(3)
              nLats = s(4)
              nAlts = s(5)
          endif
          if (s(0) eq 4) then begin
              nLons = s(2)
              nLats = s(3)
              nAlts = s(4)
              if (nLons eq 24 and nLats eq 24) then begin
                  if max(data(2,22:23,22:23,0)) eq 0.0 then begin
                      print, "Error in file.  Correcting"
                      data = data(*,0:21,0:21,*)
                      nLons = 22
                      nLats = 22
                  endif
              endif
          endif
          if (s(0) eq 3) then begin
              nLons = s(2)
              nLats = s(3)
              nAlts = 1
          endif
          nBLK = 1
          nTimes = n_elements(time)
          if (nTimes eq 1) then $
             c_r_to_a, itime, time(0) $
          else begin
             itime = intarr(nTimes,6)
             for iT=0,nTimes-1 do begin
                c_r_to_a, itimeSingle, time(iT)
                itime(iT,*) = itimeSingle
             endfor
          endelse

          return
      endif
  endelse

  f = filelist(0)

  if (strpos(filelist(0),"b0") eq 0) then begin
      all = findfile('b0*'+strmid(filelist(0),6,18)+'*')
  endif else begin
      all = filelist(0)
  endelse

  nfiles = n_elements(all)

  if (nBLKlat eq 0 and nBLKlon eq 0) then begin

    if (nfiles eq 1) then begin
      nBLKlon = 1
      nBLKlat = 1
    endif

 endif

end

