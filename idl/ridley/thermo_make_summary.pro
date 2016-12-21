
iRho  = 3
iTemp = 15
iVV   = 18
iAlt_  = 2

if (n_elements(dir) eq 0) then dir = 'data.01'
dir = ask('directory',dir)

filelist = findfile(dir+'/3DALL*.bin')

nFiles = n_elements(filelist)

for iFile = 0, nFiles-1 do begin

    file = filelist(iFile)
    print, "Reading File : ",file
    gitm_read_bin, file, data, time, nVars, Vars, version

    if (iFile eq 0) then begin

       lats = reform(data(1,*,*,0))
       lons = reform(data(0,*,*,0))
       nLons = n_elements(lons(*,0))
       nLats = n_elements(lats(0,*))
       dlon = lons
       dlon(0:nlons-2,*) = lons(1:nLons-1,*) - lons(0:nLons-2,*)
       dlon(nlons-1,*) = dlon(0,*)

       dlat = lats
       dlat(*,1:nlats-2) = (lats(*,2:nLats-1) - lats(*,0:nLats-3))/2.0
       dlat(*,0) = lats(*,1) - lats(*,0)
       dlat(*,nLats-1) = lats(*,nLats-1) - lats(*,nLats-2)
       area = 6372000.0*6372000.0 * cos(lats) * dlon * dlat

       nAlts = n_elements(data(0,0,0,*))

       global_rho  = fltarr(4, nFiles, nAlts-4)
       global_temp = fltarr(4, nFiles, nAlts-4)
       global_vv   = fltarr(4, nFiles, nAlts-4)
       global_alt  = fltarr(4, nFiles, nAlts-4)

       hem_temp    = fltarr(4, 2, nFiles, nAlts-4)
       hem_rho     = fltarr(4, 2, nFiles, nAlts-4)
       hem_vv      = fltarr(4, 2, nFiles, nAlts-4)
       hem_alt     = fltarr(4, 2, nFiles, nAlts-4)

       time_all = dblarr(nFiles)

    endif

    time_all(iFile) = time

    area_global = area(2:nLons-3,2:nLats-3)
    mag = mean(area_global)

    for iAlt = 2, nAlts-3 do begin

        global_rho(0, iFile, iAlt-2) = $
          mean(data(iRho, 2:nLons-3, 2:nLats-3, iAlt)*area_global)/mag
        global_temp(0, iFile, iAlt-2) = $
          mean(data(iTemp, 2:nLons-3, 2:nLats-3, iAlt)*area_global)/mag
        global_vv(0, iFile, iAlt-2) = $
          mean(data(iVV, 2:nLons-3, 2:nLats-3, iAlt)*area_global)/mag
        global_alt(0, iFile, iAlt-2) = $
          mean(data(iAlt_, 2:nLons-3, 2:nLats-3, iAlt)*area_global)/mag

        global_rho(1, iFile, iAlt-2) = $
          min(data(iRho, 2:nLons-3, 2:nLats-3, iAlt))
        global_temp(1, iFile, iAlt-2) = $
          min(data(iTemp, 2:nLons-3, 2:nLats-3, iAlt))
        global_vv(1, iFile, iAlt-2) = $
          min(data(iVV, 2:nLons-3, 2:nLats-3, iAlt))
        global_alt(1, iFile, iAlt-2) = $
          min(data(iAlt_, 2:nLons-3, 2:nLats-3, iAlt))

        global_rho(2, iFile, iAlt-2) = $
          max(data(iRho, 2:nLons-3, 2:nLats-3, iAlt))
        global_temp(2, iFile, iAlt-2) = $
          max(data(iTemp, 2:nLons-3, 2:nLats-3, iAlt))
        global_vv(2, iFile, iAlt-2) = $
          max(data(iVV, 2:nLons-3, 2:nLats-3, iAlt))
        global_alt(2, iFile, iAlt-2) = $
          max(data(iAlt_, 2:nLons-3, 2:nLats-3, iAlt))

        global_rho(3, iFile, iAlt-2) = $
          stddev(data(iRho, 2:nLons-3, 2:nLats-3, iAlt))
        global_temp(3, iFile, iAlt-2) = $
          stddev(data(iTemp, 2:nLons-3, 2:nLats-3, iAlt))
        global_vv(3, iFile, iAlt-2) = $
          stddev(data(iVV, 2:nLons-3, 2:nLats-3, iAlt))
        global_alt(3, iFile, iAlt-2) = $
          stddev(data(iAlt_, 2:nLons-3, 2:nLats-3, iAlt))

    endfor

    iLs = 2
    iLe = (nLats-4)/4+2-1

    print, 'south : ',iLs, iLe

    area_south = area(2:nLons-3,iLs:iLe)
    mas = mean(area_south)

    ; South -----------------
    for iAlt = 2, nAlts-3 do begin
        hem_rho(0,0,iFile, iAlt-2) = $
          mean(data(iRho, 2:nLons-3, iLs:iLe, iAlt)*area_south)/mas
        hem_temp(0,0,iFile, iAlt-2) = $
          mean(data(iTemp, 2:nLons-3, iLs:iLe, iAlt)*area_south)/mas
        hem_vv(0,0,iFile, iAlt-2) = $
          mean(data(iVV, 2:nLons-3, iLs:iLe, iAlt)*area_south)/mas
        hem_alt(0,0,iFile, iAlt-2) = $
          mean(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt)*area_south)/mas

        hem_rho(1,0,iFile, iAlt-2) = $
          min(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(1,0,iFile, iAlt-2) = $
          min(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(1,0,iFile, iAlt-2) = $
          min(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(1,0,iFile, iAlt-2) = $
          min(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

        hem_rho(2,0,iFile, iAlt-2) = $
          max(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(2,0,iFile, iAlt-2) = $
          max(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(2,0,iFile, iAlt-2) = $
          max(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(2,0,iFile, iAlt-2) = $
          max(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

        hem_rho(3,0,iFile, iAlt-2) = $
          stddev(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(3,0,iFile, iAlt-2) = $
          stddev(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(3,0,iFile, iAlt-2) = $
          stddev(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(3,0,iFile, iAlt-2) = $
          stddev(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))
    endfor

    iLs = 3*(nLats-4)/4+2
    iLe = nLats-3

    print, 'north : ',iLs, iLe

    area_north = area(2:nLons-3,iLs:iLe)
    man = mean(area_north)

    ; North -----------------
    for iAlt = 2, nAlts-3 do begin

        hem_rho(0,1,iFile, iAlt-2) = $
          mean(data(iRho, 2:nLons-3, iLs:iLe, iAlt)*area_north)/man
        hem_temp(0,1,iFile, iAlt-2) = $
          mean(data(iTemp, 2:nLons-3, iLs:iLe, iAlt)*area_north)/man
        hem_vv(0,1,iFile, iAlt-2) = $
          mean(data(iVV, 2:nLons-3, iLs:iLe, iAlt)*area_north)/man
        hem_alt(0,1,iFile, iAlt-2) = $
          mean(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt)*area_north)/man

        hem_rho(1,1,iFile, iAlt-2) = $
          min(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(1,1,iFile, iAlt-2) = $
          min(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(1,1,iFile, iAlt-2) = $
          min(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(1,1,iFile, iAlt-2) = $
          min(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

        hem_rho(2,1,iFile, iAlt-2) = $
          max(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(2,1,iFile, iAlt-2) = $
          max(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(2,1,iFile, iAlt-2) = $
          max(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(2,1,iFile, iAlt-2) = $
          max(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

        hem_rho(3,1,iFile, iAlt-2) = $
          stddev(data(iRho, 2:nLons-3, iLs:iLe, iAlt))
        hem_temp(3,1,iFile, iAlt-2) = $
          stddev(data(iTemp, 2:nLons-3, iLs:iLe, iAlt))
        hem_vv(3,1,iFile, iAlt-2) = $
          stddev(data(iVV, 2:nLons-3, iLs:iLe, iAlt))
        hem_alt(3,1,iFile, iAlt-2) = $
          stddev(data(iAlt_, 2:nLons-3, iLs:iLe, iAlt))

    endfor

endfor

save, file = dir+"/summary.save", hem_rho, hem_temp, hem_vv, hem_alt, $
  global_rho, global_temp, global_vv, global_alt, time_all

end

