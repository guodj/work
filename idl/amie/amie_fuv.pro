
filelist = findfile("map_energy_*.idl")
nFiles = n_elements(filelist)

MinLat = 40.0

for iFile=0,nFiles-1 do begin

    file = filelist(iFile)
    print, file

    if (iFile eq 0) then begin
        date = strmid(file, 11,7)
        openw,1,"fuv."+date+".amiein"
        printf,1,"Data in Geographic Coordinates"
    endif

    if (strpos(file, date) ne 11) then begin
        date = strmid(file, 11,7)
        close,1
        openw,1,"fuv."+date+".amiein"
        printf,1,"Data in Geographic Coordinates"
    endif

    year = fix(strmid(file, 11,4))
    doy  = fix(strmid(file, 15,3))
    hour = fix(strmid(file, 18,2))
    mini = fix(strmid(file, 20,2))

    itime = [year,1,doy, hour, mini, 0]
    c_a_to_r, itime, rtime
    c_r_to_a, itime, rtime

    restore, file

    nLons = n_elements(final_eflux(*,0))
    nLats = n_elements(final_eflux(0,*))
    dlat  = final_lat
    dlon  = final_lon

    dlon_amie = 10.0
    nlons_amie = 360.0 / dlon_amie + 1

    lats = dlat * findgen(nLats)
    lats = lats + (89.0 - lats(nLats-1))
    lons = dlon * findgen(nLons)

    lons_amie = dlon_amie * findgen(nlons_amie)

    loc = where(lats gt MinLat, count)

    for iiLat = 0, count-1 do begin

        iLat = loc(iiLat)

        for i= 0, nlons_amie-1 do begin

            locl = where(lons ge lons_amie(i)-dlon_amie/2 and $
                         lons lt lons_amie(i)+dlon_amie/2, countl)

            if (countl gt 0) then begin

                eflux  = mean(final_eflux(locl,iLat))
                efluxe = stddev(final_eflux(locl,iLat))
                eave   = mean(final_echar(locl,iLat))
                eavee  = stddev(final_echar(locl,iLat))
                ped    = 40.0 * eave / (16.0+eave^2) * sqrt(eflux)
                hal    = 0.45 * eave^0.85 * ped

;                eflux  = final_eflux(locl(0),iLat)
;                efluxe = eflux/2.0
;                eave   = final_echar(locl(0),iLat)
;                eavee  = eave/2.0
;                ped    = 40.0 * eave / (16.0+eave^2) * sqrt(eflux)
;                hal    = 0.45 * eave^0.85 * ped

                if (eflux gt 0.01 and hal gt 0.01) then $
                  printf,1, $
                  format="(i4,5i3,f6.2,f7.2,f6.2,i5,5f7.2,f6.2,4f7.2,1x,1a)",$
                  itime, lats(iLat), 0.0, lons_amie(i), $
                  countl, 0.0, $
                  ped, ped/2, hal, hal/2, 0.0, $
                  eflux, efluxe+eflux/2.0, eave, eavee+eave/2.0, " "

            endif

        endfor

    endfor
                  
endfor

close,1

end


