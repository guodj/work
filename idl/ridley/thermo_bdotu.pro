
VarsToGet=[0,1,2,3,4,5,6,7]
magfile = '/raid3/idl/extras/3DMAG_t110201_000000.bin'
gitm_read_bin_1var, magfile, magdata, time, nVars, Vars, version, $
                        VarsToGet = VarsToGet

lons = reform(magdata(0,*,*,0))/!dtor
lats = reform(magdata(1,*,*,0))/!dtor
alts = reform(magdata(2,0,0,*))/1000.0

filelist = findfile('3DALL*.bin')
nFiles = n_elements(filelist)

VarsToGet = [2]
file = filelist(0)
gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                    VarsToGet = VarsToGet
alts = data(0,0,0,*)/1000.0
l = where(alts gt 300)
iAlt = l(0)

mLats = reform(magdata(3,*,*,iAlt))
mLons = reform(magdata(4,*,*,iAlt))

bx = reform(magdata(5,*,*,iAlt))
by = reform(magdata(6,*,*,iAlt))
bz = reform(magdata(7,*,*,iAlt))
bm = sqrt(bx^2 + by^2 + bz^2)
bx = bx/bm
by = by/bm
bz = bz/bm

VarsToGet = [16,17,18, 36,37,38]

locNorth = where(lats gt 50.0 and lats lt 90.0 and $
                 lons gt 0.0 and lons lt 360.0)
locSouth = where(lats lt -50.0 and lats gt -90.0 and $
                 lons gt 0.0 and lons lt 360.0)

locMagNorth = where(mLats gt 50.0 and lats lt 90.0 and $
                 lons gt 0.0 and lons lt 360.0)
locMagSouth = where(mLats lt -50.0 and lats gt -90.0 and $
                 lons gt 0.0 and lons lt 360.0)

north_bdun = fltarr(nFiles)
south_bdun = fltarr(nFiles)
north_bdui = fltarr(nFiles)
south_bdui = fltarr(nFiles)

north_day_bdun = fltarr(nFiles)
south_day_bdun = fltarr(nFiles)
north_day_bdui = fltarr(nFiles)
south_day_bdui = fltarr(nFiles)

north_mag_sza = fltarr(nFiles)
south_mag_sza = fltarr(nFiles)
north_mag_bdun = fltarr(nFiles)
south_mag_bdun = fltarr(nFiles)
north_mag_bdui = fltarr(nFiles)
south_mag_bdui = fltarr(nFiles)

alltimes = dblarr(nFiles)

for iFile=0,nFiles-1 do begin

   file = filelist(iFile)
   print, 'reading file : ',file
   gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                       VarsToGet = VarsToGet
   alltimes(iFile) = time
   unx = reform(data(0,*,*,iAlt))
   uny = reform(data(1,*,*,iAlt))
   unz = reform(data(2,*,*,iAlt))
   uix = reform(data(3,*,*,iAlt))
   uiy = reform(data(4,*,*,iAlt))
   uiz = reform(data(5,*,*,iAlt))

   bdun = (bx*unx + by*uny + bz*unz)*sign(bz)
   bdui = (bx*uix + by*uiy + bz*uiz)*sign(bz)

   thermo_sza, lats, lons, time, LocalTime, sza

   locNorthDay = where(lats gt 40.0 and lats lt 80.0 and $
                    LocalTime gt 10 and LocalTime lt 14 and $
                    lons gt 0.0 and lons lt 360.0)
   locSouthDay = where(lats lt -40.0 and lats gt -80.0 and $
                    LocalTime gt 10 and LocalTime lt 14 and $
                    lons gt 0.0 and lons lt 360.0)

   north_mag_bdun(iFile) = mean(bdun(locMagNorth))
   south_mag_bdun(iFile) = mean(bdun(locMagSouth))
   north_mag_bdui(iFile) = mean(bdui(locMagNorth))
   south_mag_bdui(iFile) = mean(bdui(locMagSouth))
   north_bdun(iFile) = mean(bdun(locNorth))
   south_bdun(iFile) = mean(bdun(locSouth))
   north_bdui(iFile) = mean(bdui(locNorth))
   south_bdui(iFile) = mean(bdui(locSouth))

   north_day_bdun(iFile) = mean(bdun(locNorthDay))
   south_day_bdun(iFile) = mean(bdun(locSouthDay))
   north_day_bdui(iFile) = mean(bdui(locNorthDay))
   south_day_bdui(iFile) = mean(bdui(locSouthDay))
   
   north_mag_sza(iFile) = mean(sza(locMagNorth))
   south_mag_sza(iFile) = mean(sza(locMagSouth))

endfor

midtime = mean(alltimes)
c_r_to_a, itime, midtime
c_a_to_ymd, itime, ymd
save, file = 'mag_winds_'+ymd+'.save', $
      north_mag_bdun, south_mag_bdun, north_mag_bdui, south_mag_bdui, $
      north_bdun, south_bdun, north_bdui, south_bdui, $
      north_day_bdun, south_day_bdun, north_day_bdui, south_day_bdui, $
      north_mag_sza, south_mag_sza, alltimes

end
