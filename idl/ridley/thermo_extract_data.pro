
filelist = findfile("-t *.bin")
nFiles = n_elements(filelist)

iVars = [0,1,2]
gitm_read_bin_1var,filelist(0), data, time, nVars, Vars, version, VarsToGet = iVars
lons = reform(data(0,*,0,0))/!dtor
lats = reform(data(1,0,*,0))/!dtor
alts = reform(data(2,0,0,*))/1000.0
nLons = n_elements(lons)
nLats = n_elements(lats)

display, Vars
if (n_elements(iVar) eq 0) then iVar = 3
iVar = fix(ask('Variable to extract',tostr(iVar)))

display, alts
if (n_elements(iAlt) eq 0) then iAlt = 20
iAlt = fix(ask('Altitude to extract',tostr(iAlt)))

alt = alts(iAlt)

iVars = [iVar]

for iFile=0,nFiles-1 do begin

   infile = filelist(iFile)
   print,'Reading '+infile
   
   gitm_read_bin_1var, infile, data, time, nVars, Vars, version, $
                       VarsToGet = iVars

   c_r_to_a, itime, time
   c_a_to_ymdhms, itime, sdate
   sdate = strmid(sdate,0,8)+'_'+strmid(sdate,8,6)
   file = 'gitm_'+tostr(iVar,2)+'_'+sdate+'.txt'

   print,'Writing '+file
   openw,2,file

   printf,2,''
   printf,2,'GITM output at a single altitude'
   printf,2,''
   printf,2,'Variable'
   printf,2,Vars(iVar)
   printf,2,''
   printf,2,'Altitude (km)'
   printf,2,alt
   printf,2,''
   printf,2,'Longitudes (deg) - outer loop'
   printf,2,nLons
   printf,2,lons
   printf,2,''
   printf,2,'Latitudes (deg) - inner loop'
   printf,2,nLats
   printf,2,lats
   printf,2,''
   printf,2,'Time'
   printf,2,iTime
   printf,2,''
   printf,2,'Start'
   format = "("+tostr(nLats)+"e11.3)"
   for iLon = 0, nLons-1 do $
      printf,2,data(0,iLon,*,iAlt), format = format
   close,2

endfor

end

